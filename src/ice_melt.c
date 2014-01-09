/*
 * SUMMARY:      IceMelt.c - Calculate snow accumulation and melt for the lake model
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Mon Apr 21 17:07:05 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate snow accumulation and melt using an energy balance
 *               approach for a two layer snow model
 * DESCRIP-END.
 * FUNCTIONS:    SnowMelt()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: IceMelt()

  Purpose      : Calculate snow accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    double delta_t               - Model timestep (hours)
    double z2           - Reference height (m) 
    double displacement          - Displacement height (m)
    double aero_resist           - Aerodynamic resistance (uncorrected for
                                   stability) (s/m)
    double *aero_resist_used     - Aerodynamic resistance (corrected for
                                   stability) (s/m)
    double atmos->density        - Density of air (kg/m3)
    double atmos->vp             - Actual vapor pressure of air (Pa) 
    double Le           - Latent heat of vaporization (J/kg3)
    double atmos->net_short      - Net exchange of shortwave radiation (W/m2)
    double atmos->longwave       - Incoming long wave radiation (W/m2)
    double atmos->pressure       - Air pressure (Pa)
    double RainFall              - Amount of rain (m)
    double Snowfall              - Amount of snow (m)
    double atmos->air_temp       - Air temperature (C)
    double atmos->vpd            - Vapor pressure deficit (Pa)
    double wind                  - Wind speed (m/s)
    double snow->pack_water      - Liquid water content of snow pack 
    double snow->surf_water	 - Liquid water content of surface layer 
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Total mass flux of water vapor to or from
                                   snow (m/time step)
    double snow->blowing_flux;   - Mass flux of water vapor to or from
                                   blowing snow (m/time step)
    double snow->surface_flux;   - Mass flux of water vapor to or from
                                   snow pack (m/time step)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Modifies     :
    double atmos->melt           - Amount of snowpack outflow (m)
    double snow->pack_water      - Liquid water content of snow pack 
    double snow->surf_water	 - Liquid water content of surface layer 
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Total mass flux of water vapor to or from
                                   snow (m/timestep)
    double snow->surface_flux;   - Mass flux of water vapor to or from
                                   snow pack (m/timestep)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Comments     :

  Modifications:
  11-18-02 Modified method by which lake coverage fraction and ice height 
           are updated.                                                		LCB
  04-Jun-04 Added descriptive error message to beginning of screen dump in
	    ErrorPrintIcePackEnergyBalance.					TJB
  16-Jul-04 Changed VaporMassFlux to vapor_flux, to make it consistent with
	    IceEnergyBalance(), in which VaporMassFlux is in (kg/m2s) and
	    vapor_flux is in (m/timestep).  Changed calculations involving
	    vapor_flux to reflect these new units.				TJB
  25-Aug-04 Added calculations for surface_flux and blowing_flux.  Note that
	    blowing_flux is currently set to 0.  This can be replaced in the
	    future with a call to CalcBlowingSnow().				TJB
  27-Aug-04 Replaced *DeltaColdContent with DeltaColdContent in parameter
	    list for ErrorIcePackEnergyBalance() and
	    ErrorPrintIcePackEnergyBalance() so that types match in the caller
	    and callee.								TJB
  21-Sep-04 Added ErrorString to store error messages from
	    root_brent.								TJB
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance used
	    in flux calculations.						TJB
  04-Oct-04 Merged with Laura Bowling's updated lake model code.  Now
	    sublimation from blowing snow is calculated for lakes.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-03 Return ERROR value on error from CalcBlowingSnow and root_brent.	GCT/KAC. 
  2007-Aug-31 Checked root_brent return value against -998 rather than -9998.	JCA
  2007-Sep-25 Return ERROR instead of exiting, if ice_melt could not converge to
	      a solution in root_brent.						JCA
  2007-Nov-06 New parameters for CalcBlowingSnow().  Replaced lake.fraci,
	      lake.hice with lake.areai and lake.ice_water_eq.  More accurate
	      accounting of lake_snow->surf_water.				LCB via TJB
  2008-Jan-23 Modified lake snow pack to have 2 layers, similar to the modeling
	      of snow pack on uplands.						LCB via TJB
  2008-Sep-09 Reduced the fetch used with blowing snow calculations over lakes
	      from 2000 m to 100 m.						LCB via TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.			TJB
  2009-Dec-11 Replaced "assert" statements with "if" statements.		TJB
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.								TJB
  2013-Dec-27 Moved SPATIAL_SNOW from compile-time to run-time options.	TJB
*****************************************************************************/
int ice_melt(double            z2,
	      double            aero_resist,
	      double            *aero_resist_used,
	      double            Le,
	      snow_data_struct *snow,
	      lake_var_struct  *lake,
	      int               delta_t,
	      double            displacement,
	      double            Z0,
	      double            surf_atten,
	      double            rainfall,
	      double            snowfall,
	      double            wind,
	      double            Tcutoff,
	      double            air_temp,
	      double            net_short,
	      double            longwave,
	      double            density,
	      double            pressure,
	      double            vpd,
	      double            vp,
	      double           *melt,
	      double           *save_advection,
	      double           *save_deltaCC,
	      double           *save_SnowFlux,
	      double           *save_latent,
	      double           *save_sensible,
	      double           *save_Qnet,
	      double           *save_refreeze_energy,
	      double           *save_LWnet,
	      double            fracprv)
{
  extern option_struct   options;

  int    Twidth;

  double DeltaPackCC;            /* Change in cold content of the pack */
  double DeltaPackSwq;           /* Change in snow water equivalent of the pack (m) */
  double InitialSwq;             /* Initial snow water equivalent (m) */
  double InitialIce;
  double MassBalanceError;       /* Mass balance error (m) */
  double MaxLiquidWater;         /* Maximum liquid water content of pack (m) */
  double OldTSurf;               /* Old snow surface temperature (C) */
  double Qnet;                   /* Net energy exchange at the surface (W/m2) */
  double PackRefreezeEnergy;     /* refreeze/melt energy in pack layer (W/m2) */
  double RefreezeEnergy;         /* refreeze energy (W/m2) */
  double RefrozenWater;          /* Amount of refrozen water (m) */
  double SnowFallCC;             /* Cold content of new snowfall (J) */
  double SurfaceCC;
  double PackCC;
  double SurfaceSwq;
  double PackSwq;
  double PackIce;
  double SnowMelt;               /* Amount of snow melt during time interval (m water equivalent) */
  double IceMelt;
  double LWnet;
  double avgcond;
  double SWconducted;
  double SnowIce;
  double LakeIce;
  double Ice;
  double SnowFall;
  double RainFall;
  double vapor_flux;
  double blowing_flux;
  double surface_flux;
  double advection;
  double deltaCC;
  double SnowFlux;		/* thermal flux through snowpack from ground */
  double latent_heat;		
  double latent_heat_sub;		
  double sensible_heat;		
  double Ls;
  double melt_energy = 0.;

  char ErrorString[MAXSTRING];

  SnowFall = snowfall / 1000.; /* convert to m */
  RainFall = rainfall / 1000.; /* convert to m */
  IceMelt = 0.0;
  RefrozenWater = 0.0;
  
  InitialSwq = snow->swq;
  OldTSurf = snow->surf_temp;

  /* Initialize snowpack variables */
  SnowIce  = snow->swq - snow->pack_water - snow->surf_water;
  LakeIce = lake->ice_water_eq / lake->areai;        /* meters of water equivalent based on average ice thickness. */
  InitialIce = LakeIce;
  Ice = SnowIce + LakeIce;

  /* Reconstruct snow pack */
  if (Ice > MAX_SURFACE_SWE)
    SurfaceSwq = MAX_SURFACE_SWE;
  else
    SurfaceSwq = Ice;
  if(SurfaceSwq <= SnowIce) {
    PackSwq = SnowIce - SurfaceSwq;
    PackIce = LakeIce;
  }
  else {
    PackSwq = 0.;
    PackIce = Ice - SurfaceSwq;
  }

  /* Calculate cold contents */
  SurfaceCC = CH_ICE * SurfaceSwq * snow->surf_temp;
  PackCC = CH_ICE * (PackSwq+PackIce) * snow->pack_temp;
  if (air_temp > 0.0)
    SnowFallCC = 0.0;
  else
    SnowFallCC = CH_ICE * SnowFall * air_temp;

  /* Distribute fresh snowfall */
  /* Surface layer was not already full, snow will exceed space.  What happens to snow if SurfaceSwq = MAX_SURFACE_SWE?*/
  //  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq) && (MAX_SURFACE_SWE - SurfaceSwq) > SMALL) {
  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq)) {
    DeltaPackSwq = SurfaceSwq + SnowFall - MAX_SURFACE_SWE;
    if (DeltaPackSwq > SurfaceSwq)
      DeltaPackCC = SurfaceCC + (SnowFall - MAX_SURFACE_SWE)/SnowFall * SnowFallCC;
    else
      DeltaPackCC = DeltaPackSwq/SurfaceSwq * SurfaceCC;
    SurfaceSwq = MAX_SURFACE_SWE;
    SurfaceCC += SnowFallCC - DeltaPackCC;
    PackSwq += DeltaPackSwq;
    PackCC += DeltaPackCC;
  }
  else {
    SurfaceSwq += SnowFall;
    SurfaceCC += SnowFallCC;
    DeltaPackCC = 0;
  }
  if (SurfaceSwq > 0.0)
    snow->surf_temp = SurfaceCC/(CH_ICE * SurfaceSwq);
  else
    snow->surf_temp = 0.0;
  if (PackSwq+PackIce > 0.0)
    snow->pack_temp = PackCC/(CH_ICE * (PackSwq+PackIce));
  else
    snow->pack_temp = 0.0;

  /* Adjust ice and snow->surf_water */
  SnowIce += SnowFall;
  Ice += SnowFall;
  snow->surf_water += RainFall;
  
  icerad (net_short, lake->hice, SnowIce*RHO_W/RHOSNOW, &avgcond, &SWconducted, &deltaCC);

  /* Calculate blowing snow sublimation (m/timestep) */
  // Currently have hard-wired parameters that are approximated for ice-covered area:
  // lag-one autocorrelation = 0.95, sigma_slope = .005 (both appropriate for
  // flat terrain. Fetch = 2000 m (i.e. unlimited fetch), roughness and displacement
  // calculated assuming 10 cm high protrusions on frozen ponds.

  if(options.BLOWING && snow->swq > 0.) {
    Ls = (677. - 0.07 * snow->surf_temp) * JOULESPCAL * GRAMSPKG;
    snow->blowing_flux = CalcBlowingSnow((double) delta_t, air_temp,
					 snow->last_snow, snow->surf_water,
					 wind, Ls, density,
					 pressure, vp, Z0,
					 z2, snow->depth, .95, 0.005,
					 snow->surf_temp, 0, 1, 100.,
					 .067, .0123, &snow->transport);
    if ( (int)snow->blowing_flux == ERROR ) {
      fprintf( stderr, "ERROR: ice_melt.c has an error from the call to CalcBlowingSnow\n");
      fprintf( stderr, "Exiting module\n" );
      return ( ERROR );
    }

    snow->blowing_flux *= delta_t*SECPHOUR/RHO_W;
  }
  else
    snow->blowing_flux = 0.0;

  /* Store sublimation terms in temporary variables */
  vapor_flux = snow->vapor_flux;
  blowing_flux = snow->blowing_flux;
  surface_flux = snow->surface_flux;

  /* Calculate the surface energy balance for snow_temp = 0.0 */

  Qnet = CalcIcePackEnergyBalance((double)0.0, (double)delta_t, aero_resist,
				  aero_resist_used, z2, displacement, Z0, wind, net_short, 
				  longwave, density, Le, air_temp,
				  pressure * 1000., vpd * 1000., vp * 1000.,
				  RainFall, SurfaceSwq, 
				  snow->surf_water, OldTSurf, &RefreezeEnergy,
				  &vapor_flux, &blowing_flux, &surface_flux,
				  &advection, deltaCC, Tcutoff, avgcond, 
				  SWconducted, snow->swq*RHO_W/RHOSNOW, 
				  RHOSNOW, surf_atten, &SnowFlux, 
				  &latent_heat, &latent_heat_sub, &sensible_heat, &LWnet);

  snow->vapor_flux = vapor_flux;
  snow->surface_flux = surface_flux;
  save_refreeze_energy[0] = RefreezeEnergy;

  /* Check that snow swq exceeds minimum value for model stability */
  //if ( !UNSTABLE_SNOW ) {

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  if (Qnet == 0.0) {
    snow->surf_temp = 0.0;
    if (RefreezeEnergy >= 0.0) {                     /* Surface is freezing. */
      RefrozenWater = RefreezeEnergy/(Lf * RHO_W) * delta_t * SECPHOUR; 
      if (RefrozenWater > snow->surf_water) {
        RefrozenWater = snow->surf_water;
        RefreezeEnergy = RefrozenWater * Lf * RHO_W/(delta_t * SECPHOUR);
      } 
      melt_energy  += RefreezeEnergy;
      SurfaceSwq   += RefrozenWater;
      SnowIce      += RefrozenWater;
      Ice          += RefrozenWater;
      snow->surf_water   -= RefrozenWater;
      if (snow->surf_water < 0.0) snow->surf_water = 0.0;
      SnowMelt      = 0.0;
    }
    else {

      /* Calculate snow melt if refreeze energy is negative */      
      SnowMelt = fabs(RefreezeEnergy)/(Lf * RHO_W) * delta_t * SECPHOUR;
      melt_energy += RefreezeEnergy;
    }

    /* Adjust snow->surf_water for vapor_flux */
    if (snow->surf_water < -(snow->vapor_flux)) {
      // if vapor_flux exceeds stored water, we not only need to
      // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
      snow->blowing_flux *= -(snow->surf_water) / snow->vapor_flux;
      snow->vapor_flux    = -(snow->surf_water);
      snow->surface_flux  = -(snow->surf_water) - snow->blowing_flux;
      snow->surf_water    = 0.0;
    }  
    else {
      snow->surf_water   += snow->vapor_flux;
    }
 
    /* If SnowMelt < Ice, there was incomplete melting of the snow/ice */
    if (SnowMelt < Ice) {
      /* Subtract melt from snow pack and lake ice */
      /* Since we redefine the snow surface layer to always encompass the topmost
         portion of the snow, we essentially reduce the bottom snow layer first.
         This is only for accounting purposes and may not reflect where in the
         pack the melt actually occurs. */
      if (SnowMelt <= PackSwq) { /* Only melt part of the pack (bottom) layer. */
        snow->surf_water += SnowMelt;
        PackSwq          -= SnowMelt;
        Ice              -= SnowMelt;
        SnowIce          -= SnowMelt;
      }
      else if (SnowMelt <= SnowIce) { /* Melt all of pack layer and part of surface layer. */
        snow->surf_water += SnowMelt + snow->pack_water;
        snow->pack_water  = 0.0;
        SurfaceSwq       -= (SnowMelt - PackSwq);
        PackSwq           = 0.0;
        SnowIce          -= SnowMelt;
        Ice              -= SnowMelt;
      }
      else  { /* Melt snow pack completely and also part of the ice */
        snow->surf_water += SnowIce + snow->pack_water;
        snow->pack_water  = 0.0;
        PackSwq           = 0.0;
        Ice              -= SnowMelt;
        LakeIce          -= SnowMelt-SnowIce;
        IceMelt           = SnowMelt-SnowIce;
        if(SurfaceSwq > SnowMelt) {
          SurfaceSwq -= SnowMelt;
        }
        else {
          SurfaceSwq = 0.0;
          PackIce -= (SnowMelt-SurfaceSwq-PackSwq);
        }
        SnowIce = 0.0;
      }
    }
    /* Else, SnowMelt > TotalIce and there was complete melting of the snow and ice */
    else {
      snow->surf_water += SnowIce + snow->pack_water;
      snow->pack_water  = 0.0;
      PackSwq           = 0.0;
      SurfaceSwq        = 0.0;
      SnowIce           = 0.0;
      SnowMelt          = Ice;
      IceMelt           = LakeIce;
      LakeIce           = 0.0;
      PackIce           = 0.0;
      Ice               = 0.0;
      snow->surf_temp   = 0.0;
      snow->pack_temp   = 0.0;
      /* readjust melt energy to account for melt only of available snow */
      melt_energy -= RefreezeEnergy;
      RefreezeEnergy = RefreezeEnergy / fabs(RefreezeEnergy) * SnowMelt * Lf * RHO_W / ( delta_t );
      melt_energy += RefreezeEnergy;
    } 
  }
  
  /* Else, IceEnergyBalance(T=0.0) <= 0.0 */
  else  {
    /* Calculate surface layer temperature using "Brent method" */
    if (SurfaceSwq > MIN_SWQ_EB_THRES) {
      snow->surf_temp = root_brent((double)(snow->surf_temp-SNOW_DT), 
				   (double)(snow->surf_temp+SNOW_DT), ErrorString,
				   IceEnergyBalance, (double)delta_t, 
				   aero_resist, aero_resist_used, z2, 
				   displacement, Z0, wind, net_short, longwave,
				   density, Le, air_temp, pressure * 1000.,
				   vpd * 1000., vp * 1000., RainFall, 
				   SurfaceSwq,
				   snow->surf_water, OldTSurf, &RefreezeEnergy, 
				   &vapor_flux, &blowing_flux, &surface_flux,
				   &advection, deltaCC, Tcutoff, 
				   avgcond, SWconducted,
				   snow->swq*RHO_W/RHOSNOW,
				   RHOSNOW,surf_atten,&SnowFlux,
				   &latent_heat, &latent_heat_sub, &sensible_heat, &LWnet);

      if (snow->surf_temp <= -998) {
        if (options.TFALLBACK) {
          snow->surf_temp = OldTSurf;
          snow->surf_temp_fbflag = 1;
          snow->surf_temp_fbcount++;
        }
        else {
          ErrorIcePackEnergyBalance(snow->surf_temp, (double)delta_t, aero_resist,
				    aero_resist_used, z2, displacement, Z0, wind, net_short,
				    longwave, density, Le, air_temp,
				    pressure * 1000., vpd * 1000., vp * 1000.,
				    RainFall, SurfaceSwq, 
				    snow->surf_water, OldTSurf, &RefreezeEnergy,
				    &vapor_flux, &blowing_flux, &surface_flux,
				    &advection, deltaCC, Tcutoff, 
				    avgcond, SWconducted,		
				    snow->swq*RHO_W/RHOSNOW, RHOSNOW,surf_atten,
				    &SnowFlux, &latent_heat, &latent_heat_sub,
				    &sensible_heat, &LWnet, ErrorString);
          return( ERROR );
        }
      }
    }
    else {
//        fprintf(stderr,"Snow/Ice layer is too thin to solve separately \n");
      snow->surf_temp = 999;
    }
    if (snow->surf_temp > -998 && snow->surf_temp < 999) {
      Qnet = CalcIcePackEnergyBalance(snow->surf_temp, (double)delta_t, aero_resist,
				      aero_resist_used, z2, displacement, Z0, wind, net_short,
				      longwave, density, Le, air_temp,
				      pressure * 1000., vpd * 1000., 
				      vp * 1000.,RainFall, SurfaceSwq, 
				      snow->surf_water, OldTSurf, &RefreezeEnergy, 
				      &vapor_flux, &blowing_flux, &surface_flux,
				      &advection, deltaCC, Tcutoff, 
				      avgcond, SWconducted, snow->swq*RHO_W/RHOSNOW, 
				      RHOSNOW,surf_atten, &SnowFlux, &latent_heat,
				      &latent_heat_sub, &sensible_heat, &LWnet);

      snow->vapor_flux = vapor_flux;
      snow->surface_flux = surface_flux;
      save_refreeze_energy[0] = RefreezeEnergy;

      /* since we iterated, the surface layer is below freezing and no snowmelt */ 

      SnowMelt = 0.0;
      IceMelt  = 0.0;

      /* Since updated snow_temp < 0.0, all of the liquid water in the surface
         layer has been frozen */ 

      SnowIce         += snow->surf_water;
      Ice             += snow->surf_water;
      melt_energy     += snow->surf_water * Lf * RHO_W/(delta_t * SECPHOUR);
      RefrozenWater    = snow->surf_water;
      snow->surf_water = 0.0;

      /* Adjust SurfaceSwq for vapor_flux */
      if (SurfaceSwq < -(snow->vapor_flux)) {
        // if vapor_flux exceeds stored snow/ice, we not only need to
        // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
        if(SurfaceSwq > SnowIce) {
          snow->blowing_flux *= -(SurfaceSwq) / snow->vapor_flux;
          snow->vapor_flux    = -SurfaceSwq;
          snow->surface_flux  = -SurfaceSwq - snow->blowing_flux;
          LakeIce            -= SurfaceSwq - SnowIce;
          Ice                 = PackIce;
          SnowIce             = 0.0;
        }
        else {
          snow->blowing_flux *= -(SurfaceSwq) / snow->vapor_flux;
          snow->vapor_flux    = -SurfaceSwq;
          snow->surface_flux  = -SurfaceSwq - snow->blowing_flux;
          SurfaceSwq          = 0.0;
          Ice                 = PackSwq + PackIce;
        }
      }
      else {
        SurfaceSwq += snow->vapor_flux;
        if(SnowIce > -(snow->vapor_flux))
          SnowIce += snow->vapor_flux;
        else {
          LakeIce += (snow->vapor_flux + SnowIce);
          SnowIce = 0.;
        }
        Ice += snow->vapor_flux;
      }

    }

    else {
      snow->surf_temp = 999;
    }

  }

  /* Done with iteration etc, now Update the liquid water content of the
     surface layer */ 
  
  if(SnowIce > SurfaceSwq)
    MaxLiquidWater = LIQUID_WATER_CAPACITY * SurfaceSwq;
  else
    MaxLiquidWater = LIQUID_WATER_CAPACITY * SnowIce;
  if (snow->surf_water > MaxLiquidWater) {
    melt[0]          = snow->surf_water - MaxLiquidWater;
    snow->surf_water = MaxLiquidWater;
  }
  else
    melt[0] = 0.0;

  /* Refreeze liquid water in the pack.
     variable 'RefreezeEnergy' is the heat released to the snow pack
     if all liquid water were refrozen.
     if RefreezeEnergy < PackCC then all water IS refrozen
     PackCC always <=0.0
     WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does
     not involve energy transported to the pixel.  Instead heat from the snow
     pack is used to refreeze water */
  snow->pack_water += melt[0]; /* add surface layer outflow to pack liquid water*/
  PackRefreezeEnergy = snow->pack_water * Lf * RHO_W;

  /* calculate energy released to freeze*/
  if (PackCC < -PackRefreezeEnergy) { /* cold content not fully depleted*/
    PackSwq         += snow->pack_water;    /* refreeze all water and update*/
    Ice             += snow->pack_water;
    SnowIce         += snow->pack_water;
    snow->pack_water = 0.0;
    if (PackSwq + PackIce > 0.0) {
      PackCC          = (PackSwq+PackIce) * CH_ICE * snow->pack_temp + PackRefreezeEnergy;
      snow->pack_temp = PackCC / (CH_ICE * (PackSwq+PackIce));
      if(snow->pack_temp > 0.) snow->pack_temp = 0.;
    }
    else
      snow->pack_temp = 0.0;
  }
  else {
    /* cold content has been either exactly satisfied or exceeded. If
       PackCC = refreeze then pack is ripe and all pack water is
       refrozen, else if energy released in refreezing exceeds PackCC
       then exactly the right amount of water is refrozen to satify PackCC.
       The refrozen water is added to PackSwq and Ice */
    snow->pack_temp   = 0.0;
    DeltaPackSwq      = -PackCC/(Lf * RHO_W);
    snow->pack_water -= DeltaPackSwq;
    PackSwq          += DeltaPackSwq;
    Ice              += DeltaPackSwq;
    SnowIce          += DeltaPackSwq;
  }

  /* Update the liquid water content of the pack */
  MaxLiquidWater = LIQUID_WATER_CAPACITY * PackSwq;
  if (snow->pack_water > MaxLiquidWater) {
    melt[0]          = snow->pack_water - MaxLiquidWater;
    snow->pack_water = MaxLiquidWater;
  }
  else
    melt[0] = 0.0;

  /* Update snow properties */
  Ice  = PackIce + PackSwq + SurfaceSwq;
  if (Ice > MAX_SURFACE_SWE) {
    SurfaceCC   = CH_ICE * snow->surf_temp * SurfaceSwq;
    PackCC      = CH_ICE * snow->pack_temp * (PackSwq + PackIce);
    if (SurfaceSwq > MAX_SURFACE_SWE) {
      PackCC     += SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      SurfaceCC  -= SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq;
      PackSwq    += SurfaceSwq - MAX_SURFACE_SWE;
      SurfaceSwq -= SurfaceSwq - MAX_SURFACE_SWE;
    }
    else if ( SurfaceSwq < MAX_SURFACE_SWE) {
      PackCC     -= PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / (PackSwq+PackIce);
      SurfaceCC  += PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / (PackSwq+PackIce);
      PackSwq    -= MAX_SURFACE_SWE - SurfaceSwq;
      SurfaceSwq += MAX_SURFACE_SWE - SurfaceSwq;
    }
    snow->pack_temp      = PackCC / (CH_ICE * (PackSwq+PackIce));
    snow->surf_temp      = SurfaceCC / (CH_ICE * SurfaceSwq);
  }
  else {
    PackSwq = 0.0;
    PackCC  = 0.0;
    PackIce = 0.0;
    snow->pack_temp  = 0.0;
  }
  snow->swq = SnowIce + snow->surf_water + snow->pack_water;
  lake->ice_water_eq = LakeIce * lake->areai;
  lake->volume -= (InitialIce - LakeIce - IceMelt) * lake->areai;
  if (lake->ice_water_eq <= 0.0) {
    lake->ice_water_eq = 0.0;
  }
  if (options.SPATIAL_SNOW) {
  /*  snow->coverage = calc_snow_coverage(&snow->store_snow,
    soil_con->max_snow_distrib_slope,
    old_coverage, snow->swq,
    old_swq, snow->depth, old_depth,
    melt + snow->vapor_flux,
    &snow->max_snow_depth, snowfall,
    &snow->store_swq,
    &snow->snow_distrib_slope,
    &snow->store_coverage);
    */
  }
  else {
    if ( snow->swq > 0 ) snow->coverage = 1.;
    else snow->coverage = 0.;
  }

  /* Mass balance test */
  MassBalanceError = (InitialSwq - snow->swq) + (InitialIce - LakeIce)
                     + (RainFall + SnowFall) - IceMelt - melt[0] + snow->vapor_flux;
  //   if(fabs(MassBalanceError) > SMALL)
  //  fprintf(stderr, "MassBalanceError = %g %e %e %e, %e, \n", MassBalanceError, InitialIce - LakeIce, IceMelt, InitialSwq, snow->swq);

  melt[0] *= 1000.; /* converts back to mm */
  snow->mass_error = MassBalanceError;
  snow->coldcontent        = SurfaceCC;
  snow->vapor_flux *= -1.;
  *save_LWnet = LWnet;
  *save_advection = advection;
  *save_deltaCC = deltaCC;
  *save_SnowFlux = SnowFlux;
  *save_latent = latent_heat + latent_heat_sub;
  *save_sensible = sensible_heat;
  *save_refreeze_energy    = RefreezeEnergy;
  *save_Qnet = Qnet;

  return (0);

}

/*****************************************************************************
  Function name: CalcIcePackEnergyBalance()

  Purpose      : Dummy function to make a direct call to
                 IceEnergyBalance() possible.

  Required     : 
    double TSurf - IcePack surface temperature (C)
    other arguments required by IcePackEnergyBalance()

  Returns      :
    double Qnet - Net energy exchange at the IcePack snow surface (W/m^2)

  Modifies     : none

  Comments     : function is local to this module
*****************************************************************************/
double CalcIcePackEnergyBalance(double Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  double Qnet;                   /* Net energy exchange at the IcePack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = IceEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

double ErrorIcePackEnergyBalance(double Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  double Qnet;                   /* Net energy exchange at the IcePack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = ErrorPrintIcePackEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

double ErrorPrintIcePackEnergyBalance(double TSurf, va_list ap)
{


  double Dt;                     /* Model time step (hours) */
  double Ra;                     /* Aerodynamic resistance (s/m) */
  double *Ra_used;               /* Aerodynamic resistance (s/m) after stability correction */
  double Z;                      /* Reference height (m) */
  double Displacement;           /* Displacement height (m) */
  double Z0;                     /* surface roughness height (m) */
  double Wind;                   /* Wind speed (m/s) */
  double ShortRad;               /* Net incident shortwave radiation (W/m2) */
  double LongRadIn;              /* Incoming longwave radiation (W/m2) */
  double AirDens;                /* Density of air (kg/m3) */
  double Lv;                     /* Latent heat of vaporization (J/kg3) */
  double Tair;                   /* Air temperature (C) */
  double Press;                  /* Air pressure (Pa) */
  double Vpd;			/* Vapor pressure deficit (Pa) */
  double EactAir;                /* Actual vapor pressure of air (Pa) */
  double Rain;                   /* Rain fall (m/timestep) */
  double SweSurfaceLayer;        /* Snow water equivalent in surface layer (m)
                                 */ 
  double SurfaceLiquidWater;     /* Liquid water in the surface layer (m) */
  double OldTSurf;               /* Surface temperature during previous time
                                   step */ 
  double *RefreezeEnergy;        /* Refreeze energy (W/m2) */
  double *vapor_flux;            /* Total mass flux of water vapor to or from
                                    snow (m/timestep) */
  double *blowing_flux;          /* Mass flux of water vapor to or from
                                    blowing snow (m/timestep) */
  double *surface_flux;          /* Mass flux of water vapor to or from
                                    snow pack (m/timestep) */
  double *AdvectedEnergy;        /* Energy advected by precipitation (W/m2) */
  double DeltaColdContent;       /* Change in cold content (W/m2) */
  double Tfreeze;
  double AvgCond;
  double SWconducted;
  double SWabsorbed;
  double SnowDepth;  
  double SnowDensity; 
  double SurfAttenuation; 
  
  /* end of list of arguments in variable argument list */
  double *GroundFlux;
  double *LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;	/* Latent heat exchange at surface (W/m2) due to sublimation */
  double *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */
  double *LWnet;

  char *ErrorString;

  /* initialize variables */
  Dt                 = (double) va_arg(ap, double);
  Ra                 = (double) va_arg(ap, double);
  Ra_used            = (double *) va_arg(ap, double *);
  Z                  = (double) va_arg(ap, double);
  Displacement       = (double) va_arg(ap, double);
  Z0                 = (double) va_arg(ap, double);
  Wind               = (double) va_arg(ap, double);
  ShortRad           = (double) va_arg(ap, double);
  LongRadIn          = (double) va_arg(ap, double);
  AirDens            = (double) va_arg(ap, double);
  Lv                 = (double) va_arg(ap, double);
  Tair               = (double) va_arg(ap, double);
  Press              = (double) va_arg(ap, double);
  Vpd                = (double) va_arg(ap, double);
  EactAir            = (double) va_arg(ap, double);
  Rain               = (double) va_arg(ap, double);
  SweSurfaceLayer    = (double) va_arg(ap, double);
  SurfaceLiquidWater = (double) va_arg(ap, double);
  OldTSurf           = (double) va_arg(ap, double);
  RefreezeEnergy     = (double *) va_arg(ap, double *);
  vapor_flux         = (double *) va_arg(ap, double *);
  blowing_flux       = (double *) va_arg(ap, double *);
  surface_flux       = (double *) va_arg(ap, double *);
  AdvectedEnergy     = (double *) va_arg(ap, double *);
  DeltaColdContent   = (double) va_arg(ap, double);
  Tfreeze            = (double) va_arg(ap, double);
  AvgCond            = (double) va_arg(ap, double);
  SWconducted        = (double) va_arg(ap, double);
  SnowDepth          = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfAttenuation    = (double) va_arg(ap, double);
  GroundFlux         = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  LatentHeatSub      = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  LWnet              = (double *) va_arg(ap, double *);
  ErrorString        = (char *) va_arg(ap, double *);
  
  /* print variables */
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: ice_melt failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  fprintf(stderr,"Dt = %f\n",Dt);
  fprintf(stderr,"Ra = %f\n",Ra);
  fprintf(stderr,"Ra_used = %f\n",*Ra_used);
  fprintf(stderr,"Z = %f\n",Z);
  fprintf(stderr,"Displacement = %f\n",Displacement);
  fprintf(stderr,"Z0 = %f\n",Z0);
  fprintf(stderr,"Wind = %f\n",Wind);
  fprintf(stderr,"ShortRad = %f\n",ShortRad);
  fprintf(stderr,"LongRadIn = %f\n",LongRadIn);
  fprintf(stderr,"AirDens = %f\n",AirDens);
  fprintf(stderr,"Lv = %f\n",Lv);
  fprintf(stderr,"Tair = %f\n",Tair);
  fprintf(stderr,"Press = %f\n",Press);
  fprintf(stderr,"Vpd = %f\n",Vpd);
  fprintf(stderr,"EactAir = %f\n",EactAir);
  fprintf(stderr,"Rain = %f\n",Rain);
  fprintf(stderr,"SweSurfaceLayer = %f\n",SweSurfaceLayer);
  fprintf(stderr,"SurfaceLiquidWater = %f\n",SurfaceLiquidWater);
  fprintf(stderr,"OldTSurf = %f\n",OldTSurf);
  fprintf(stderr,"RefreezeEnergy = %f\n",RefreezeEnergy[0]);
  fprintf(stderr,"vapor_flux = %f\n",*vapor_flux);
  fprintf(stderr,"blowing_flux = %f\n",*blowing_flux);
  fprintf(stderr,"surface_flux = %f\n",*surface_flux);
  fprintf(stderr,"AdvectedEnergy = %f\n",AdvectedEnergy[0]);
  fprintf(stderr,"DeltaColdContent = %f\n",DeltaColdContent);
  fprintf(stderr,"Tfreeze = %f\n",Tfreeze);
  fprintf(stderr,"AvgCond = %f\n",AvgCond);
  fprintf(stderr,"SWconducted = %f\n",SWconducted);
  fprintf(stderr,"SnowDepth = %f\n",SnowDepth);
  fprintf(stderr,"SnowDensity = %f\n",SnowDensity);
  fprintf(stderr,"SurfAttenuation = %f\n",SurfAttenuation);
  fprintf(stderr,"GroundFlux = %f\n",GroundFlux[0]);
  fprintf(stderr,"LatentHeat = %f\n",LatentHeat[0]);
  fprintf(stderr,"LatentHeatSub = %f\n",LatentHeatSub[0]);
  fprintf(stderr,"SensibleHeat = %f\n",SensibleHeat[0]);
  fprintf(stderr,"LWnet = %f\n",*LWnet);
  
  fprintf(stderr,"Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities.\n");

  return(ERROR);

}
