/*
 * SUMMARY:      SnowMelt.c - Calculate snow accumulation and melt
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Tue Apr 22 09:48:43 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
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
  Function name: SnowMelt()

  Purpose      : Calculate snow accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    double delta_t               - Model timestep (secs)
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
    double snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m/time step)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Modifies     :
    double *melt                 - Amount of snowpack outflow (initially is m, but converted to mm for output)
    double snow->pack_water      - Liquid water content of snow pack 
    double snow->surf_water	 - Liquid water content of surface layer 
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m/time step)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Comments     :

  Modifications:
  10-06-00 modified to handle partial snow cover                KAC
  10-31-00 modified to assure that ground heat flux is used
           properly in the snow surface energy balance as well
           as imporving the handling of energy fluxes for 
           partial snow cover.                                  KAC
  11-18-02 modified to handle blowing snow.                     LCB
  04-Jun-04 For the case in which snowpack is too thin to solve
	    separately, added message explaining that root_brent's
	    error is not fatal and that snow pack will be solved in
	    conjunction with surface energy balance.			TJB
  16-Jul-04 Added "month" to parameter list to allow this function to
	    pass month to latent_heat_from_snow().  Changed calculations
	    involving vapor_flux to make it consistent with convention
	    that vapor_flux has units of m/timestep.			TJB
  16-Jul-04 Changed the type of the last few variables (lag_one, iveg,
	    etc) to be double in the parameter lists of root_brent
	    and CalcSnowPackEnergyBalance.  For some reason, passing
	    them as float or int caused them to become garbage.  This may
	    have to do with the fact that they followed variables of type
	    (double *) in va_list, which may have caused memory alignment
	    problems.							TJB
  16-Jul-04 Modified cap on vapor_flux to re-scale values of blowing_flux
	    and surface_flux so that blowing_flux and surface_flux still
	    add up to the new value of vapor_flux.			TJB
  05-Aug-04 Removed overstory, lag_one, sigma_slope, fetch, iveg, Nveg,
	    and month from argument list, since these were only used in
	    call to SnowPackEnergyBalance (and ErrorSnowPackEnergyBalance),
	    which no longer needs them.					TJB
  25-Aug-04 Modified re-scaling of surface_flux to reduce round-off
	    error.							TJB
  21-Sep-04 Added ErrorString to store error messages from
	    root_brent.							TJB
	    Removed message explaining non-fatal root_brent warning.
	    These warnings were not a sign of any failures, and served
	    only to confuse users and take up valuable space in the
	    output display.						TJB
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    used in flux calculations.					TJB
  2007-Apr-11 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.	KAC via TJB
  2007-Jul-03 Corrected the units of melt in the comment section.	TJB
  2007-Aug-31 Checked root_brent return value against -998 rather than -9998.    JCA
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
*****************************************************************************/
int  snow_melt(double            Le, 
               double            NetShortSnow,  // net SW at absorbed by snow
	       double            Tcanopy,
	       double            Tgrnd,
	       double            *Z0,  // roughness
               double            aero_resist,  // aerodynamic resistance
               double            *aero_resist_used,  // stability-corrected aerodynamic resistance
               double            air_temp,  // air temperature
	       double            coverage, // snowpack cover fraction
	       double            delta_t,  // time step in secs
               double            density,  // atmospheric density
	       double            displacement,  // surface displacement
               double            grnd_flux,  // ground heat flux
	       double            LongSnowIn,  // incoming longwave radiation
	       double            pressure,
	       double            rainfall,
               double            snowfall,
               double            vp,
	       double            vpd,
	       double            wind,
               double            z2,
	       double           *NetLongSnow, 
               double           *OldTSurf,
	       double           *melt,
               double           *save_Qnet, 
               double           *save_advected_sensible,
               double           *save_advection,
               double           *save_deltaCC,
               double           *save_grnd_flux,
               double           *save_latent,
               double           *save_latent_sub,
               double           *save_refreeze_energy,
               double           *save_sensible,
	       int               UNSTABLE_SNOW,
	       int               rec,
	       int               iveg,
	       int               band,
               snow_data_struct *snow,
	       soil_con_struct  *soil_con)
{
  extern option_struct   options;
  int    Twidth;
  double error;
  double DeltaPackCC;            /* Change in cold content of the pack */
  double DeltaPackSwq;           /* Change in snow water equivalent of the
                                   pack (m) */
  double Ice;                    /* Ice content of snow pack (m)*/
  double InitialSwq;             /* Initial snow water equivalent (m) */
  double MassBalanceError;       /* Mass balance error (m) */
  double MaxLiquidWater;         /* Maximum liquid water content of pack (m) */
  double PackCC;                 /* Cold content of snow pack (J) */
  double PackSwq;                /* Snow pack snow water equivalent (m) */
  double Qnet;                   /* Net energy exchange at the surface (W/m2) */
  double RefreezeEnergy;         /* refreeze/melt energy in surface layer (W/m2) */
  double PackRefreezeEnergy;     /* refreeze/melt energy in pack layer (W/m2) */
  double RefrozenWater;          /* Amount of refrozen water (m) */
  double SnowFallCC;             /* Cold content of new snowfall (J) */
  double SnowMelt;               /* Amount of snow melt during time interval
                                   (m water equivalent) */
  double SurfaceCC;              /* Cold content of snow pack (J) */
  double SurfaceSwq;             /* Surface layer snow water equivalent (m) */
  double SnowFall;
  double RainFall;
  double advection;
  double deltaCC;
  double latent_heat;		
  double latent_heat_sub;		
  double sensible_heat;		
  double advected_sensible_heat;		
  double melt_energy = 0.;

  char ErrorString[MAXSTRING];

  SnowFall = snowfall / 1000.; /* convet to m */
  RainFall = rainfall / 1000.; /* convet to m */

  InitialSwq  = snow->swq;
  (*OldTSurf) = snow->surf_temp;

  /* Initialize snowpack variables */
  
  Ice  = snow->swq - snow->pack_water - snow->surf_water;
  
  /* Reconstruct snow pack */
  if (Ice > MAX_SURFACE_SWE)
    SurfaceSwq = MAX_SURFACE_SWE;
  else
    SurfaceSwq = Ice;
  PackSwq = Ice - SurfaceSwq;
  
  /* Calculate cold contents */
  SurfaceCC = CH_ICE * SurfaceSwq * snow->surf_temp;
  PackCC = CH_ICE * PackSwq * snow->pack_temp;
  if (air_temp > 0.0)
    SnowFallCC = 0.0;
  else
    SnowFallCC = CH_ICE * SnowFall * air_temp;
  
  /* Distribute fresh snowfall */
  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq) && (MAX_SURFACE_SWE - SurfaceSwq) > SMALL) {
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
  }
  if (SurfaceSwq > 0.0)
    snow->surf_temp = SurfaceCC/(CH_ICE * SurfaceSwq);
  else 
    snow->surf_temp = 0.0;
  if (PackSwq > 0.0)    
    snow->pack_temp = PackCC/(CH_ICE * PackSwq);
  else
    snow->pack_temp = 0.0;

  /* Adjust ice and snow->surf_water */
  Ice += SnowFall;
  snow->surf_water += RainFall;
  
  /* Calculate the surface energy balance for snow_temp = 0.0 */
  
  Qnet = CalcSnowPackEnergyBalance((double)0.0, delta_t, aero_resist, aero_resist_used,
				   displacement, z2, Z0, 
				   density, vp, LongSnowIn, Le, pressure,
				   RainFall, NetShortSnow, vpd, 
				   wind, (*OldTSurf), coverage, 
				   snow->depth, snow->density,
				   snow->surf_water, SurfaceSwq, 
				   Tcanopy, Tgrnd, 
				   &advection, &advected_sensible_heat, 
				   &deltaCC, 
				   &grnd_flux, &latent_heat, 
				   &latent_heat_sub, NetLongSnow, 
				   &RefreezeEnergy, &sensible_heat, 
				   &snow->vapor_flux, &snow->blowing_flux,
				   &snow->surface_flux);

  /* Check that snow swq exceeds minimum value for model stability */
//  if ( SurfaceSwq > MIN_SWQ_EB_THRES && !UNSTABLE_SNOW ) {
  if ( !UNSTABLE_SNOW ) {

    /* If Qnet == 0.0, then set the surface temperature to 0.0 */
    if (Qnet == 0.0) {
      snow->surf_temp = 0.0;
      if (RefreezeEnergy >= 0.0) {
	RefrozenWater = RefreezeEnergy / (Lf * RHO_W) * delta_t; 
	if (RefrozenWater > snow->surf_water) {
	  RefrozenWater = snow->surf_water;
	  RefreezeEnergy = RefrozenWater * Lf * RHO_W/ (delta_t);
	} 
	melt_energy  += RefreezeEnergy;
	SurfaceSwq   += RefrozenWater;
	Ice          += RefrozenWater;
	snow->surf_water   -= RefrozenWater;
	if (snow->surf_water < 0.0) snow->surf_water = 0.0;
	SnowMelt      = 0.0;
      }
      else {
	
	/* Calculate snow melt */      
	SnowMelt = fabs(RefreezeEnergy)/(Lf * RHO_W) * delta_t;
	melt_energy += RefreezeEnergy;
      }
      
      /* Adjust snow->surf_water for vapor_flux */
      if (snow->surf_water < -(snow->vapor_flux)) {
        // if vapor_flux exceeds surf_water, we not only need to
        // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
//        snow->surface_flux *= -( snow->surf_water / snow->vapor_flux );
        snow->blowing_flux *= -( snow->surf_water / snow->vapor_flux );
	snow->vapor_flux = -(snow->surf_water);
        snow->surface_flux = -(snow->surf_water) - snow->blowing_flux;
	snow->surf_water = 0.0;
      }
      else
	snow->surf_water += snow->vapor_flux;
      
      /* If SnowMelt < Ice, there was incomplete melting of the pack */
      
      if (SnowMelt < Ice) {
	if (SnowMelt <= PackSwq) {
	  snow->surf_water += SnowMelt;
	  PackSwq    -= SnowMelt;
	  Ice        -= SnowMelt;
	}
	else {
	  snow->surf_water += SnowMelt + snow->pack_water;
	  snow->pack_water  = 0.0;
	  PackSwq     = 0.0;
	  Ice        -= SnowMelt;
	  SurfaceSwq  = Ice;
	}
      }
      
      /* Else, SnowMelt > Ice and there was complete melting of the pack */
      else {
	SnowMelt    = Ice;
	snow->surf_water += Ice;
	SurfaceSwq  = 0.0;
	snow->surf_temp      = 0.0;
	PackSwq     = 0.0;
	snow->pack_temp      = 0.0;
	Ice         = 0.0;
	/* readjust melt energy to account for melt only of available snow */
	melt_energy -= RefreezeEnergy;
	RefreezeEnergy = RefreezeEnergy / fabs(RefreezeEnergy) 
	  * SnowMelt * Lf * RHO_W / ( delta_t );
	melt_energy += RefreezeEnergy;
      }
    }
  
    /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
    else  {
      /* Calculate surface layer temperature using "Brent method" */
      if (SurfaceSwq > MIN_SWQ_EB_THRES) {
	snow->surf_temp = root_brent((double)(snow->surf_temp-SNOW_DT), 
				     (double)(snow->surf_temp+SNOW_DT),
				     ErrorString, SnowPackEnergyBalance, 
				     delta_t, aero_resist, aero_resist_used,
				     displacement, z2, Z0, 
				     density, vp, LongSnowIn, Le, pressure,
				     RainFall, NetShortSnow, vpd, 
				     wind, (*OldTSurf), coverage, 
				     snow->depth, snow->density, 
				     snow->surf_water, SurfaceSwq, 
				     Tcanopy, Tgrnd, 
				     &advection, &advected_sensible_heat, 
				     &deltaCC, 
				     &grnd_flux, &latent_heat, 
				     &latent_heat_sub, NetLongSnow, 
				     &RefreezeEnergy, &sensible_heat, 
				     &snow->vapor_flux, &snow->blowing_flux,
				     &snow->surface_flux);
      
        if (snow->surf_temp <= -998) {
          if (options.TFALLBACK) {
            snow->surf_temp = *OldTSurf;
            snow->surf_temp_fbflag = 1;
            snow->surf_temp_fbcount++;
          }
          else {
	    error = ErrorSnowPackEnergyBalance(snow->surf_temp, rec, iveg, band,
					 delta_t, aero_resist, aero_resist_used,
					 displacement, z2, Z0, 
					 density, vp, LongSnowIn, Le, pressure,
					 RainFall, NetShortSnow, vpd, 
					 wind, (*OldTSurf), coverage, 
					 snow->depth, snow->density,
					 snow->surf_water, SurfaceSwq, 
					 Tcanopy, Tgrnd, 
					 &advection, &advected_sensible_heat, 
					 &deltaCC, 
					 &grnd_flux, &latent_heat, 
					 &latent_heat_sub, 
					 NetLongSnow, &RefreezeEnergy, 
					 &sensible_heat, &snow->vapor_flux,
					 &snow->blowing_flux, &snow->surface_flux);
            return(ERROR);
          }
        }
      }
      else {
	/* Thin snowpack must be solved in conjunction with ground surface energy balance */
//        fprintf(stderr,"Snowpack is too thin to solve separately; it will be solved in conjunction with ground surface energy balance\n");
	snow->surf_temp = 999;
      }
      if (snow->surf_temp > -998 && snow->surf_temp < 999) {
	Qnet = CalcSnowPackEnergyBalance(snow->surf_temp, 
					 delta_t, aero_resist, aero_resist_used,
					 displacement, z2, Z0, 
					 density, vp, LongSnowIn, Le, pressure,
					 RainFall, NetShortSnow, vpd, 
					 wind, (*OldTSurf), coverage, 
					 snow->depth, snow->density,
					 snow->surf_water, SurfaceSwq, 
					 Tcanopy, Tgrnd, 
					 &advection, &advected_sensible_heat, 
					 &deltaCC, 
					 &grnd_flux, &latent_heat, 
					 &latent_heat_sub, 
					 NetLongSnow, &RefreezeEnergy, 
					 &sensible_heat, &snow->vapor_flux,
					 &snow->blowing_flux, &snow->surface_flux);
	
	/* since we iterated, the surface layer is below freezing and no snowmelt */ 
	
	SnowMelt = 0.0;
	
	/* Since updated snow_temp < 0.0, all of the liquid water in the surface
	   layer has been frozen */ 
	
	SurfaceSwq += snow->surf_water;
	Ice        += snow->surf_water;
	snow->surf_water  = 0.0;
	melt_energy += snow->surf_water * Lf * RHO_W/(delta_t);
	
	/* Adjust SurfaceSwq for vapor flux */
	if (SurfaceSwq < -(snow->vapor_flux)) {
          // if vapor_flux exceeds SurfaceSwq, we not only need to
          // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
//          snow->surface_flux *= -( SurfaceSwq / snow->vapor_flux );
          snow->blowing_flux *= -( SurfaceSwq / snow->vapor_flux );
	  snow->vapor_flux = -SurfaceSwq;
          snow->surface_flux = -SurfaceSwq - snow->blowing_flux;
	  SurfaceSwq    = 0.0;
	  Ice           = PackSwq;
	}
	else {
	  SurfaceSwq += snow->vapor_flux;
	  Ice += snow->vapor_flux;
	}
      }
    }
  }
  else {
    /* Snow solution is unstable as independent layer */
/*     RefreezeEnergy     = -Ice * Lf * RHO_W / delta_t; */
/*     melt_energy        = -RefreezeEnergy; */
/*     SnowMelt           = Ice; */
/*     snow->surf_water  += Ice; */
/*     SurfaceSwq         = 0.0; */
/*     snow->surf_temp    = 0.0; */
/*     PackSwq            = 0.0; */
/*     snow->pack_temp    = 0.0; */
/*     Ice                = 0.0; */
/*     latent_heat        = 0.0; */
/*     latent_heat_sub    = 0.0; */
/*     sensible_heat      = 0.0; */
/*     snow->vapor_flux   = 0.0; */
    snow->surf_temp = 999;
  }

  /* Done with iteration etc, now Update the liquid water content of the
     surface layer */ 
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * SurfaceSwq;
  if  (snow->surf_water > MaxLiquidWater) {
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
  
  snow->pack_water += melt[0]; /* add surface layer outflow to pack 
                                      liquid water*/
  PackRefreezeEnergy = snow->pack_water * Lf * RHO_W;

  /* calculate energy released to freeze*/
  
  if (PackCC < -PackRefreezeEnergy) { /* cold content not fully depleted*/
    PackSwq   += snow->pack_water;    /* refreeze all water and update*/
    Ice       += snow->pack_water;
    snow->pack_water = 0.0;
    if (PackSwq > 0.0) {
      PackCC = PackSwq * CH_ICE * snow->pack_temp + PackRefreezeEnergy;
      snow->pack_temp = PackCC / (CH_ICE * PackSwq);
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
  
  Ice  = PackSwq + SurfaceSwq;

  if (Ice > MAX_SURFACE_SWE) {
    SurfaceCC   = CH_ICE * snow->surf_temp * SurfaceSwq;
    PackCC      = CH_ICE * snow->pack_temp * PackSwq;
    if (SurfaceSwq > MAX_SURFACE_SWE) {
      PackCC     += SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq; 
      SurfaceCC  -= SurfaceCC * (SurfaceSwq - MAX_SURFACE_SWE) / SurfaceSwq; 
      PackSwq    += SurfaceSwq - MAX_SURFACE_SWE;
      SurfaceSwq -= SurfaceSwq - MAX_SURFACE_SWE;
    }
    else if ( SurfaceSwq < MAX_SURFACE_SWE) {
      PackCC     -= PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq; 
      SurfaceCC  += PackCC * (MAX_SURFACE_SWE - SurfaceSwq) / PackSwq; 
      PackSwq    -= MAX_SURFACE_SWE - SurfaceSwq;
      SurfaceSwq += MAX_SURFACE_SWE - SurfaceSwq;
    }
    snow->pack_temp      = PackCC / (CH_ICE * PackSwq);
    snow->surf_temp      = SurfaceCC / (CH_ICE * SurfaceSwq);
  }
  else {
    PackSwq = 0.0;
    PackCC  = 0.0;
    snow->pack_temp  = 0.0;
  }

  snow->swq = Ice + snow->pack_water + snow->surf_water;

  if (snow->swq == 0.0) {
    snow->surf_temp = 0.0;
    snow->pack_temp = 0.0;
  }
    
  /* Mass balance test */
  
  MassBalanceError = (InitialSwq - snow->swq) + (RainFall + SnowFall) 
    - melt[0] + snow->vapor_flux; 
  
/*  printf("%d %d %g\n", y, x, MassBalanceError);*/

  melt[0]                 *= 1000.; /* converts back to mm */
  snow->mass_error         = MassBalanceError;
  snow->coldcontent        = SurfaceCC;
  snow->vapor_flux        *= -1.;
  *save_advection          = advection;
  *save_deltaCC            = deltaCC;
  *save_grnd_flux          = grnd_flux;
  *save_latent             = latent_heat;
  *save_latent_sub         = latent_heat_sub;
  *save_sensible           = sensible_heat;
  *save_advected_sensible  = advected_sensible_heat;
  *save_refreeze_energy    = RefreezeEnergy;
  *save_Qnet               = Qnet;

  return ( 0 );
}

/*****************************************************************************
  Function name: CalcSnowPackEnergyBalance()

  Purpose      : Dummy function to make a direct call to
                 SnowEnergyBalance() possible.

  Required     : 
    double TSurf - SnowPack surface temperature (C)
    other arguments required by SnowPackEnergyBalance()

  Returns      :
    double Qnet - Net energy exchange at the SnowPack snow surface (W/m^2)

  Modifies     : none

  Comments     : function is local to this module
*****************************************************************************/
double CalcSnowPackEnergyBalance(double Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  double Qnet;                   /* Net energy exchange at the SnowPack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = SnowPackEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

double ErrorSnowPackEnergyBalance(double Tsurf, ...)
{

  va_list ap;                   /* Used in traversing variable argument list
                                 */ 
  double Qnet;                   /* Net energy exchange at the SnowPack snow
                                   surface (W/m^2) */

  va_start(ap, Tsurf);
  Qnet = ErrorPrintSnowPackEnergyBalance(Tsurf, ap);
  va_end(ap);
  
  return Qnet;
}

double ErrorPrintSnowPackEnergyBalance(double TSurf, va_list ap)
{

  /* Define Variable Argument List */

  /* General Model Parameters */
  int rec, iveg, band;
  double Dt;                      /* Model time step (sec) */

  /* Vegetation Parameters */
  double Ra;                      /* Aerodynamic resistance (s/m) */
  double Displacement;            /* Displacement height (m) */
  double Z;                       /* Reference height (m) */
  double Z0;                      /* surface roughness height (m) */

  /* Atmospheric Forcing Variables */
  double AirDens;                 /* Density of air (kg/m3) */
  double EactAir;                 /* Actual vapor pressure of air (Pa) */
  double LongSnowIn;              /* Incoming longwave radiation (W/m2) */
  double Lv;                      /* Latent heat of vaporization (J/kg3) */
  double Press;                   /* Air pressure (Pa) */
  double Rain;                    /* Rain fall (m/timestep) */
  double ShortRad;                /* Net incident shortwave radiation 
				     (W/m2) */ 
  double Vpd;			  /* Vapor pressure deficit (Pa) */
  double Wind;                    /* Wind speed (m/s) */

  /* Snowpack Variables */
  double OldTSurf;                /* Surface temperature during previous time
				     step */ 
  double SnowCoverFract;          /* Fraction of area covered by snow */
  double SnowDensity;             /* Density of snowpack (kg/m^3) */
  double SurfaceLiquidWater;      /* Liquid water in the surface layer (m) */
  double SweSurfaceLayer;         /* Snow water equivalent in surface layer 
				     (m) */ 

  /* Energy Balance Components */
  double Tair;                    /* Canopy surface temperature (C) */
  double TGrnd;                   /* Ground surface temperature (C) */

  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *AdvectedSensibleHeat;	  /* Sensible heat advected from snow-free
				     area into snow covered area (W/m^2) */
  double *DeltaColdContent;       /* Change in cold content of surface 
				     layer (W/m2) */
  double *DeltaPackColdContent;   /* Change in sold content of pack 
				     layer (W/m^2) */
  double *GroundFlux;		  /* Ground Heat Flux (W/m2) */
  double *LatentHeat;		  /* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;          /* Latent heat of sub exchange at 
				     surface (W/m2) */
  double *NetLongSnow;            /* Net longwave radiation at snowpack 
				     surface (W/m^2) */
  double *RefreezeEnergy;         /* Refreeze energy (W/m2) */
  double *SensibleHeat;		  /* Sensible heat exchange at surface 
				     (W/m2) */
  double *VaporMassFlux;          /* Mass flux of water vapor to or from the
				     intercepted snow */
  double *BlowingMassFlux;          /* Mass flux of water vapor to or from the
				     intercepted snow */
  double *SurfaceMassFlux;          /* Mass flux of water vapor to or from the
					 intercepted snow */

  char *ErrorString;

  /* Read Variable Argument List */

  /* General Model Parameters */
  rec       = (int) va_arg(ap, int);
  iveg      = (int) va_arg(ap, int);
  band      = (int) va_arg(ap, int);
  Dt        = (double) va_arg(ap, double);

  /* Vegetation Parameters */
  Ra           = (double) va_arg(ap, double);
  Displacement = (double) va_arg(ap, double);
  Z            = (double) va_arg(ap, double);
  Z0           = (double) va_arg(ap, double);

  /* Atmospheric Forcing Variables */
  AirDens    = (double) va_arg(ap, double);
  EactAir    = (double) va_arg(ap, double);
  LongSnowIn = (double) va_arg(ap, double);
  Lv         = (double) va_arg(ap, double);
  Press      = (double) va_arg(ap, double);
  Rain       = (double) va_arg(ap, double);
  ShortRad   = (double) va_arg(ap, double);
  Vpd        = (double) va_arg(ap, double);
  Wind       = (double) va_arg(ap, double);

  /* Snowpack Variables */
  OldTSurf           = (double) va_arg(ap, double);
  SnowCoverFract     = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfaceLiquidWater = (double) va_arg(ap, double);
  SweSurfaceLayer    = (double) va_arg(ap, double);

  /* Energy Balance Components */
  Tair         = (double) va_arg(ap, double);
  TGrnd           = (double) va_arg(ap, double);

  AdvectedEnergy        = (double *) va_arg(ap, double *);
  AdvectedSensibleHeat  = (double *)va_arg(ap, double *);
  DeltaColdContent      = (double *) va_arg(ap, double *);
  DeltaPackColdContent  = (double *) va_arg(ap, double *);
  GroundFlux            = (double *) va_arg(ap, double *);
  LatentHeat            = (double *) va_arg(ap, double *);
  LatentHeatSub         = (double *) va_arg(ap, double *);
  NetLongSnow           = (double *) va_arg(ap, double *);
  RefreezeEnergy        = (double *) va_arg(ap, double *);
  SensibleHeat          = (double *) va_arg(ap, double *);
  VaporMassFlux         = (double *) va_arg(ap, double *);
  BlowingMassFlux       = (double *) va_arg(ap, double *);
  SurfaceMassFlux       = (double *) va_arg(ap, double *);
  ErrorString           = (char *) va_arg(ap, char *);

  /* print variables */
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: snow_melt failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  /* general model terms */
  fprintf(stderr, "rec = %i\n", rec);
  fprintf(stderr, "iveg = %i\n", iveg);
  fprintf(stderr, "band = %i\n", band);
  fprintf(stderr, "Dt = %f\n",Dt);

  /* land surface parameters */
  fprintf(stderr,"Ra = %f\n",Ra);
  fprintf(stderr,"Displacement = %f\n",Displacement);
  fprintf(stderr,"Z = %f\n",Z);
  fprintf(stderr,"Z0 = %f\n",Z0);

  /* meteorological terms */
  fprintf(stderr,"AirDens = %f\n",AirDens);
  fprintf(stderr,"EactAir = %f\n",EactAir);
  fprintf(stderr,"LongSnowIn = %f\n",LongSnowIn);
  fprintf(stderr,"Lv = %f\n",Lv);
  fprintf(stderr,"Press = %f\n",Press);
  fprintf(stderr,"Rain = %f\n",Rain);
  fprintf(stderr,"ShortRad = %f\n",ShortRad);
  fprintf(stderr,"Vpd = %f\n",Vpd);
  fprintf(stderr,"Wind = %f\n",Wind);

  /* snow pack terms */
  fprintf(stderr,"OldTSurf = %f\n",OldTSurf);
  fprintf(stderr,"SnowCoverFract = %f\n", SnowCoverFract);
  fprintf(stderr,"SnowDensity = %f\n",SnowDensity);
  fprintf(stderr,"SurfaceLiquidWater = %f\n",SurfaceLiquidWater);
  fprintf(stderr,"SweSurfaceLayer = %f\n",SweSurfaceLayer);
  fprintf(stderr,"Tair = %f\n",Tair);
  fprintf(stderr,"TGrnd = %f\n",TGrnd);
  fprintf(stderr,"AdvectedEnergy = %f\n",AdvectedEnergy[0]);
  fprintf(stderr,"AdvectedSensibleHeat = %f\n",AdvectedSensibleHeat[0]);
  fprintf(stderr,"DeltaColdContent = %f\n",DeltaColdContent[0]);
  fprintf(stderr,"DeltaPackColdContent = %f\n",DeltaPackColdContent[0]);
  fprintf(stderr,"GroundFlux = %f\n",GroundFlux[0]);
  fprintf(stderr,"LatentHeat = %f\n",LatentHeat[0]);
  fprintf(stderr,"LatentHeatSub = %f\n",LatentHeatSub[0]);
  fprintf(stderr,"NetLongSnow = %f\n",NetLongSnow[0]);
  fprintf(stderr,"RefreezeEnergy = %f\n",RefreezeEnergy[0]);
  fprintf(stderr,"SensibleHeat = %f\n",SensibleHeat[0]);
  fprintf(stderr,"VaporMassFlux = %f\n",VaporMassFlux[0]);
  fprintf(stderr,"BlowingMassFlux = %f\n",BlowingMassFlux[0]);
  fprintf(stderr,"SurfaceMassFlux = %f\n",SurfaceMassFlux[0]);

  fprintf(stderr,"Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThencheck output for instabilities.\n");

  return(ERROR);

}

