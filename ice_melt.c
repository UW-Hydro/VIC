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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#if LAKE_MODEL

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: IceMelt()

  Purpose      : Calculate snow accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    double delta_t               - Model timestep (hours)
    double z2           - Reference height (m) 
    double displacement          - Displacement height (m)
    double aero_resist  - Aerodynamic resistance (uncorrected for
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
    double atmos->melt           - Amount of snowpack outflow (m)
    double snow->pack_water      - Liquid water content of snow pack 
    double snow->surf_water	 - Liquid water content of surface layer 
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Comments     :

  Modifications:
  11-18-02 Modified method by which lake coverage fraction and ice height 
           are updated.                                                LCB
  04-Jun-04 Added descriptive error message to beginning of screen dump in
	    ErrorPrintIcePackEnergyBalance.				TJB

*****************************************************************************/
void ice_melt(double            z2,
	      double            aero_resist,
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
  int    Twidth;

  double DeltaPackCC;            /* Change in cold content of the pack */
  double DeltaPackSwq;           /* Change in snow water equivalent of the
                                   pack (m) */
  double InitialSwq;             /* Initial snow water equivalent (m) */
  double InitialIce;
  double MassBalanceError;       /* Mass balance error (m) */
  double MaxLiquidWater;         /* Maximum liquid water content of pack (m) */
  double OldTSurf;               /* Old snow surface temperature (C) */
  
  double Qnet;                   /* Net energy exchange at the surface (W/m2) */
  double RefreezeEnergy;         /* refreeze energy (W/m2) */
  double RefrozenWater;          /* Amount of refrozen water (m) */
  double SnowFallCC;             /* Cold content of new snowfall (J) */
  double SnowMelt;               /* Amount of snow melt during time interval
                                   (m water equivalent) */
  double IceMelt;
  double LWnet;
  double avgcond;
  double SWconducted;
  double SnowIce;
  double LakeIce;
  double SnowFall;
  double RainFall;
  double vapor_flux;
  double advection;
  double deltaCC;
  double SnowFlux;		/* thermal flux through snowpack from ground */
  double latent_heat;		
  double sensible_heat;		
  double melt_energy = 0.;

  SnowFall = snowfall / 1000.; /* convet to m */
  RainFall = rainfall / 1000.; /* convet to m */
  IceMelt = 0.0;
  RefrozenWater = 0.0;
  
  InitialSwq = snow->swq;
  OldTSurf = snow->surf_temp;

  /* Initialize snowpack variables, differs from snow_melt.c because 
  *  lake model only uses one layer */
 
  SnowIce  = snow->swq  - snow->surf_water;
  LakeIce = lake->hice * RHOICE/RHO_W;        /* meters of water equivalent. */
  InitialIce = LakeIce;
    
  /* Distribute fresh snowfall */
  SnowIce += SnowFall;
  snow->surf_water += RainFall;
  
  icerad (net_short, lake->hice, SnowIce*RHO_W/RHOSNOW, &avgcond, &SWconducted, &deltaCC);

  /* Calculate the surface energy balance for snow_temp = 0.0 */
  
  vapor_flux = snow->vapor_flux;

  Qnet = CalcIcePackEnergyBalance((double)0.0, (double)delta_t, aero_resist,
				  z2, displacement, Z0, wind, net_short, 
				  longwave, density, Le, air_temp,
				  pressure * 1000., vpd * 1000., vp * 1000.,
				  RainFall, snow->swq+LakeIce, 
				  snow->surf_water, 
				  OldTSurf, &RefreezeEnergy, &vapor_flux, 
				  &advection, deltaCC, Tcutoff, avgcond, 
				  SWconducted, snow->swq*RHO_W/RHOSNOW, 
				  RHOSNOW, surf_atten, &SnowFlux, 
				  &latent_heat, &sensible_heat, &LWnet);

  snow->vapor_flux = vapor_flux;
  save_refreeze_energy[0] = RefreezeEnergy;

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  if (Qnet == 0.0) {
    snow->surf_temp = 0.0;
    if (RefreezeEnergy >= 0.0) {                     /* Surface is freezing. */
      RefrozenWater = RefreezeEnergy/(Lf * RHO_W) 
          * delta_t * SECPHOUR; 
      if (RefrozenWater > snow->surf_water) {
        RefrozenWater = snow->surf_water;
        RefreezeEnergy = RefrozenWater * Lf * RHO_W/
          (delta_t * SECPHOUR);
      } 
      melt_energy  += RefreezeEnergy;
      snow->swq    += RefrozenWater;
      SnowIce      += RefrozenWater;
      snow->surf_water   -= RefrozenWater;
      assert(snow->surf_water >= 0.0);
      SnowMelt      = 0.0;
    }
    else {
      
      /* Calculate snow melt if refreeze energy is negative */      
      SnowMelt = fabs(RefreezeEnergy)/(Lf * RHO_W) * 
        delta_t * SECPHOUR;
      melt_energy += RefreezeEnergy;
    }

  
    /* Convert vapor mass flux to a depth per timestep and adjust snow->surf_water */
     snow->vapor_flux *= delta_t * SECPHOUR;
    if ((SnowIce + snow->surf_water + LakeIce) < -(snow->vapor_flux)) {
      snow->vapor_flux = -SnowIce - LakeIce;
      lake->volume -= LakeIce*fracprv*lake->surface[0];	    
      LakeIce=0.0;
      SnowIce=0.0; 
    }  
    else if ((SnowIce + snow->surf_water) < -(snow->vapor_flux) && (SnowIce + snow->surf_water +LakeIce) > -(snow->vapor_flux)) {
      LakeIce += (snow->vapor_flux+SnowIce);
      lake->volume += lake->surface[0]*fracprv*(snow->vapor_flux+SnowIce);
      SnowIce    = 0.0;
    }
    else {
      if(-(snow->vapor_flux) > snow->surf_water)
	{
	  SnowIce += (snow->vapor_flux + snow->surf_water);
	  snow->surf_water = 0.0;
	}
      else
	snow->surf_water += snow->vapor_flux;
    }
 
    /* If SnowMelt < SnowIce, there was incomplete melting of the pack */
    
    if (SnowMelt < SnowIce) {
        snow->surf_water += SnowMelt;
        SnowIce        -= SnowMelt;
    }
    /* Else if there was complete melting of the pack, but not of the ice. */
    else if(SnowMelt >= SnowIce && SnowMelt < (SnowIce + LakeIce)) {
      snow->surf_water += SnowIce;
      LakeIce -= (SnowMelt-SnowIce);
      IceMelt = (SnowMelt-SnowIce);
      SnowIce = 0.0;
    }
    /* Else, SnowMelt > TotalIce and there was complete melting of the ice */
    else {
      SnowMelt    = SnowIce + LakeIce;
      snow->surf_water += SnowIce;
      IceMelt = LakeIce;
      SnowIce         = 0.0;
      LakeIce         =0.0;
    } 
  }
  
  /* Else, IceEnergyBalance(T=0.0) <= 0.0 */
  else  {
    /* Calculate surface layer temperature using "Brent method" */

    vapor_flux = snow->vapor_flux;

    snow->surf_temp = root_brent((double)(snow->surf_temp-SNOW_DT), 
				 (double)0.0,
				 IceEnergyBalance, (double)delta_t, 
				 aero_resist, z2, 
				 displacement, Z0, wind, net_short, longwave,
				 density, Le, air_temp, pressure * 1000.,
				 vpd * 1000., vp * 1000., RainFall, 
				 snow->swq+LakeIce,
				 snow->surf_water, OldTSurf, &RefreezeEnergy, 
				 &vapor_flux, &advection, deltaCC, Tcutoff, 
				 avgcond, SWconducted,
				 snow->swq*RHO_W/RHOSNOW,
				 RHOSNOW,surf_atten,&SnowFlux,
				 &latent_heat, &sensible_heat, &LWnet);

    if(snow->surf_temp <= -9998)
      ErrorIcePackEnergyBalance(snow->surf_temp, (double)delta_t, aero_resist,
				 z2, displacement, Z0, wind, net_short,
				 longwave, density, Le, air_temp,
				 pressure * 1000., vpd * 1000., vp * 1000.,
				 RainFall, snow->swq+LakeIce, 
				snow->surf_water, 
				 OldTSurf,
				 &RefreezeEnergy, &vapor_flux, &advection, 
				 deltaCC, Tcutoff, 
				 avgcond, SWconducted,		
				 snow->swq*RHO_W/RHOSNOW, RHOSNOW,surf_atten,
				 &SnowFlux,&latent_heat,
				 &sensible_heat, &LWnet);

    Qnet = CalcIcePackEnergyBalance(snow->surf_temp, (double)delta_t, aero_resist,
				     z2, displacement, Z0, wind, net_short,
				     longwave, density, Le, air_temp,
				     pressure * 1000., vpd * 1000., 
				     vp * 1000.,RainFall, snow->swq+LakeIce, 
				     snow->surf_water, OldTSurf, &RefreezeEnergy, 
				     &vapor_flux, &advection, deltaCC, Tcutoff, 
				     avgcond, SWconducted, snow->swq*RHO_W/RHOSNOW, 
				     RHOSNOW,surf_atten, &SnowFlux,&latent_heat,
				     &sensible_heat, &LWnet);
    

    snow->vapor_flux = vapor_flux;
    save_refreeze_energy[0] = RefreezeEnergy;

    /* since we iterated, the surface layer is below freezing and no snowmelt
     */ 
    
    SnowMelt = 0.0;
    
    /* Since updated snow_temp < 0.0, all of the liquid water in the surface
       layer has been frozen */ 
    
    SnowIce        += snow->surf_water;
    melt_energy += snow->surf_water * Lf 
                       * RHO_W/(delta_t * SECPHOUR);
    RefrozenWater = snow->surf_water;
    snow->surf_water  = 0.0;
    
    /* Convert mass flux to a depth per timestep and adjust SurfaceSwq */
    
    snow->vapor_flux *= delta_t * SECPHOUR;
    if ((SnowIce + LakeIce) < -(snow->vapor_flux)) {
      snow->vapor_flux = -SnowIce - LakeIce;
      lake->volume -= lake->surface[0]*fracprv*LakeIce;	    
      LakeIce=0.0;
      SnowIce=0.0; 
    }  
    else if (SnowIce < -(snow->vapor_flux) && (SnowIce + LakeIce) > -(snow->vapor_flux)) {
      LakeIce += (snow->vapor_flux+SnowIce);
      lake->volume += lake->surface[0]*fracprv*(snow->vapor_flux+SnowIce);
      SnowIce    = 0.0;
    }
    else {
      if(SnowIce > 0.0)
	SnowIce += snow->vapor_flux;
      else
	lake->volume += lake->surface[0]*fracprv*snow->vapor_flux;
    }
  }

 
  /* Done with iteration etc, now Update the liquid water content of theN
     surface layer */ 
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * SnowIce;
  if  (snow->surf_water > MaxLiquidWater) {
    melt[0]    = snow->surf_water - MaxLiquidWater;
    snow->surf_water = MaxLiquidWater;
  }
  else
    melt[0] = 0.0;

  /* Update snow properties */
  
  snow->swq = SnowIce + snow->surf_water;
  lake->hice = LakeIce * RHO_W/RHOICE;
  if(lake->hice <= 0.0)
    {
      lake->hice = 0.0;
      lake->fraci = 0.0;
    }
  //  if (lake->hice < FRACMIN) {
      /* ....................................................................
       * If the ice height is lower than the minimum ice height increase the
       * height and decrease the fractional cover in order to conserve
       * numerical stability.  Ice covered fraction changes linearly up to FRACMIN.
       * ....................................................................*/

  //    lake->fraci=(lake->hice*lake->fraci)/FRACMIN;

  //  if(lake->fraci > 0.0)
  //lake->hice=FRACMIN;
  //  else {
  //lake->hice=0.0;
  //lake->fraci = 0.0;
  //  }
  //}

  /* Mass balance test */

    MassBalanceError = (InitialSwq - snow->swq) + (InitialIce - LakeIce) + (RainFall + SnowFall) 
    - IceMelt - melt[0] + snow->vapor_flux; 
    
    //printf("MassBalanceError = %g\n", MassBalanceError);

  melt[0] *= 1000.; /* converts back to mm */
  snow->mass_error = MassBalanceError;
  snow->vapor_flux *= -1.;
  *save_LWnet = LWnet;
  *save_advection = advection;
  *save_deltaCC = deltaCC;
  *save_SnowFlux = SnowFlux;
  *save_latent = latent_heat;
  *save_sensible = sensible_heat;
  *save_Qnet = Qnet;

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
  double *VaporMassFlux;          /* Mass flux of water vapor to or from the
                                   intercepted snow */
  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *DeltaColdContent;       /* Change in cold content (W/m2) */
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
  double *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */

  /* initialize variables */
    Dt                 = (double) va_arg(ap, double);
  Ra                 = (double) va_arg(ap, double);
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
  VaporMassFlux      = (double *) va_arg(ap, double *);
  AdvectedEnergy     = (double *) va_arg(ap, double *);
  DeltaColdContent   = (double *) va_arg(ap, double *);
  Tfreeze              = (double) va_arg(ap, double);
  AvgCond              = (double) va_arg(ap, double);
  SWconducted              = (double) va_arg(ap, double);
  SWabsorbed              = (double) va_arg(ap, double);
  SnowDepth          = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfAttenuation    = (double) va_arg(ap, double);
  GroundFlux         = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  
  /* print variables */
  fprintf(stderr, "ERROR: ice_melt failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  fprintf(stderr,"Dt = %f\n",Dt);
  fprintf(stderr,"Ra = %f\n",Ra);
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
  fprintf(stderr,"VaporMassFlux = %f\n",VaporMassFlux[0]);
  fprintf(stderr,"AdvectedEnergy = %f\n",AdvectedEnergy[0]);
  fprintf(stderr,"DeltaColdContent = %f\n",DeltaColdContent[0]);
  fprintf(stderr,"Tfreeze = %f\n",Tfreeze);
  fprintf(stderr,"AvgCond = %f\n",AvgCond);
  fprintf(stderr,"SWconducted = %f\n",SWconducted);
  fprintf(stderr,"SWabsorbed = %f\n",SWabsorbed);
  fprintf(stderr,"SnowDepth = %f\n",SnowDepth);
  fprintf(stderr,"SnowDensity = %f\n",SnowDensity);
  fprintf(stderr,"SurfAttenuation = %f\n",SurfAttenuation);
  fprintf(stderr,"GroundFlux = %f\n",GroundFlux[0]);
  fprintf(stderr,"LatentHeat = %f\n",LatentHeat[0]);
  fprintf(stderr,"SensibleHeat = %f\n",SensibleHeat[0]);
  
  fprintf(stderr,"Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities.\n");
  exit(0);

  return(0.0);

}

#endif // LAKE_MODEL
