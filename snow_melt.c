/*
 * SUMMARY:      SnowMelt.c - Calculate snow accumulation and melt
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Fri Oct 16 13:41:13 1998 by VIC Administrator <vicadmin@u.washington.edu>
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

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: SnowMelt()

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
*****************************************************************************/
void snow_melt(soil_con_struct   soil_con, 
	       int               rec,
	       int               iveg,
               double            z2,
               double            aero_resist,
               double            Le,
               snow_data_struct *snow, 
               double            delta_t,
               double            displacement,
               double            Z0,
               double            surf_atten,
               double            rainfall,
               double            snowfall,
               double            wind,
               double            grnd_surf_temp,
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
               double           *save_grnd_flux,
               double           *save_latent,
               double           *save_sensible,
               double           *save_Qnet,
	       double           *save_Trad,
	       double           *save_refreeze_energy)
{
  double DeltaPackCC;            /* Change in cold content of the pack */
  double DeltaPackSwq;           /* Change in snow water equivalent of the
                                   pack (m) */
  double Ice;                    /* Ice content of snow pack (m)*/
  double InitialSwq;             /* Initial snow water equivalent (m) */
  double MassBalanceError;       /* Mass balance error (m) */
  double MaxLiquidWater;         /* Maximum liquid water content of pack (m) */
  double OldTSurf;               /* Old snow surface temperature (C) */
  double PackCC;                 /* Cold content of snow pack (J) */
  double PackSwq;                /* Snow pack snow water equivalent (m) */
  double Qnet;                   /* Net energy exchange at the surface (W/m2) */
  double RefreezeEnergy;         /* refreeze energy (W/m2) */
  double RefrozenWater;          /* Amount of refrozen water (m) */
  double SnowFallCC;             /* Cold content of new snowfall (J) */
  double SnowMelt;               /* Amount of snow melt during time interval
                                   (m water equivalent) */
  double SurfaceCC;              /* Cold content of snow pack (J) */
  double SurfaceSwq;             /* Surface layer snow water equivalent (m) */
  double SnowFall;
  double RainFall;
  double vapor_flux;
  double advection;
  double deltaCC;
  double grnd_flux;		/* thermal flux through snowpack from ground */
  double latent_heat;		/* thermal flux through snowpack from ground */
  double sensible_heat;		/* thermal flux through snowpack from ground */
  double melt_energy = 0.;

  SnowFall = snowfall / 1000.; /* convet to m */
  RainFall = rainfall / 1000.; /* convet to m */

  InitialSwq = snow->swq;
  OldTSurf = snow->surf_temp;

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
  if (SnowFall > (MAX_SURFACE_SWE - SurfaceSwq)) {
    DeltaPackSwq = SurfaceSwq + SnowFall - MAX_SURFACE_SWE;
    if (DeltaPackSwq > SurfaceSwq)
      DeltaPackCC = SurfaceCC + (SnowFall - MAX_SURFACE_SWE)/SnowFall *
        SnowFallCC;
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
  
  vapor_flux = snow->vapor_flux;

  Qnet = CalcSnowPackEnergyBalance((double)0.0, delta_t, aero_resist,
         z2, displacement, Z0, wind, net_short,
         longwave, density, Le, air_temp,
         pressure * 1000., vpd * 1000., vp * 1000.,
         RainFall, SurfaceSwq, snow->surf_water, OldTSurf,
         &RefreezeEnergy, &vapor_flux, &advection, &deltaCC, grnd_surf_temp,
         snow->depth, snow->density,surf_atten,&grnd_flux,&latent_heat,
         &sensible_heat,save_Trad);

  snow->vapor_flux = vapor_flux;
  save_refreeze_energy[0] = RefreezeEnergy;

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  if (Qnet == 0.0) {
    snow->surf_temp = 0.0;
    if (RefreezeEnergy >= 0.0) {
      RefrozenWater = RefreezeEnergy/(Lf * RHO_W) 
          * delta_t * SECPHOUR; 
      if (RefrozenWater > snow->surf_water) {
        RefrozenWater = snow->surf_water;
        RefreezeEnergy = RefrozenWater * Lf * RHO_W/
          (delta_t * SECPHOUR);
      } 
      melt_energy  += RefreezeEnergy;
      SurfaceSwq   += RefrozenWater;
      Ice          += RefrozenWater;
      snow->surf_water   -= RefrozenWater;
      assert(snow->surf_water >= 0.0);
      SnowMelt      = 0.0;
    }
    else {
      
      /* Calculate snow melt */      
      SnowMelt = fabs(RefreezeEnergy)/(Lf * RHO_W) * 
        delta_t * SECPHOUR;
      melt_energy += RefreezeEnergy;
    }

    /* Convert vapor mass flux to a depth per timestep and adjust snow->surf_water */
    snow->vapor_flux *= delta_t * SECPHOUR;
    
    if (snow->surf_water < -(snow->vapor_flux)) {
      snow->vapor_flux = -(snow->surf_water);
      snow->surf_water    = 0.0;
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
    }
  }
  
  /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
  else  {
    /* Calculate surface layer temperature using "Brent method" */

    vapor_flux = snow->vapor_flux;

    snow->surf_temp = root_brent((double)(snow->surf_temp-DELTAT), 
				 (double)0.0,
				 SnowPackEnergyBalance, delta_t, 
				 aero_resist, z2, 
				 displacement, Z0, wind, net_short, longwave,
				 density, Le, air_temp, pressure * 1000.,
				 vpd * 1000., vp * 1000., RainFall, 
				 SurfaceSwq,
				 snow->surf_water, OldTSurf, &RefreezeEnergy, 
				 &vapor_flux, &advection, &deltaCC, 
				 grnd_surf_temp, snow->depth,
				 snow->density,surf_atten,&grnd_flux,
				 &latent_heat, &sensible_heat,
				 save_Trad);

    if(snow->surf_temp <= -9998)
      ErrorSnowPackEnergyBalance(snow->surf_temp, delta_t, aero_resist,
				 z2, displacement, Z0, wind, net_short,
				 longwave, density, Le, air_temp,
				 pressure * 1000., vpd * 1000., vp * 1000.,
				 RainFall, SurfaceSwq, snow->surf_water, 
				 OldTSurf,
				 &RefreezeEnergy, &vapor_flux, &advection, 
				 &deltaCC, grnd_surf_temp,
				 snow->depth, snow->density,surf_atten,
				 &grnd_flux,&latent_heat,
				 &sensible_heat,save_Trad);

    Qnet = CalcSnowPackEnergyBalance(snow->surf_temp, delta_t, aero_resist,
				     z2, displacement, Z0, wind, net_short,
				     longwave, density, Le, air_temp,
				     pressure * 1000., vpd * 1000., 
				     vp * 1000.,
				     RainFall, SurfaceSwq, snow->surf_water, 
				     OldTSurf,
				     &RefreezeEnergy, &vapor_flux, 
				     &advection, &deltaCC, grnd_surf_temp,
				     snow->depth, snow->density,surf_atten,
				     &grnd_flux,&latent_heat,
				     &sensible_heat,save_Trad);

    snow->vapor_flux = vapor_flux;
    save_refreeze_energy[0] = RefreezeEnergy;

    /* since we iterated, the surface layer is below freezing and no snowmelt
     */ 
    
    SnowMelt = 0.0;
    
    /* Since updated snow_temp < 0.0, all of the liquid water in the surface
       layer has been frozen */ 
    
    SurfaceSwq += snow->surf_water;
    Ice        += snow->surf_water;
    snow->surf_water  = 0.0;
    melt_energy += snow->surf_water * Lf 
                       * RHO_W/(delta_t * SECPHOUR);
    
    /* Convert mass flux to a depth per timestep and adjust SurfaceSwq */
    
    snow->vapor_flux *= delta_t * SECPHOUR;

    if (SurfaceSwq < -(snow->vapor_flux)) {
      snow->vapor_flux = -SurfaceSwq;
      SurfaceSwq    = 0.0;
      Ice           = PackSwq;
    }
    else {
      SurfaceSwq += snow->vapor_flux;
      Ice += snow->vapor_flux;
    }

  }
  
  /* Done with iteration etc, now Update the liquid water content of the
     surface layer */ 
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * SurfaceSwq;
  if  (snow->surf_water > MaxLiquidWater) {
    melt[0]    = snow->surf_water - MaxLiquidWater;
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
  RefreezeEnergy = snow->pack_water * Lf * RHO_W;

  /* calculate energy released to freeze*/
  
  if (PackCC < -RefreezeEnergy) { /* cold content not fully depleted*/
    PackSwq   += snow->pack_water;    /* refreeze all water and update*/
    Ice       += snow->pack_water;
    snow->pack_water = 0.0;
    if (PackSwq > 0.0) {
      PackCC = PackSwq * CH_ICE * snow->pack_temp + RefreezeEnergy;
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

    snow->pack_temp      = 0.0;    
    DeltaPackSwq = -PackCC/(Lf * RHO_W); 
    snow->pack_water -= DeltaPackSwq;
    PackSwq += DeltaPackSwq;
    Ice += DeltaPackSwq;
  }
  
  /* Update the liquid water content of the pack */
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * PackSwq;
  if (snow->pack_water > MaxLiquidWater) {
    melt[0]    = snow->pack_water - MaxLiquidWater;
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

  melt[0] *= 1000.; /* converts back to mm */
  snow->mass_error = MassBalanceError;
  snow->coldcontent = SurfaceCC;
  snow->vapor_flux *= -1.;
  *save_advection = advection;
  *save_deltaCC = deltaCC;
  *save_grnd_flux = grnd_flux;
  *save_latent = latent_heat;
  *save_sensible = sensible_heat;
  *save_Qnet = Qnet;

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
  char *Routine = "CalcSnowPackEnergyBalance";

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
  char *Routine = "CalcSnowPackEnergyBalance";

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
  double *GroundFlux;		/* Ground Heat Flux (W/m2) */
  double *RefreezeEnergy;        /* Refreeze energy (W/m2) */
  double *VaporMassFlux;          /* Mass flux of water vapor to or from the
                                   intercepted snow */
  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *DeltaColdContent;       /* Change in cold content (W/m2) */
  double TGrnd;     
  double SnowDepth;  
  double SnowDensity; 
  double SurfAttenuation; 
  
  /* end of list of arguments in variable argument list */

  double *LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  double *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */
  double *TMean;                /* Average temperature for time step (C) */

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
  TGrnd              = (double) va_arg(ap, double);
  SnowDepth          = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfAttenuation    = (double) va_arg(ap, double);
  GroundFlux         = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  TMean              = (double *) va_arg(ap, double *);
  
  /* print variables */
  fprintf(stderr,"Dt = %lf\n",Dt);
  fprintf(stderr,"Ra = %lf\n",Ra);
  fprintf(stderr,"Z = %lf\n",Z);
  fprintf(stderr,"Displacement = %lf\n",Displacement);
  fprintf(stderr,"Z0 = %lf\n",Z0);
  fprintf(stderr,"Wind = %lf\n",Wind);
  fprintf(stderr,"ShortRad = %lf\n",ShortRad);
  fprintf(stderr,"LongRadIn = %lf\n",LongRadIn);
  fprintf(stderr,"AirDens = %lf\n",AirDens);
  fprintf(stderr,"Lv = %lf\n",Lv);
  fprintf(stderr,"Tair = %lf\n",Tair);
  fprintf(stderr,"Press = %lf\n",Press);
  fprintf(stderr,"Vpd = %lf\n",Vpd);
  fprintf(stderr,"EactAir = %lf\n",EactAir);
  fprintf(stderr,"Rain = %lf\n",Rain);
  fprintf(stderr,"SweSurfaceLayer = %lf\n",SweSurfaceLayer);
  fprintf(stderr,"SurfaceLiquidWater = %lf\n",SurfaceLiquidWater);
  fprintf(stderr,"OldTSurf = %lf\n",OldTSurf);
  fprintf(stderr,"RefreezeEnergy = %lf\n",RefreezeEnergy[0]);
  fprintf(stderr,"VaporMassFlux = %lf\n",VaporMassFlux[0]);
  fprintf(stderr,"AdvectedEnergy = %lf\n",AdvectedEnergy[0]);
  fprintf(stderr,"DeltaColdContent = %lf\n",DeltaColdContent[0]);
  fprintf(stderr,"TGrnd = %lf\n",TGrnd);
  fprintf(stderr,"SnowDepth = %lf\n",SnowDepth);
  fprintf(stderr,"SnowDensity = %lf\n",SnowDensity);
  fprintf(stderr,"SurfAttenuation = %lf\n",SurfAttenuation);
  fprintf(stderr,"GroundFlux = %lf\n",GroundFlux[0]);
  fprintf(stderr,"LatentHeat = %lf\n",LatentHeat[0]);
  fprintf(stderr,"SensibleHeat = %lf\n",SensibleHeat[0]);
  fprintf(stderr,"TMean = %lf\n",TMean[0]);
  
  vicerror("Finished dumping snow_melt variables");

  return(0.0);

}

