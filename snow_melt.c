/*
 * SUMMARY:      SnowMelt.c - Calculate snow accumulation and melt
 * USAGE:        
 *
 * AUTHOR:       Mark Wigmosta and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Sat Apr 18 15:47:39 1998 by Keith Aric Cherkauer <cherkaue@u.washington.edu>
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
void snow_melt(atmos_data_struct *atmos, 
               soil_con_struct soil_con, 
               double z2,
               double aero_resist,
               double Le,
               snow_data_struct *snow, 
               double delta_t,
               double displacement,
               double Z0,
               double surf_atten,
               double rainfall,
               double snowfall,
               double wind,
               double grnd_surf_temp,
               double *save_advection,
               double *save_deltaCC,
               double *save_grnd_flux,
               double *save_latent,
               double *save_sensible,
               double *save_Qnet,
	       double *save_Trad,
	       double *save_refreeze_energy)
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
  if (atmos->air_temp > 0.0)
    SnowFallCC = 0.0;
  else
    SnowFallCC = CH_ICE * SnowFall * atmos->air_temp;
  
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
         z2, displacement, Z0, wind, atmos->net_short,
         atmos->longwave, atmos->density, Le, atmos->air_temp,
         atmos->pressure * 1000., atmos->vpd * 1000., atmos->vp * 1000.,
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
      RefrozenWater = RefreezeEnergy/(Lf * WaterDensity(snow->surf_temp)) 
          * delta_t * SECPHOUR; 
      if (RefrozenWater > snow->surf_water) {
        RefrozenWater = snow->surf_water;
        RefreezeEnergy = RefrozenWater * Lf * WaterDensity(snow->surf_temp)/
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
      SnowMelt = fabs(RefreezeEnergy)/(Lf * WaterDensity(snow->surf_temp)) * 
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

    snow->surf_temp = root_brent((double)(snow->surf_temp-DELTAT), (double)0.0,
        SnowPackEnergyBalance, delta_t, aero_resist, z2, 
        displacement, Z0, wind, atmos->net_short, atmos->longwave,
        atmos->density, Le, atmos->air_temp, atmos->pressure * 1000.,
        atmos->vpd * 1000., atmos->vp * 1000., RainFall, SurfaceSwq,
        snow->surf_water, OldTSurf, &RefreezeEnergy, 
        &vapor_flux, &advection, &deltaCC, grnd_surf_temp, snow->depth,
        snow->density,surf_atten,&grnd_flux,&latent_heat,&sensible_heat,
				 save_Trad);

    Qnet = CalcSnowPackEnergyBalance(snow->surf_temp, delta_t, aero_resist,
           z2, displacement, Z0, wind, atmos->net_short,
           atmos->longwave, atmos->density, Le, atmos->air_temp,
           atmos->pressure * 1000., atmos->vpd * 1000., atmos->vp * 1000.,
           RainFall, SurfaceSwq, snow->surf_water, OldTSurf,
           &RefreezeEnergy, &vapor_flux, &advection, &deltaCC, grnd_surf_temp,
           snow->depth, snow->density,surf_atten,&grnd_flux,&latent_heat,
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
                       * WaterDensity(snow->surf_temp)/(delta_t * SECPHOUR);
    
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
    atmos->melt    = snow->surf_water - MaxLiquidWater;
    snow->surf_water = MaxLiquidWater;
  }
  else
    atmos->melt = 0.0;
  
  /* Refreeze liquid water in the pack.                                   
     variable 'RefreezeEnergy' is the heat released to the snow pack
     if all liquid water were refrozen.                                   
     if RefreezeEnergy < PackCC then all water IS refrozen           
     PackCC always <=0.0 

     WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does 
     not involve energy transported to the pixel.  Instead heat from the snow 
     pack is used to refreeze water */
  
  snow->pack_water += atmos->melt; /* add surface layer outflow to pack 
                                      liquid water*/
  RefreezeEnergy = snow->pack_water * Lf * WaterDensity(snow->pack_temp);

  /* calculate energy released to freeze*/
  
  if (PackCC < -RefreezeEnergy) { /* cold content not fully depleted*/
    PackSwq   += snow->pack_water;    /* refreeze all water and update*/
    Ice       += snow->pack_water;
    snow->pack_water = 0.0;
    if (PackSwq > 0.0)
      snow->pack_temp = (PackCC + RefreezeEnergy) / (CH_ICE * PackSwq);
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
    DeltaPackSwq = -PackCC/(Lf * WaterDensity(snow->pack_temp)); 
    snow->pack_water -= DeltaPackSwq;
    PackSwq += DeltaPackSwq;
    Ice += DeltaPackSwq;
  }
  
  /* Update the liquid water content of the pack */
  
  MaxLiquidWater = LIQUID_WATER_CAPACITY * PackSwq;
  if (snow->pack_water > MaxLiquidWater) {
    atmos->melt    = snow->pack_water - MaxLiquidWater;
    snow->pack_water = MaxLiquidWater;
  }
  else
    atmos->melt = 0.0;

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
    - atmos->melt + snow->vapor_flux; 
  
/*  printf("%d %d %g\n", y, x, MassBalanceError);*/

  atmos->melt *= 1000.; /* converts back to mm */
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
