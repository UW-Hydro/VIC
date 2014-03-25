/*
 * SUMMARY:      IceEnergyBalance.c - Calculate lake ice energy balance
 * USAGE:        Part of the lake algorithm
 *
 * AUTHOR:       Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     March 16, 2001
 * LAST-MOD: Mon Apr 21 15:51:12 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate ice energy balance
 * DESCRIP-END.
 * FUNCTIONS:    IceEnergyBalance()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: IceEnergyBalance()

  Purpose      : Calculate the surface energy balance for the snow pack

  Required     :
    double TSurf           - new estimate of effective surface temperature
    va_list ap            - Argument list initialized by va_start().  For
                            elements of list and order, see beginning of
                            routine

  Returns      :
    double RestTerm        - Rest term in the energy balance

  Modifies     : 
    double *RefreezeEnergy - Refreeze energy (W/m2) 
    double *VaporMassFlux  - Mass flux of water vapor to or from the
                            intercepted snow (m/s)

  Comments     :
    Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
                science, Addisson Wesley, Inc., Reading, etc., 1990.

  Modifications:
  16-Jul-04 Renamed VaporMassFlux to vapor_flux to denote fact that its
	    units are m/timestep rather than kg/m2s.  Created new variable
	    VaporMassFlux with units of kg/m2s.  After VaporMassFlux is
	    computed, vapor_flux is derived from it by unit conversion.
	    vapor_flux is the variable that is passed in/out of this
	    function.							TJB
  24-Aug-04 Added logic to handle blowing_flux and vapor_flux.		TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in flux
	    calculations.						TJB
  04-Oct-04 Merged with Laura Bowling's updated lake model code.  Now
	    blowing snow sublimation is calculated for lakes.		TJB
  2006-Sep-23 Replaced redundant STEFAN constant with STEFAN_B.  TJB
  2006-Nov-07 Removed LAKE_MODEL option. TJB

*****************************************************************************/
double IceEnergyBalance(double TSurf, va_list ap)
{

  extern option_struct options;

  const char *Routine = "IceEnergyBalance";

  /* start of list of arguments in variable argument list */

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
  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double DeltaColdContent;       /* Change in cold content (W/m2) */
  double Tfreeze;     
  double AvgCond;
  double SWconducted;
  double SnowDepth;  
  double SnowDensity; 
  double SurfAttenuation; 
  double *qf;		/* Ground Heat Flux (W/m2) */
  double *LatentHeat;		/* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;	/* Latent heat exchange at surface (W/m2) due to sublimation */
  double *SensibleHeat;		/* Sensible heat exchange at surface (W/m2) */
  double *LongRadOut;

  /* end of list of arguments in variable argument list */

  double Density;                /* Density of water/ice at TMean (kg/m3) */
  double EsSnow;                 /* saturated vapor pressure in the snow pack
                                   (Pa)  */  
  
  double Ls;                     /* Latent heat of sublimation (J/kg) */
  double NetRad;			/* Net radiation exchange at surface (W/m2) */
  double RestTerm;		/* Rest term in surface energy balance
				   (W/m2) */
  double  TMean;                /* Average temperature for time step (C) */
  double qnull;
  double VaporMassFlux;          /* Total mass flux of water vapor to or from
                                    snow (kg/m2s) */
  double BlowingMassFlux;        /* Mass flux of water vapor to or from
                                    blowing snow (kg/m2s) */
  double SurfaceMassFlux;        /* Mass flux of water vapor to or from
                                    snow pack (kg/m2s) */

  /* Assign the elements of the array to the appropriate variables.  The list
     is traversed as if the elements are doubles, because:

     In the variable-length part of variable-length argument lists, the old
     ``default argument promotions'' apply: arguments of type double are
     always promoted (widened) to type double, and types char and short int
     are promoted to int. Therefore, it is never correct to invoke
     va_arg(argp, double); instead you should always use va_arg(argp,
     double). 

     (quoted from the comp.lang.c FAQ list)
     */
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
  DeltaColdContent   = (double) va_arg(ap, double );
  Tfreeze            = (double) va_arg(ap, double);
  AvgCond            = (double) va_arg(ap, double);
  SWconducted        = (double) va_arg(ap, double);
  SnowDepth          = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfAttenuation    = (double) va_arg(ap, double);
  qf                 = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  LatentHeatSub      = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  LongRadOut         = (double *) va_arg(ap, double *);
  
  /* Calculate active temp for energy balance as average of old and new  */
  
/*   TMean = 0.5 * (OldTSurf + TSurf); */
  TMean = TSurf;
  Density = RHO_W;
  
  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */                        
  
  /* Apply the stability correction to the aerodynamic resistance 
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */


  if (Wind > 0.0)
    *Ra_used = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0);
    /* *Ra_used = Ra / StabilityCorrection(2.f, 0.f, TMean, Tair, Wind, Z0);*/
  else
    *Ra_used = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */
 
  *LongRadOut = LongRadIn - STEFAN_B * (TMean+273.15) * (TMean+273.15) 
    * (TMean+273.15) * (TMean+273.15);
  NetRad = ShortRad + *LongRadOut;
  
  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * CP_PM * (Tair - TMean) / *Ra_used;

  /* Calculate the mass flux of ice to or from the surface layer */
 
  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  EsSnow = svp(TMean) /* * 1000. */;

/*   EsSnow = 610.78 * exp((double)((17.269 * TMean) / (237.3 + TMean))); */

/*   if (TMean < 0.0) */
/*     EsSnow *= 1.0 + .00972 * TMean + .000042 * pow((double)TMean,(double)2.0); */
  
  /* Calculate sublimation terms and latent heat flux */

  /* blowing_flux was calculated outside of the root_brent iteration */
  BlowingMassFlux = *blowing_flux * Density / (Dt * SECPHOUR);

  latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra, TMean, Vpd,
                        LatentHeat, LatentHeatSub, &VaporMassFlux, &BlowingMassFlux,
                        &SurfaceMassFlux);

  /* Convert sublimation terms from kg/m2s to m/timestep */
  *vapor_flux = VaporMassFlux * Dt * SECPHOUR / Density;
  *surface_flux = SurfaceMassFlux * Dt * SECPHOUR / Density;

  /* Calculate advected heat flux from rain 
     WORK IN PROGRESS:  Should the following read (Tair - Tsurf) ?? */
  
  // Temporary fix for lake model.
  *AdvectedEnergy = (CH_WATER * Tair * Rain) / (Dt*SECPHOUR);
  //*AdvectedEnergy = 0.0;
  
  /* Calculate change in cold content */
  /* No change in cold content in lake model */

  /* Changes for lake model start here. Actually, equals qo-Io (P&H eq. 7)*/
    /* because Io (SWnet) is included in NetRad below. */
  qnull = (1/AvgCond)*(Tfreeze - TMean + SWconducted);
  *qf   = qnull;

  /* Changes for lake model end here. */
  
  /* Calculate net energy exchange at the snow surface */
  
  RestTerm = ( NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy 
	       + qnull );

  *RefreezeEnergy = (SurfaceLiquidWater * Lf * Density)/(Dt * SECPHOUR);

  /* Melting, or partially refreeze surface water. */
  if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
    *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
                                    used to melt, i.e. Qrf is negative value
                                    (energy out of pack)*/ 
    RestTerm = 0.0;
  }
  else {                         /* Pack is getting colder. */
    RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
  }
  
  return RestTerm;
}
