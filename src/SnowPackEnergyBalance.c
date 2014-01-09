/*
 * SUMMARY:      SnowPackEnergyBalance.c - Calculate snow pack energy balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 09:09:29
 * LAST-MOD: Mon Apr 21 16:39:04 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate snow pack energy balance
 * DESCRIP-END.
 * FUNCTIONS:    SnowPackEnergyBalance()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: SnowPackEnergyBalance()

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
    double *vapor_flux     - Mass flux of water vapor to or from the
                             intercepted snow (m/timestep)
    double *blowing_flux   - Mass flux of water vapor from blowing snow (m/timestep)
    double *surface_flux   - Mass flux of water vapor from pack snow (m/timestep)

  Comments     :
    Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
                science, Addisson Wesley, Inc., Reading, etc., 1990.

  Modifications:
  10-6-2000 modified to handle partial snow cover including the
    advection of sensible heat from bare patches to the edges
    of the remaining snow cover.                               KAC
  11-18-02 Modified to include the effects of blowing snow on the
           accumulation and ablation of the snowpack.          LCB
  04-21-03 Removed constant variable declarations.  Those still
           required by the VIC model are now defined in 
           vicNl_def.h.                                        KAC
  16-Jul-04 Renamed VaporMassFlux, BlowingMassFlux, and SurfaceMassFlux
	    to vapor_flux, blowing_flux, and surface_flux, respectively,
	    to denote fact that their units are m/timestep rather than
	    kg/m2s.  Created new variables VaporMassFlux, BlowingMassFlux,
	    and SurfaceMassFlux with units of kg/m2s.  The addresses of
	    the *MassFlux variables are passed to latent_heat_from_snow()
	    where values for the variables are computed.  After these
	    values are computed, vapor_flux, blowing_flux and surface_flux
	    are derived from them by unit conversion.  vapor_flux,
	    blowing_flux, and surface_flux are the variables that are
	    passed in/out of this function.				TJB
  16-Jul-04 Changed the type of the last few variables (lag_one, iveg,
	    etc) in the va_list to be double.  For some reason, passing
	    them as float or int caused them to become garbage.  This may
	    have to do with the fact that they followed variables of type
	    (double *) in va_list, which may have caused memory alignment
	    problems.							TJB
  05-Aug-04 Removed lag_one, sigma_slope, fetch, iveg, Nveg, month, overstory,
	    LastSnow, and SnowDepth from argument list, since they were only
	    needed to pass to latent_heat_from_snow(), which no longer needs
	    them.							TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in flux
	    calculations.						TJB
  2006-Sep-23 Replaced redundant STEFAN_B constant with STEFAN_B_B.	TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-Sep-19 Added Added ground flux computation consistent with 4.0.6.	TJB
  2013-Dec-27 Moved SPATIAL_SNOW from compile-time to run-time options.	TJB

*****************************************************************************/
double SnowPackEnergyBalance(double TSurf, va_list ap)
{

  extern option_struct options;

  const char *Routine = "SnowPackEnergyBalance";

  /* Define Variable Argument List */

  /* General Model Parameters */
  double Dt;                      /* Model time step (sec) */
  double Ra;                      /* Aerodynamic resistance (s/m) */
  double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction */

  /* Vegetation Parameters */
  double Displacement;            /* Displacement height (m) */
  double Z;                       /* Reference height (m) */
  double *Z0;                      /* surface roughness height (m) */

  /* Atmospheric Forcing Variables */
  double AirDens;                 /* Density of air (kg/m3) */
  double EactAir;                 /* Actual vapor pressure of air (Pa) */
  double LongSnowIn;               /* Incoming longwave radiation (W/m2) */
  double Lv;                      /* Latent heat of vaporization (J/kg3) */
  double Press;                   /* Air pressure (Pa) */
  double Rain;                    /* Rain fall (m/timestep) */
  double NetShortUnder;           /* Net incident shortwave radiation 
				     (W/m2) */ 
  double Vpd;			  /* Vapor pressure deficit (Pa) */
  double Wind;                    /* Wind speed (m/s) */

  /* Snowpack Variables */
  double OldTSurf;                /* Surface temperature during previous time
				     step */ 
  double SnowCoverFract;          /* Fraction of area covered by snow */
  double SnowDepth;               /* Depth of snowpack (m) */
  double SnowDensity;             /* Density of snowpack (kg/m^3) */
  double SurfaceLiquidWater;      /* Liquid water in the surface layer (m) */
  double SweSurfaceLayer;         /* Snow water equivalent in surface layer 
				     (m) */ 

  /* Energy Balance Components */
  double Tair;                    /* Canopy air / Air temperature (C) */
  double TGrnd;                   /* Ground surface temperature (C) */

  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *AdvectedSensibleHeat;	  /* Sensible heat advected from snow-free
				     area into snow covered area (W/m^2) */
  double *DeltaColdContent;       /* Change in cold content of surface 
				     layer (W/m2) */
  double *GroundFlux;		  /* Ground Heat Flux (W/m2) */
  double *LatentHeat;		  /* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;  /* Latent heat of sublimation exchange at 
				     surface (W/m2) */
  double *NetLongUnder;           /* Net longwave radiation at snowpack 
				     surface (W/m^2) */
  double *RefreezeEnergy;         /* Refreeze energy (W/m2) */
  double *SensibleHeat;		  /* Sensible heat exchange at surface 
				     (W/m2) */
  double *vapor_flux;             /* Mass flux of water vapor to or from the
				     intercepted snow (m/timestep) */
  double *blowing_flux;           /* Mass flux of water vapor from blowing snow. (m/timestep) */
  double *surface_flux;           /* Mass flux of water vapor from pack snow. (m/timestep) */

  /* Internal Routine Variables */

  double Density;                 /* Density of water/ice at TMean (kg/m3) */
  /* double LongRadOut; */		  /* long wave radiation emitted by surface
				     (W/m2) */
  double NetRad;		  /* Net radiation exchange at surface 
				     (W/m2) */
  double RestTerm;		  /* Rest term in surface energy balance
				     (W/m2) */
  double TMean;                   /* Average temperature for time step (C) */
  double Tmp;
  double VaporMassFlux;           /* Mass flux of water vapor to or from the
				     intercepted snow (kg/m2s) */
  double BlowingMassFlux;         /* Mass flux of water vapor from blowing snow. (kg/m2s) */
  double SurfaceMassFlux;         /* Mass flux of water vapor from pack snow. (kg/m2s) */

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

  /* General Model Parameters */
  Dt           = (double) va_arg(ap, double);
  Ra           = (double) va_arg(ap, double);
  Ra_used      = (double *) va_arg(ap, double *);

  /* Vegetation Parameters */
  Displacement = (double) va_arg(ap, double);
  Z            = (double) va_arg(ap, double);
  Z0           = (double *) va_arg(ap, double *);

  /* Atmospheric Forcing Variables */
  AirDens       = (double) va_arg(ap, double);
  EactAir       = (double) va_arg(ap, double);
  LongSnowIn    = (double) va_arg(ap, double);
  Lv            = (double) va_arg(ap, double);
  Press         = (double) va_arg(ap, double);
  Rain          = (double) va_arg(ap, double);
  NetShortUnder = (double) va_arg(ap, double);
  Vpd           = (double) va_arg(ap, double);
  Wind          = (double) va_arg(ap, double);

  /* Snowpack Variables */
  OldTSurf           = (double) va_arg(ap, double);
  SnowCoverFract     = (double) va_arg(ap, double);
  SnowDepth          = (double) va_arg(ap, double);
  SnowDensity        = (double) va_arg(ap, double);
  SurfaceLiquidWater = (double) va_arg(ap, double);
  SweSurfaceLayer    = (double) va_arg(ap, double);

  /* Energy Balance Components */
  Tair = (double) va_arg(ap, double);
  TGrnd   = (double) va_arg(ap, double);

  AdvectedEnergy        = (double *) va_arg(ap, double *);
  AdvectedSensibleHeat  = (double *)va_arg(ap, double *);
  DeltaColdContent      = (double *) va_arg(ap, double *);
  GroundFlux            = (double *) va_arg(ap, double *);
  LatentHeat            = (double *) va_arg(ap, double *);
  LatentHeatSub         = (double *) va_arg(ap, double *);
  NetLongUnder          = (double *) va_arg(ap, double *);
  RefreezeEnergy        = (double *) va_arg(ap, double *);
  SensibleHeat          = (double *) va_arg(ap, double *);
  vapor_flux            = (double *) va_arg(ap, double *);
  blowing_flux          = (double *) va_arg(ap, double *); 
  surface_flux          = (double *) va_arg(ap, double *);   

  /* Calculate active temp for energy balance as average of old and new  */
  
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
    Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0[2]); 
  else
    Ra_used[0] = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */

  Tmp = TMean + KELVIN;
  (*NetLongUnder) = LongSnowIn - STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  NetRad = NetShortUnder + (*NetLongUnder);
  
  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * Cp * (Tair - TMean) / Ra_used[0];

  if (options.SPATIAL_SNOW) {
    /* Add in Sensible heat flux turbulent exchange from surrounding 
       snow free patches - if present */
    if ( SnowCoverFract > 0 ) {
      *(AdvectedSensibleHeat) = advected_sensible_heat(SnowCoverFract, 
						       AirDens, Tair, TGrnd, 
						       Ra_used[0]);
    }
    else (*AdvectedSensibleHeat) = 0;
  }
  else {
    (*AdvectedSensibleHeat) = 0;
  }

  /* Convert sublimation terms from m/timestep to kg/m2s */
  VaporMassFlux = *vapor_flux * Density / Dt;
  BlowingMassFlux = *blowing_flux * Density / Dt;
  SurfaceMassFlux = *surface_flux * Density / Dt;

  /* Calculate the mass flux of ice to or from the surface layer */
 
  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra_used[0], TMean, Vpd,
			LatentHeat, LatentHeatSub, &VaporMassFlux, &BlowingMassFlux, 
			&SurfaceMassFlux);

  /* Convert sublimation terms from kg/m2s to m/timestep */
  *vapor_flux = VaporMassFlux * Dt / Density;
  *blowing_flux = BlowingMassFlux * Dt / Density;
  *surface_flux = SurfaceMassFlux * Dt / Density;
  
  /* Calculate advected heat flux from rain 
     Equation 7.3.12 from H.B.H. for rain falling on melting snowpack */
  
  if ( TMean == 0 )
    *AdvectedEnergy = (CH_WATER * (Tair) * Rain) / (Dt);
  else
    *AdvectedEnergy = 0.;
  
  /* Calculate change in cold content */
  *DeltaColdContent = CH_ICE * SweSurfaceLayer * (TSurf - OldTSurf) /
    (Dt);

  /* Calculate Ground Heat Flux */
  if(SnowDepth>0.) {
    *GroundFlux = 2.9302e-6 * SnowDensity * SnowDensity
        * (TGrnd - TMean) / SnowDepth / (Dt);
  }
  else *GroundFlux=0;
  *DeltaColdContent -= *GroundFlux;

  /* Calculate energy balance error at the snowpack surface */
  RestTerm = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub 
    + *AdvectedEnergy + *GroundFlux - *DeltaColdContent 
    + *AdvectedSensibleHeat;

  *RefreezeEnergy = (SurfaceLiquidWater * Lf * Density)/(Dt);

  if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
    *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
                                    used to melt, i.e. Qrf is negative value
                                    (energy out of pack)*/ 
    RestTerm = 0.0;
  }
  else {
    RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
  }
  
  return RestTerm;
}

