/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow pack energy balance
 *
 * Based on the SnowPackEnergyBalance function in DHSVM
 * Reference: Bras, R.A., Hydrology, an introduction to hydrologic science,
 * Addison Wesley, Inc., Reading, etc., 1990.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate the surface energy balance for the snow pack.
 *****************************************************************************/
double
SnowPackEnergyBalance(double  TSurf,
                      va_list ap)
{
    extern option_struct     options;
    extern parameters_struct param;

    /* Define Variable Argument List */

    /* General Model Parameters */
    double  Dt;                   /* Model time step (sec) */
    double  Ra;                   /* Aerodynamic resistance (s/m) */
    double *Ra_used;              /* Aerodynamic resistance (s/m) after stability correction */

    /* Vegetation Parameters */
    double  Z;                    /* Reference height (m) */
    double *Z0;                    /* surface roughness height (m) */

    /* Atmospheric Forcing Variables */
    double  AirDens;              /* Density of air (kg/m3) */
    double  EactAir;              /* Actual vapor pressure of air (Pa) */
    double  LongSnowIn;            /* Incoming longwave radiation (W/m2) */
    double  Lv;                   /* Latent heat of vaporization (J/kg3) */
    double  Press;                /* Air pressure (Pa) */
    double  Rain;                 /* Rain fall (m/timestep) */
    double  NetShortUnder;        /* Net incident shortwave radiation
                                     (W/m2) */
    double  Vpd;                  /* Vapor pressure deficit (Pa) */
    double  Wind;                 /* Wind speed (m/s) */

    /* Snowpack Variables */
    double  OldTSurf;             /* Surface temperature during previous time
                                     step */
    double  SnowCoverFract;       /* Fraction of area covered by snow */
    double  SnowDepth;            /* Depth of snowpack (m) */
    double  SnowDensity;          /* Density of snowpack (kg/m^3) */
    double  SurfaceLiquidWater;   /* Liquid water in the surface layer (m) */
    double  SweSurfaceLayer;      /* Snow water equivalent in surface layer
                                     (m) */

    /* Energy Balance Components */
    double  Tair;                 /* Canopy air / Air temperature (C) */
    double  TGrnd;                /* Ground surface temperature (C) */

    double *AdvectedEnergy;       /* Energy advected by precipitation (W/m2) */
    double *AdvectedSensibleHeat; /* Sensible heat advected from snow-free
                                     area into snow covered area (W/m^2) */
    double *DeltaColdContent;     /* Change in cold content of surface
                                     layer (W/m2) */
    double *GroundFlux;           /* Ground Heat Flux (W/m2) */
    double *LatentHeat;           /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub; /* Latent heat of sublimation exchange at
                                      surface (W/m2) */
    double *NetLongUnder;         /* Net longwave radiation at snowpack
                                     surface (W/m^2) */
    double *RefreezeEnergy;       /* Refreeze energy (W/m2) */
    double *SensibleHeat;         /* Sensible heat exchange at surface
                                     (W/m2) */
    double *vapor_flux;           /* Mass flux of water vapor to or from the
                                     intercepted snow (m/timestep) */
    double *blowing_flux;         /* Mass flux of water vapor from blowing snow. (m/timestep) */
    double *surface_flux;         /* Mass flux of water vapor from pack snow. (m/timestep) */

    /* Internal Routine Variables */

    double Density;               /* Density of water/ice at TMean (kg/m3) */
    double NetRad;                /* Net radiation exchange at surface
                                     (W/m2) */
    double RestTerm;              /* Rest term in surface energy balance
                                     (W/m2) */
    double TMean;                 /* Average temperature for time step (C) */
    double Tmp;
    double VaporMassFlux;         /* Mass flux of water vapor to or from the
                                     intercepted snow (kg/m2s) */
    double BlowingMassFlux;       /* Mass flux of water vapor from blowing snow. (kg/m2s) */
    double SurfaceMassFlux;       /* Mass flux of water vapor from pack snow. (kg/m2s) */

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
    Dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    Ra_used = (double *) va_arg(ap, double *);

    /* Vegetation Parameters */
    Z = (double) va_arg(ap, double);
    Z0 = (double *) va_arg(ap, double *);

    /* Atmospheric Forcing Variables */
    AirDens = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    LongSnowIn = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    NetShortUnder = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);

    /* Snowpack Variables */
    OldTSurf = (double) va_arg(ap, double);
    SnowCoverFract = (double) va_arg(ap, double);
    SnowDepth = (double) va_arg(ap, double);
    SnowDensity = (double) va_arg(ap, double);
    SurfaceLiquidWater = (double) va_arg(ap, double);
    SweSurfaceLayer = (double) va_arg(ap, double);

    /* Energy Balance Components */
    Tair = (double) va_arg(ap, double);
    TGrnd = (double) va_arg(ap, double);

    AdvectedEnergy = (double *) va_arg(ap, double *);
    AdvectedSensibleHeat = (double *)va_arg(ap, double *);
    DeltaColdContent = (double *) va_arg(ap, double *);
    GroundFlux = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    NetLongUnder = (double *) va_arg(ap, double *);
    RefreezeEnergy = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);
    blowing_flux = (double *) va_arg(ap, double *);
    surface_flux = (double *) va_arg(ap, double *);

    /* Calculate active temp for energy balance as average of old and new  */

    TMean = TSurf;
    Density = CONST_RHOFW;

    /* Correct aerodynamic conductance for stable conditions
       Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
       velocity (vel_2m) is expected to be in m/sec */

    /* Apply the stability correction to the aerodynamic resistance
       NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
       think that it is more correct to calculate ALL fluxes at the same
       reference level */


    if (Wind > 0.0) {
        Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0[2]);
    }
    else {
        Ra_used[0] = param.HUGE_RESIST;
    }

    /* Calculate longwave exchange and net radiation */

    Tmp = TMean + CONST_TKFRZ;
    (*NetLongUnder) = LongSnowIn - calc_outgoing_longwave(Tmp,
                                                          param.EMISS_SNOW);
    NetRad = NetShortUnder + (*NetLongUnder);

    /* Calculate the sensible heat flux */

    *SensibleHeat = calc_sensible_heat(AirDens, Tair, TMean, Ra_used[0]);

    if (options.SPATIAL_SNOW) {
        /* Add in Sensible heat flux turbulent exchange from surrounding
           snow free patches - if present */
        if (SnowCoverFract > 0.) {
            *(AdvectedSensibleHeat) = advected_sensible_heat(SnowCoverFract,
                                                             AirDens, Tair,
                                                             TGrnd,
                                                             Ra_used[0]);
        }
        else {
            (*AdvectedSensibleHeat) = 0.;
        }
    }
    else {
        (*AdvectedSensibleHeat) = 0.;
    }

    /* Convert sublimation terms from m/timestep to kg/m2s */
    VaporMassFlux = *vapor_flux * Density / Dt;
    BlowingMassFlux = *blowing_flux * Density / Dt;
    SurfaceMassFlux = *surface_flux * Density / Dt;

    /* Calculate the mass flux of ice to or from the surface layer */

    /* Calculate the saturated vapor pressure in the snow pack,
       (Equation 3.32, Bras 1990) */

    latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra_used[0], TMean, Vpd,
                          LatentHeat, LatentHeatSub, &VaporMassFlux,
                          &BlowingMassFlux,
                          &SurfaceMassFlux);

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * Dt / Density;
    *blowing_flux = BlowingMassFlux * Dt / Density;
    *surface_flux = SurfaceMassFlux * Dt / Density;

    /* Calculate advected heat flux from rain
       Equation 7.3.12 from H.B.H. for rain falling on melting snowpack */

    if (TMean == 0.) {
        *AdvectedEnergy = (CONST_CPFW * CONST_RHOFW * (Tair) * Rain) / Dt;
    }
    else {
        *AdvectedEnergy = 0.;
    }

    /* Calculate change in cold content */
    *DeltaColdContent = CONST_VCPICE_WQ * SweSurfaceLayer *
                        (TSurf - OldTSurf) / Dt;

    /* Calculate Ground Heat Flux */
    if (SnowDepth > param.SNOW_DEPTH_THRES) {
        *GroundFlux = param.SNOW_CONDUCT * pow(SnowDensity, 2.) *
                      (TGrnd - TMean) / SnowDepth / Dt;
    }
    else {
        *GroundFlux = 0;
    }
    *DeltaColdContent -= *GroundFlux;

    /* Calculate energy balance error at the snowpack surface */
    RestTerm = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub +
               *AdvectedEnergy + *GroundFlux - *DeltaColdContent +
               *AdvectedSensibleHeat;

    *RefreezeEnergy = (SurfaceLiquidWater * CONST_LATICE * Density) / Dt;

    if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
        *RefreezeEnergy = -RestTerm; /* available energy input over cold content
                                        used to melt, i.e. Qrf is negative value
                                        (energy out of pack)*/
        RestTerm = 0.0;
    }
    else {
        RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
    }

    return RestTerm;
}
