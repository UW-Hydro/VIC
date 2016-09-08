/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate lake ice energy balance
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
IceEnergyBalance(double  TSurf,
                 va_list ap)
{
    extern parameters_struct param;

    /* start of list of arguments in variable argument list */

    double  Dt;                  /* Model time step (seconds) */
    double  Ra;                  /* Aerodynamic resistance (s/m) */
    double *Ra_used;             /* Aerodynamic resistance (s/m) after stability correction */
    double  Z;                   /* Reference height (m) */
    double  Z0;                  /* surface roughness height (m) */
    double  Wind;                /* Wind speed (m/s) */
    double  ShortRad;            /* Net incident shortwave radiation (W/m2) */
    double  LongRadIn;           /* Incoming longwave radiation (W/m2) */
    double  AirDens;             /* Density of air (kg/m3) */
    double  Lv;                  /* Latent heat of vaporization (J/kg3) */
    double  Tair;                /* Air temperature (C) */
    double  Press;               /* Air pressure (Pa) */
    double  Vpd;                 /* Vapor pressure deficit (Pa) */
    double  EactAir;             /* Actual vapor pressure of air (Pa) */
    double  Rain;                /* Rain fall (m/timestep) */
    double  SurfaceLiquidWater;  /* Liquid water in the surface layer (m) */
    double *RefreezeEnergy;      /* Refreeze energy (W/m2) */
    double *vapor_flux;          /* Total mass flux of water vapor to or from
                                    snow (m/timestep) */
    double *blowing_flux;        /* Mass flux of water vapor to or from
                                    blowing snow (m/timestep) */
    double *surface_flux;        /* Mass flux of water vapor to or from
                                    snow pack (m/timestep) */
    double *AdvectedEnergy;       /* Energy advected by precipitation (W/m2) */
    double  Tfreeze;
    double  AvgCond;
    double  SWconducted;
    double *qf;         /* Ground Heat Flux (W/m2) */
    double *LatentHeat;         /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;      /* Latent heat exchange at surface (W/m2) due to sublimation */
    double *SensibleHeat;       /* Sensible heat exchange at surface (W/m2) */
    double *LongRadOut;

    /* end of list of arguments in variable argument list */

    double Density;              /* Density of water/ice at TMean (kg/m3) */
    double NetRad;                      /* Net radiation exchange at surface (W/m2) */
    double RestTerm;            /* Rest term in surface energy balance
                                   (W/m2) */
    double TMean;               /* Average temperature for time step (C) */
    double qnull;
    double VaporMassFlux;        /* Total mass flux of water vapor to or from
                                    snow (kg/m2s) */
    double BlowingMassFlux;      /* Mass flux of water vapor to or from
                                    blowing snow (kg/m2s) */
    double SurfaceMassFlux;      /* Mass flux of water vapor to or from
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
    Dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    Ra_used = (double *) va_arg(ap, double *);
    Z = (double) va_arg(ap, double);
    Z0 = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);
    ShortRad = (double) va_arg(ap, double);
    LongRadIn = (double) va_arg(ap, double);
    AirDens = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Tair = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    SurfaceLiquidWater = (double) va_arg(ap, double);
    RefreezeEnergy = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);
    blowing_flux = (double *) va_arg(ap, double *);
    surface_flux = (double *) va_arg(ap, double *);
    AdvectedEnergy = (double *) va_arg(ap, double *);
    Tfreeze = (double) va_arg(ap, double);
    AvgCond = (double) va_arg(ap, double);
    SWconducted = (double) va_arg(ap, double);
    qf = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    LongRadOut = (double *) va_arg(ap, double *);

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
        *Ra_used = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0);
    }
    else {
        *Ra_used = param.HUGE_RESIST;
    }

    /* Calculate longwave exchange and net radiation */

    *LongRadOut = LongRadIn - calc_outgoing_longwave(TMean + CONST_TKFRZ,
                                                     param.EMISS_SNOW);
    NetRad = ShortRad + *LongRadOut;

    /* Calculate the sensible heat flux */

    *SensibleHeat = calc_sensible_heat(AirDens, Tair, TMean, *Ra_used);

    /* Calculate the mass flux of ice to or from the surface layer */

    /* Calculate sublimation terms and latent heat flux */

    /* blowing_flux was calculated outside of the root_brent iteration */
    BlowingMassFlux = *blowing_flux * Density / Dt;

    latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra, TMean, Vpd,
                          LatentHeat, LatentHeatSub, &VaporMassFlux,
                          &BlowingMassFlux,
                          &SurfaceMassFlux);

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * Dt / Density;
    *surface_flux = SurfaceMassFlux * Dt / Density;

    /* Calculate advected heat flux from rain */

    // Temporary fix for lake model.
    *AdvectedEnergy = (CONST_CPFW * CONST_RHOFW * Tair * Rain) / Dt;

    /* Calculate change in cold content */
    /* No change in cold content in lake model */

    /* Changes for lake model start here. Actually, equals qo-Io (P&H eq. 7)*/
    /* because Io (SWnet) is included in NetRad below. */
    qnull = (1 / AvgCond) * (Tfreeze - TMean + SWconducted);
    *qf = qnull;

    /* Changes for lake model end here. */

    /* Calculate net energy exchange at the snow surface */

    RestTerm =
        (NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub +
         *AdvectedEnergy +
         qnull);

    *RefreezeEnergy = (SurfaceLiquidWater * CONST_LATICE * Density) / Dt;

    /* Melting, or partially refreeze surface water. */
    if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
        *RefreezeEnergy = -RestTerm; /* available energy input over cold content
                                        used to melt, i.e. Qrf is negative value
                                        (energy out of pack)*/
        RestTerm = 0.0;
    }
    else {                       /* Pack is getting colder. */
        RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
    }

    return RestTerm;
}
