/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate the aerodynamic resistances.
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
 * @brief    Calculate the aerodynamic resistance for each vegetation layer,
 *           and the wind 2m above the layer boundary.  In case of an overstory,
 *           also calculate the wind in the overstory. The values are
 *           normalized based on a reference height wind speed, Uref, of 1 m/s.
 *           To get wind speeds and aerodynamic resistances for other values of
 *           Uref, you need to multiply the here calculated wind speeds by Uref
 *           and divide the here calculated aerodynamic resistances by Uref.
 *****************************************************************************/
int
CalcAerodynamic(bool    OverStory,          /* overstory flag */
                double  Height,             /* vegetation height */
                double  Trunk,              /* trunk ratio parameter */
                double  Z0_SNOW,            /* snow roughness */
                double  Z0_SOIL,            /* soil roughness */
                double  n,                  /* wind attenuation parameter */
                double *Ra,                 /* aerodynamic resistances */
                double *U,                  /* adjusted wind speed */
                double *displacement,       /* vegetation displacement */
                double *ref_height,         /* vegetation reference height */
                double *roughness)          /* vegetation roughness */
{
    extern parameters_struct param;

    double                   d_Lower;
    double                   d_Upper;
    double                   K2;
    double                   Uh;
    double                   Ut;
    double                   Uw;
    double                   Z0_Lower;
    double                   Z0_Upper;
    double                   Zt;
    double                   Zw;
    double                   tmp_wind;

    tmp_wind = U[0];

    K2 = CONST_KARMAN * CONST_KARMAN;

    /* No OverStory, thus maximum one soil layer */

    if (!OverStory) {
        /* vegetation cover */
        Z0_Lower = roughness[0];
        d_Lower = displacement[0];

        /* No snow */
        U[0] = log((2. + Z0_Lower) / Z0_Lower) / log(
            (ref_height[0] - d_Lower) / Z0_Lower);

        /* Old VIC, may not match DHSVM */
        Ra[0] = log((2. + (1.0 / 0.63 - 1.0) * d_Lower) / Z0_Lower) *
                log((2. +
                     (1.0 / 0.63 - 1.0) * d_Lower) / (0.1 * Z0_Lower)) / K2;

        /* Copy bare parameters into canopy top parameters */
        U[1] = U[0];
        Ra[1] = Ra[0];

        /** Set aerodynamic resistance terms for canopy */
        /* not currently used */
        ref_height[1] = ref_height[0];
        roughness[1] = roughness[0];
        displacement[1] = displacement[0];

        /* Snow */
        U[2] = log((2. + Z0_SNOW) / Z0_SNOW) / log(ref_height[0] / Z0_SNOW);
        Ra[2] =
            log((2. + Z0_SNOW) / Z0_SNOW) * log(ref_height[0] / Z0_SNOW) / K2;
        /** Set aerodynamic resistance terms for snow */
        ref_height[2] = 2. + Z0_SNOW;
        roughness[2] = Z0_SNOW;
        displacement[2] = 0.;
    }
    /* Overstory present, one or two vegetation layers possible */
    else {
        Z0_Upper = roughness[0];
        d_Upper = displacement[0];

        Z0_Lower = Z0_SOIL;
        d_Lower = 0;

        Zw = 1.5 * Height - 0.5 * d_Upper;
        Zt = Trunk * Height;

        if (Zt < (Z0_Lower + d_Lower)) {
            log_err("Trunk space height below \"center\" of lower boundary");
        }

        /* Resistance for overstory */
        Ra[1] = log((ref_height[0] - d_Upper) / Z0_Upper) / K2 *
                (Height / (n * (Zw - d_Upper)) *
                 (exp(n * (1 - (d_Upper + Z0_Upper) / Height)) - 1) +
                 (Zw - Height) / (Zw - d_Upper) +
                 log((ref_height[0] - d_Upper) / (Zw - d_Upper)));

        /* Wind at different levels in the profile */
        Uw = log((Zw - d_Upper) / Z0_Upper) / log(
            (ref_height[0] - d_Upper) / Z0_Upper);
        Uh = Uw - (1 - (Height - d_Upper) / (Zw - d_Upper)) /
             log((ref_height[0] - d_Upper) / Z0_Upper);
        U[1] = Uh * exp(n * ((Z0_Upper + d_Upper) / Height - 1.));
        Ut = Uh * exp(n * (Zt / Height - 1.));

        /* resistance at the lower boundary */
        U[0] = log((2. + Z0_Upper) / Z0_Upper) / log(
            (ref_height[0] - d_Upper) / Z0_Upper);
        Ra[0] = log((2. + (1.0 / 0.63 - 1.0) * d_Upper) / Z0_Upper) *
                log((2. +
                     (1.0 / 0.63 - 1.0) * d_Upper) / (0.1 * Z0_Upper)) / K2;

        /* Snow */

        /* case 1: the wind profile to a height of 2m above the lower boundary is
           entirely logarithmic */
        if (Zt > (2. + Z0_SNOW)) {
            U[2] = Ut * log((2. + Z0_SNOW) / Z0_SNOW) / log(Zt / Z0_SNOW);
            Ra[2] =
                log((2. + Z0_SNOW) / Z0_SNOW) * log(Zt / Z0_SNOW) / (K2 * Ut);
        }

        /* case 2: the wind profile to a height of 2m above the lower boundary
           is part logarithmic and part exponential, but the top of the overstory
           is more than 2 m above the lower boundary */
        else if (Height > (2. + Z0_SNOW)) {
            U[2] = Uh * exp(n * ((2. + Z0_SNOW) / Height - 1.));
            Ra[2] = log(Zt / Z0_SNOW) * log(Zt / Z0_SNOW) /
                    (K2 * Ut) +
                    Height *
                    log((ref_height[0] -
                         d_Upper) / Z0_Upper) / (n * K2 * (Zw - d_Upper)) *
                    (exp(n *
                         (1 - Zt /
                          Height)) - exp(n * (1 - (Z0_SNOW + 2.) / Height)));
        }

        /* case 3: the top of the overstory is less than 2 m above the lower
           boundary.  The wind profile above the lower boundary is part
           logarithmic and part exponential, but only extends to the top of the
           overstory */
        else {
            U[2] = Uh;
            Ra[2] = log(Zt / Z0_SNOW) * log(Zt / Z0_SNOW) /
                    (K2 * Ut) +
                    Height *
                    log((ref_height[0] -
                         d_Upper) / Z0_Upper) / (n * K2 * (Zw - d_Upper)) *
                    (exp(n * (1 - Zt / Height)) - 1);
            log_warn("Top of overstory is less than 2 meters above the lower "
                     "boundary");
        }

        /** Set aerodynamic resistance terms for canopy */
        /* not currently used */
        ref_height[1] = ref_height[0];
        roughness[1] = roughness[0];
        displacement[1] = displacement[0];
        ref_height[0] = 2.;
        roughness[0] = Z0_Lower;
        displacement[0] = d_Lower;

        /** Set aerodynamic resistance terms for snow */
        ref_height[2] = 2. + Z0_SNOW;
        roughness[2] = Z0_SNOW;
        displacement[2] = 0.;
    }

    if (tmp_wind > 0.) {
        U[0] *= tmp_wind;
        Ra[0] /= tmp_wind;
        if (U[1] != -999) {
            U[1] *= tmp_wind;
            Ra[1] /= tmp_wind;
        }
        if (U[2] != -999) {
            U[2] *= tmp_wind;
            Ra[2] /= tmp_wind;
        }
    }
    else {
        U[0] *= tmp_wind;
        Ra[0] = param.HUGE_RESIST;
        if (U[1] != -999) {
            U[1] *= tmp_wind;
        }
        Ra[1] = param.HUGE_RESIST;
        if (U[2] != -999) {
            U[2] *= tmp_wind;
        }
        Ra[2] = param.HUGE_RESIST;
    }
    return (0);
}
