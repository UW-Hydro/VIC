/******************************************************************************
* @section DESCRIPTION
*
* Collection of snow utilities.
*
* @section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014  The Land Surface Hydrology Group, Department of Civil
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
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Compute the snow density based on swe and snow metamorphism.
*
* @note     DENS_SNTHRM = Algorithm is taken from SNTHERM89, adjusted for an
*               essentially single-layer model.
*           DENS_BRAS   = Original algorithm, originally based on a plot of
*               seasonal variation of typical snow densities found in Bras
*               (figure 6.10, p 258).  Because this equation was developed by
*               regression against data from Southern Manitoba, it therefore is
*               limited in applicability to other regions.
******************************************************************************/
double
snow_density(snow_data_struct *snow,
             double            new_snow,
             double            sswq,
             double            Tair,
             double            dt)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   density_new;
    double                   density;
    double                   depth;
    double                   swq;
    double                   CR; /* compaction rate */
    double                   dexpf;
    double                   ddz1; /* rate of settling of snowpack due to destructive metamorphism */
    double                   ddz2; /* rate of compaction of snowpack due to overburden */
    double                   Tavg; /* average snowpack temperature */
    double                   c3, c4; /* snow densification factors */
    double                   dm; /* upper snow density limit for the settling process */
    double                   Ps; /* snow load pressure (N/m^2) */
    double                   f; /* effective compaction coefficient */
    double                   overburden;
    double                   viscosity;
    double                   delta_depth;
    double                   depth_new;

    density = 0.0;

    /* Estimate density of new snow based on air temperature */
    if (new_snow > 0.) {
        density_new = new_snow_density(Tair);
    }
    else {
        density_new = 0.0;
    }

    /* Estimate average snowpack temperature */
    Tavg = snow->surf_temp + CONST_TKFRZ;

    if (options.SNOW_DENSITY == DENS_SNTHRM) {
        if (new_snow > 0.) {
            if (snow->depth > 0.0) {
                density = snow->density;
            }
            else {
                density = density_new;
            }
        }
        else {
            density = snow->density;
        }

        dexpf = exp(-param.SNOW_DENS_C1 * (CONST_TKFRZ - Tavg));

        /* Settling due to destructive metamorphism */
        if (new_snow > 0.0 && density_new > 0.0) {
            dm =
                (param.SNOW_DENS_DMLIMIT > 1.15 *
                 density_new) ? param.SNOW_DENS_DMLIMIT : 1.15 * density_new;
        }
        else {
            dm = param.SNOW_DENS_DMLIMIT;
        }

        if (density <= dm) {
            c3 = 1.0;
            c4 = 1.0;
        }
        else {
            c3 = exp(-0.046 * (density - dm));
            c4 = 1.0;
        }
        if ((snow->surf_water + snow->pack_water) / snow->depth > 0.01) {
            c4 = 2.0; /* presence of wet snow */
        }
        ddz1 = -param.SNOW_DENS_C2 * c3 * c4 * dexpf;

        /* Compaction due to overburden */
        // swq in this context is the amount of snow whose weight contributes
        // to compaction
        f = param.SNOW_DENS_F;

        /* Currently VIC essentially has only one layer of snow, so compaction
           due to overburden will come mostly from new snowfall. */
        swq = new_snow / MM_PER_M + f * sswq;

        if (new_snow > 0.0) {
            Ps = 0.5 * CONST_G * CONST_RHOFW * swq;
            ddz2 = -Ps / param.SNOW_DENS_ETA0 *
                   exp(-(-param.SNOW_DENS_C5 *
                         (Tavg - CONST_TKFRZ) + param.SNOW_DENS_C6 * density));
        }
        else {
            ddz2 = 0.0;
        }

        /* Calculate compaction rate and new snow density */
        CR = -ddz1 - ddz2;
        density = density * (1. + CR * dt);
    }
    else if (options.SNOW_DENSITY == DENS_BRAS) {
        depth = snow->depth;
        swq = sswq;

        /** Compaction of snow pack by new snow fall **/
        /** Bras pg. 257 **/

        if (new_snow > 0) {
            if (depth > 0.) {
                /* Compact current snowpack by weight of new snowfall */
                delta_depth =
                    (((new_snow / 25.4) * (depth / 0.0254)) / (swq / 0.0254) *
                     pow((depth / 0.0254) / 10., 0.35)) * 0.0254;
                if (delta_depth > param.SNOW_DENS_MAX_CHANGE * depth) {
                    delta_depth = param.SNOW_DENS_MAX_CHANGE * depth;
                }
                depth_new = new_snow / density_new;
                depth = depth - delta_depth + depth_new;
                swq += new_snow / MM_PER_M;
                density = MM_PER_M * swq / depth;
            }
            else {
                /* no snowpack present, so snow density equals that of new snow */
                density = density_new;
                swq += new_snow / MM_PER_M;
                depth = MM_PER_M * swq / density;
            }
        }
        else {
            density = MM_PER_M * swq / snow->depth;
        }

        /** Densification of the snow pack due to aging **/
        /** based on SNTHRM89 R. Jordan 1991 - used in Bart's DHSVM code **/
        if (depth > 0.) {
            overburden = 0.5 * CONST_G * CONST_RHOFW * swq;
            viscosity = param.SNOW_DENS_ETA0 * exp(
                -param.SNOW_DENS_C5 *
                (Tavg - CONST_TKFRZ) + param.SNOW_DENS_C6 * density);
            delta_depth = overburden / viscosity * depth * dt;
            if (delta_depth > param.SNOW_DENS_MAX_CHANGE * depth) {
                delta_depth = param.SNOW_DENS_MAX_CHANGE * depth;
            }
            depth -= delta_depth;
            density = MM_PER_M * swq / depth;
        }
    }

    return (density);
}

/******************************************************************************
* @brief    Estimate the density of new snow
******************************************************************************/
double
new_snow_density(double air_temp)
{
    extern parameters_struct param;
    extern option_struct     options;

    double                   density_new;

    density_new = 0.0;

    if (options.SNOW_DENSITY == DENS_SNTHRM) {
        density_new = 67.9 + 51.3 * exp(air_temp / 2.6);
    }
    else if (options.SNOW_DENSITY == DENS_BRAS) {
        air_temp = air_temp * 9. / 5. + 32.;
        if (air_temp > 0) {
            density_new = param.SNOW_NEW_SNOW_DENSITY + 1000. *
                          (air_temp / 100.) * (air_temp / 100.);
        }
        else {
            density_new = param.SNOW_NEW_SNOW_DENSITY;
        }
    }

    return (density_new);
}

/******************************************************************************
* @brief    This subroutine computes the snow pack surface albedo.'
*
* @note     Computes albedo as a function of snow age and season, based on the
*           algorithm of the US Army Corps of Engineers.
******************************************************************************/
double
snow_albedo(double new_snow,
            double swq,
            double albedo,
            double cold_content,
            double dt,
            int    last_snow,
            bool   MELTING)
{
    extern parameters_struct param;

    /** New Snow **/
    if (new_snow > param.SNOW_TRACESNOW && cold_content < 0.0) {
        albedo = param.SNOW_NEW_SNOW_ALB;
    }
    /** Aged Snow **/
    else if (swq > 0.0) {
        /* Accumulation season */
        if (cold_content < 0.0 && !MELTING) {
            albedo = param.SNOW_NEW_SNOW_ALB * pow(param.SNOW_ALB_ACCUM_A,
                                                   pow((double) last_snow * dt /
                                                       SEC_PER_DAY,
                                                       param.SNOW_ALB_ACCUM_B));
        }
        /* Melt Season */
        else {
            albedo = param.SNOW_NEW_SNOW_ALB * pow(param.SNOW_ALB_THAW_A,
                                                   pow((double) last_snow * dt /
                                                       SEC_PER_DAY,
                                                       param.SNOW_ALB_THAW_B));
        }
    }
    else {
        /* No snow falling or present */
        albedo = 0;
    }

    return(albedo);
}
