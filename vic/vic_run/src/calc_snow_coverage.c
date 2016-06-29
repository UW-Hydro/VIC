/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the current fraction of the vegetation band that is
 * covered with snow.  The snow distribution is assumed to be uniform with a
 * slope based on the value of max_snow_distrib_slope.  The
 * original value was based on field observations from the University of Minnesota's
 * Rosemount Agricultural Experiment station (see dissertation by Keith
 * Cherkauer, 2001).

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
 * @brief    Compute the current fraction of the vegetation band that is
 *           covered with snow.
 *****************************************************************************/
double
calc_snow_coverage(bool   *store_snow,
                   double  max_snow_distrib_slope,
                   double  old_coverage,
                   double  swq,
                   double  old_swq,
                   double  depth,
                   double  old_depth,
                   double  melt,
                   double *max_snow_depth,
                   double  snowfall,
                   double *store_swq,
                   double *snow_distrib_slope,
                   double *store_coverage)
{
    double old_max_snow_depth;
    double coverage;

    /*************************************
       New snow falls on partial snowpack
    *************************************/

    if (snowfall > 0) {
        /*************************
           Continued accumulation
        *************************/

        coverage = 1;

        if (*store_snow) {
            // store coverage fraction before it is buried
            if (*store_swq == 0 && old_coverage < 1) {
                *store_coverage = old_coverage;
            }
            else if (*store_swq == 0) {
                *store_coverage = 1;
            }

            /* store snow falling over partial snowpack */
            *store_swq += swq - old_swq;

            if (depth >= max_snow_distrib_slope / 2.) {
                /* snow accumulation deep enough to remove memory of previous
                   melt distribution */
                *store_snow = false;
                *store_swq = 0;
                *snow_distrib_slope = 0;
                *store_coverage = 1;
            }
        }
        else if (old_coverage < 1) {
            /* store new snow fall over old distribution */
            *store_snow = true;
            *store_swq = swq - old_swq;
        }
    }
    else if (melt > 0) {
        /***********************************************
           Snowpack begins to melt or continues melting
        ***********************************************/

        if (*store_swq > 0 && swq < old_swq) {
            /* Melt thin snowfall off previous distribution */
            *store_swq += swq - old_swq;

            if (*store_swq <= 0) {
                /* Snowpack cover has melted - clear storage */
                *store_swq = 0;
                // restore buried cover fraction
                old_coverage = *store_coverage;
                *store_coverage = 1;
            }
        }

        if (*store_swq == 0) {
            /* compute the current coverage fraction */

            if ((*snow_distrib_slope) == 0) {
                /* compute new distribution function */
                if (old_depth > max_snow_distrib_slope / 2.) {
                    *snow_distrib_slope = -max_snow_distrib_slope;
                }
                else {
                    *snow_distrib_slope = -2. * old_depth;
                }
                *max_snow_depth = -(*snow_distrib_slope);
                *store_snow = true;
            }

            /* if currently raining, new swq may be higher than previous
               maximum swq even if melt occurs.  reset maximum swq if this
               occurs. */

            old_max_snow_depth = *max_snow_depth;
            *max_snow_depth = 2. * depth;

            if (*max_snow_depth < old_max_snow_depth || old_max_snow_depth ==
                0) {
                /* melt has occured, reduce coverage fraction */
                coverage = -(*max_snow_depth) / (*snow_distrib_slope);
                coverage = (coverage > 1) ? 1 : coverage;
            }
            else {
                /* rain or sublimation has increased swq,
                   coverage fraction does not change */
                coverage = old_coverage;
            }

            /** NOTE if (swq - old_swq) > 0, but there was no snowfall,
               there was rain so the coverage fraction should not change */
        }
        else {
            coverage = old_coverage;
        }
    }
    else {
        /* no change to snowpack */
        coverage = old_coverage;
    }

    return(coverage);
}
