/******************************************************************************
 * @section DESCRIPTION
 *
 * calculates nitrogen scaling factors for all canopy layers, following eqns
 * 106 and 107 in Knorr 1997.
 *
 * Note: this should only be applied to veg types that have a canopy, e.g.
 * trees and shrubs, but not grass or tundra vegetation.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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

#include <vic_def.h>
#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate nitrogen scaling factors for all canopy layers.
 *****************************************************************************/
void
calc_Nscale_factors(char       NscaleFlag,
                    double    *CanopLayerBnd,
                    double     LAItotal,
                    double     lat,
                    double     lng,
                    double     time_zone_lng,
                    dmy_struct dmy,
                    double    *NscaleFactor)
{
    extern option_struct     options;
    extern parameters_struct param;

    dmy_struct               dmy_tmp;
    double                   coszen_noon;
    double                   k12;
    size_t                   cidx; // canopy layer index

    /* Compute solar zenith angle at local noon */
    dmy_tmp.year = dmy.year;
    dmy_tmp.month = dmy.month;
    dmy_tmp.day = dmy.day;
    dmy_tmp.day_in_year = dmy.day_in_year;
    dmy_tmp.hour = 12;
    coszen_noon = compute_coszen(lat, lng, time_zone_lng, dmy_tmp);
    if (coszen_noon < param.PHOTO_ZENITHMINPAR) {
        coszen_noon = param.PHOTO_ZENITHMINPAR;
    }

    /* Extinction factor; eqn 119c in Knorr 1997 */
    k12 = 0.5 / coszen_noon;

    /* Condition: LAI > LaiLimit; eqns 107 and 108 in Knorr 1997 */
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
        if (NscaleFlag && LAItotal > param.PHOTO_LAILIMIT && cidx > 0) {
            NscaleFactor[cidx] = exp(-k12 * CanopLayerBnd[cidx - 1] * LAItotal);
            if (NscaleFactor[cidx] < 1e-10) {
                NscaleFactor[cidx] = 1e-10;
            }
        }
        else {
            NscaleFactor[cidx] = 1;
        }
    }
}
