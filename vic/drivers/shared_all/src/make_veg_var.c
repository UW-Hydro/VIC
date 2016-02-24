/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of vegitation variables for each vegitation type.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of vegitation variables for each vegitation type.
 *****************************************************************************/
veg_var_struct **
make_veg_var(size_t veg_type_num)
{
    extern option_struct options;

    size_t               i, j;
    veg_var_struct     **temp = NULL;

    temp = calloc(veg_type_num, sizeof(*temp));
    if (temp == NULL) {
        log_err("Memory allocation error in make_veg_var().");
    }

    for (i = 0; i < veg_type_num; i++) {
        temp[i] = calloc(options.SNOW_BAND, sizeof(*(temp[i])));
        if (temp[i] == NULL) {
            log_err("Memory allocation error in make_veg_var().");
        }

        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                temp[i][j].NscaleFactor = calloc(options.Ncanopy,
                                                 sizeof(*(temp[i][j].
                                                          NscaleFactor)));
                if (temp[i][j].NscaleFactor == NULL) {
                    log_err("Memory allocation error in make_veg_var().");
                }
                temp[i][j].aPARLayer = calloc(options.Ncanopy,
                                              sizeof(*(temp[i][j].aPARLayer)));
                if (temp[i][j].aPARLayer == NULL) {
                    log_err("Memory allocation error in make_veg_var().");
                }
                temp[i][j].CiLayer = calloc(options.Ncanopy,
                                            sizeof(*(temp[i][j].CiLayer)));
                if (temp[i][j].CiLayer == NULL) {
                    log_err("Memory allocation error in make_veg_var().");
                }
                temp[i][j].rsLayer = calloc(options.Ncanopy,
                                            sizeof(*(temp[i][j].rsLayer)));
                if (temp[i][j].rsLayer == NULL) {
                    log_err("Memory allocation error in make_veg_var().");
                }
            }
        }
    }

    return temp;
}
