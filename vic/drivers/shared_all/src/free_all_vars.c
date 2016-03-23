/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine frees all memory allocated down the all_vars data structure.
 *
 * This include all grid cell specific variables (soil, vegetation, energy,
 * snow).
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
 * @brief    Free all variables.
 *****************************************************************************/
void
free_all_vars(all_vars_struct *all_vars,
              int              Nveg)
{
    extern option_struct options;

    int                  i, j, Nitems;
    size_t               k;

    Nitems = Nveg + 1;

    for (j = 0; j < Nitems; j++) {
        free((char *) all_vars[0].cell[j]);
    }
    free((char *) all_vars[0].cell);
    for (j = 0; j < Nitems; j++) {
        if (options.CARBON) {
            for (k = 0; k < options.SNOW_BAND; k++) {
                free((char *) all_vars[0].veg_var[j][k].NscaleFactor);
                free((char *) all_vars[0].veg_var[j][k].aPARLayer);
                free((char *) all_vars[0].veg_var[j][k].CiLayer);
                free((char *) all_vars[0].veg_var[j][k].rsLayer);
            }
        }
        free((char *)(*all_vars).veg_var[j]);
    }
    free((char *)(*all_vars).veg_var);
    for (j = 0; j < Nitems; j++) {
        free((char *) all_vars[0].energy[j]);
    }
    free((char *) all_vars[0].energy);
    for (i = 0; i < Nitems; i++) {
        free((char *) all_vars[0].snow[i]);
    }
    free((char *) all_vars[0].snow);
}
