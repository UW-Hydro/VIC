/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes the forcing file parameters
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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize Forcing File Parameters
 *****************************************************************************/
void
initialize_forcing_files()
{
    extern param_set_struct param_set;

    int                     i, j;

    /** Initialize forcing file input controls **/

    for (j = 0; j < N_FORCING_TYPES; j++) {
        param_set.TYPE[j].SUPPLIED = false;
        param_set.TYPE[j].SIGNED = 1;
        param_set.TYPE[j].multiplier = 1;
    }
    for (i = 0; i < 2; i++) {
        param_set.FORCE_DT[i] = MISSING;
        param_set.force_steps_per_day[i] = 0;
        param_set.N_TYPES[i] = 0;
        param_set.FORCE_FORMAT[i] = MISSING;
        for (j = 0; j < N_FORCING_TYPES; j++) {
            param_set.FORCE_INDEX[i][j] = MISSING;
        }
    }
}
