/******************************************************************************
 * @section DESCRIPTION
 *
 * Finalize VIC run by freeing memory and closing open files.
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

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/
void
vic_cesm_finalize(void)
{
    extern x2l_data_struct *x2l_vic;
    extern l2x_data_struct *l2x_vic;

    // free VIC/CESM data structures
    free(x2l_vic);
    free(l2x_vic);

    vic_finalize();
}
