/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine creates an array of structures that contain information about a
 * cell's states and fluxes.
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
 * @brief    Creates an array of structures that contain information about a
 *           cell's states and fluxes.
 *****************************************************************************/
all_vars_struct
make_all_vars(size_t nveg)
{
    all_vars_struct temp;
    size_t          Nitems;

    Nitems = nveg + 1;

    temp.snow = make_snow_data(Nitems);
    temp.energy = make_energy_bal(Nitems);
    temp.veg_var = make_veg_var(Nitems);
    temp.cell = make_cell_data(Nitems);

    return (temp);
}
