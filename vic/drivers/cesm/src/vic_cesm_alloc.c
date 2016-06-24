/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
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
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void
vic_cesm_alloc(void)
{
    extern x2l_data_struct *x2l_vic;
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;

    debug("In vic_cesm_alloc");

    // allocate memory for x2l_vic structure
    x2l_vic = malloc(local_domain.ncells_active * sizeof(*x2l_vic));
    check_alloc_status(x2l_vic, "Memory allocation error.");
    // initialize x2l data
    initialize_x2l_data();

    // allocate memory for l2x_vic structure
    l2x_vic = malloc(local_domain.ncells_active * sizeof(*l2x_vic));
    check_alloc_status(l2x_vic, "Memory allocation error.");
    // initialize l2x data
    initialize_l2x_data();

    // allocate the rest of the image mode structures
    vic_alloc();
}
