/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine handles the startup tasks for the image driver.
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

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Wrapper function for VIC startup tasks.
 *****************************************************************************/
void
vic_image_start(void)
{
    extern filep_struct     filep;
    extern filenames_struct filenames;
    extern int              mpi_rank;

    // Initialize structures
    initialize_global_structures();

    if (mpi_rank == VIC_MPI_ROOT) {
        // Read the global parameter file
        filep.globalparam = open_file(filenames.global, "r");
        get_global_param(filep.globalparam);
    }

    // initialize image mode structures and settings
    vic_start();
}
