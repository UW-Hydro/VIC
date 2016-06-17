/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine opens files for soil, vegetation, and global parameters.
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
 * @brief    This routine opens files for soil, vegetation, and global
 *           parameters.
 *****************************************************************************/
void
check_files(filep_struct     *filep,
            filenames_struct *fnames)
{
    extern option_struct options;
    extern FILE          *open_file(char string[], char type[]);

    filep->soilparam = open_file(fnames->soil, "r");
    filep->veglib = open_file(fnames->veglib, "r");
    filep->vegparam = open_file(fnames->veg, "r");
    if (options.SNOW_BAND > 1) {
        filep->snowband = open_file(fnames->snowband, "r");
    }
    if (options.LAKES) {
        filep->lakeparam = open_file(fnames->lakeparam, "r");
    }
}
