/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine frees the list of output files.
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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine frees the memory in the out_data_files array.
 *****************************************************************************/
void
free_out_data_files(out_data_file_struct **out_data_files)
{
    extern option_struct options;
    size_t               filenum;

    for (filenum = 0; filenum < options.Noutfiles; filenum++) {
        free((char*) (*out_data_files)[filenum].varid);
    }
    free((char*) (*out_data_files));
}

