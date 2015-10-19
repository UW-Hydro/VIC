/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine makes an array of type cell, which contains soil column
 * variables for a single grid cell.
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
#include <vic_driver_shared.h>

/******************************************************************************
 * @brief    Make an array of type cell, which contains soil column variables
 *           for a single grid cell.
 *****************************************************************************/
cell_data_struct **
make_cell_data(size_t veg_type_num)
{
    extern option_struct options;

    size_t               i;
    cell_data_struct   **temp;

    temp = (cell_data_struct**) calloc(veg_type_num,
                                       sizeof(cell_data_struct*));
    for (i = 0; i < veg_type_num; i++) {
        temp[i] = (cell_data_struct*) calloc(options.SNOW_BAND,
                                             sizeof(cell_data_struct));
    }
    return temp;
}
