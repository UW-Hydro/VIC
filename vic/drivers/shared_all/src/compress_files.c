/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine compresses the file "string" using a system call.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This subroutine compresses the file "string" using a system call.
 *****************************************************************************/
void
compress_files(char      string[],
               short int level)
{
    char command[MAXSTRING];

    // Compress the file
    if (level == COMPRESSION_LVL_DEFAULT) {
        sprintf(command, "nice gzip -f %s &", string);
    }
    else if (level != COMPRESSION_LVL_UNSET) {
        sprintf(command, "nice gzip -%d -f %s &", level, string);
    }
    else if (level <= 0) {
        log_err("Invalid compression level for gzip, must be an integer 1-9");
    }

    system(command);
}
