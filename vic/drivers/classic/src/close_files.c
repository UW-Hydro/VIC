/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine closes all forcing data files, and output files.
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
 * @brief    This routine closes all forcing data files, and output files.
 *****************************************************************************/
void
close_files(filep_struct   *filep,
            stream_struct **streams)
{
    extern option_struct options;

    size_t               streamnum;

    /**********************
       Close All Input Files
    **********************/

    fclose(filep->forcing[0]);

    if (filep->forcing[1] != NULL) {
        fclose(filep->forcing[1]);
    }

    /*******************
       Close Output Files
    *******************/
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        fclose((*streams)[streamnum].fh);
        if ((*streams)[streamnum].compress) {
            compress_files((*streams)[streamnum].filename,
                           (*streams)[streamnum].compress);
        }
    }
}
