/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine resets the values of all output variables to 0.
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
 * @brief    This routine resets the values of all output variables to 0.
 *****************************************************************************/
void
zero_output_list(double **out_data)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 varid, i;

    for (varid = 0; varid < N_OUTVAR_TYPES; varid++) {
        for (i = 0; i < out_metadata[varid].nelem; i++) {
            out_data[varid][i] = 0.;
        }
    }
}
