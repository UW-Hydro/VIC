/******************************************************************************
 * @section DESCRIPTION
 *
 * Write output data.
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
 * @brief    Write output data and convert units if necessary.
 *****************************************************************************/
void
write_output(stream_struct **streams,
             dmy_struct     *dmy)
{
    extern option_struct options;

    size_t               stream_idx;

    // Write data
    for (stream_idx = 0; stream_idx < options.Noutstreams; stream_idx++) {
        if (raise_alarm(&(*streams)[stream_idx].agg_alarm, dmy)) {
            write_data(&((*streams)[stream_idx]));
            reset_stream((&(*streams)[stream_idx]), dmy);
        }
    }
}
