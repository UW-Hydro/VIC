/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize output structures.
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
 * @brief    Perform temporal aggregation on stream data
 *****************************************************************************/
void
agg_stream_data(stream_struct *stream,
                double       **out_data) {

    extern out_metadata_struct *out_metadata;

    size_t       i;
    size_t       j;
    size_t       nelem;
    unsigned int varid;

    stream->counter++;

    for (i = 0; i < stream->nvars; i++) {

        varid = stream->varid[i];
        nelem = out_metadata[varid].nelem;

        // Instantaneous at the beginning of the period
        if ((stream->aggtype[i] == AGG_TYPE_END) &&
            (stream->counter == stream->nextagg)) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] = out_data[i][j];
            }
        }
        // Instantaneous at the end of the period
        else if ((stream->aggtype[i] == AGG_TYPE_BEG) &&
                 (stream->counter == 1)) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] = out_data[i][j];
            }
        }
        // Sum over the period
        else if ((stream->aggtype[i] == AGG_TYPE_SUM) ||
                 (stream->aggtype[i] == AGG_TYPE_AVG)) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] += out_data[i][j];
            }
        }
        // Maximum over the period
        else if (stream->aggtype[i] == AGG_TYPE_MAX) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] += min(stream->aggdata[i][j][0], out_data[i][j]);
            }
        }
        // Minimum over the period
        else if (stream->aggtype[i] == AGG_TYPE_MAX) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] += max(stream->aggdata[i][j][0], out_data[i][j]);
            }
        }
        // Average over the period if counter is full
        if ((stream->aggtype[i] == AGG_TYPE_AVG) &&
            (stream->counter == stream->nextagg)) {
            for (j = 0; j < nelem; j++) {
                stream->aggdata[i][j][0] /= (double) stream->nextagg;
            }
        }
    }
}
