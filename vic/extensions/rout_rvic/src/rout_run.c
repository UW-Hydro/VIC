/******************************************************************************
 * @section DESCRIPTION
 *
 * Run routing over the domain.
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
#include <rout.h>
#include <assert.h>
#include <vic_driver_image.h>

void
rout_run(void)
{
    log_info("In Routing Lohmann Model");

    extern rout_struct       rout;
//    extern out_data_struct **out_data;
    extern double           ***out_data;

    size_t                   iSource, iOutlet, iTimestep, jTimestep;
    int                      offset; /*2d indicies*/
    size_t                   iRing, iUH; /*1d indicies*/

    // Zero out current ring
    // in python: (from variables.py) self.ring[tracer][0, :] = 0.
    for (iOutlet = 0; iOutlet < rout.rout_param.nOutlets; iOutlet++) {
        rout.ring[iOutlet] = 0.0;
    }

    // Equivalent to Fortran 90 cshift function, in python: (from variables.py) self.ring[tracer] = np.roll(self.ring[tracer], -1, axis=0)
    cshift(rout.ring, rout.rout_param.full_time_length,
           rout.rout_param.nOutlets, 0,
           1);

    /*Loop through all sources*/
    for (iSource = 0; iSource < rout.rout_param.nSources; iSource++) {
        iOutlet = rout.rout_param.source2outlet_ind[iSource];
        offset = rout.rout_param.source_time_offset[iSource];

        /* Do the convolution */
        // iTimestep is the position in the unit hydrograph
        // jTimestep is the position in the ring
        for (iTimestep = 0; iTimestep < rout.rout_param.nTimesteps;
             iTimestep++) {
            jTimestep = iTimestep + offset;

            // index locations
            iRing = (jTimestep * rout.rout_param.nOutlets) + iOutlet;
            iUH = (iTimestep * rout.rout_param.nSources) + iSource;

            rout.ring[iRing] += (rout.rout_param.unit_hydrograph[iUH] *
                                 (out_data[rout.rout_param.source_VIC_index[
                                               iSource]][
                                      OUT_RUNOFF][0] +
                                  out_data[rout.rout_param.source_VIC_index[
                                               iSource]][
                                      OUT_BASEFLOW][0]));
        }
    }

    // Write to output struct...
    for (iOutlet = 0; iOutlet < rout.rout_param.nOutlets; iOutlet++) {
        out_data[rout.rout_param.outlet_VIC_index[iOutlet]][OUT_DISCHARGE][0] = rout.ring[iOutlet];
    }
}
