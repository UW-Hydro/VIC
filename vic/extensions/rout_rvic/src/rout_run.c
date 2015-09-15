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

void rout_run(void) {
    //log_info("In Routing Lohmann Model");

    extern rout_struct rout;
    extern out_data_struct **out_data;
    size_t i, j, s;

    ///TESTS!!!
    printf("\nmapping test!!!...");

    for (j = 0; j < rout.rout_param.iSources; j++) {
        printf("\nsource: %2zu, outlet_ind: %2zu, Lon_source: , %8.4f, Lat_source: %8.4f", j, rout.rout_param.source2outlet_ind[j],
                rout.rout_param.source_lon[j],
                rout.rout_param.source_lat[j]);
    }


//    // even worse stuff here 
//    for (i = 0; i < global_domain.ncells; i++) {
//        rout.rout_param.aggrunin[i] = out_data[i][OUT_RUNOFF].data[0] + out_data[i][OUT_BASEFLOW].data[0];
//    }

    printf("\nThe hydrograph, just printed...\n");
    print_array(rout.rout_param.unit_hydrograph, rout.rout_param.iSubsetLength, rout.rout_param.iSources);

    // Filling the ring with dummy data
    for (j = 0; j < rout.rout_param.iOutlets; j++) {
        for (i = 0; i < rout.rout_param.iSubsetLength; i++) {
            //     rout.ring[i * rout.rout_param.iOutlets + j]=i * rout.rout_param.iOutlets + j;
        }
    }

    // Zero out current ring
    // in python: (from variables.py) self.ring[tracer][0, :] = 0.
    printf("\nThe Ring (ZERO)...\n");
    for (i = 0; i < rout.rout_param.iOutlets; i++) {
        rout.ring[i] = 0.0;
    }

    print_array(rout.ring, rout.rout_param.iSubsetLength, rout.rout_param.iOutlets);

    // Equivalent to Fortran 90 cshift function, in python: (from variables.py) self.ring[tracer] = np.roll(self.ring[tracer], -1, axis=0)
    cshift(rout.ring, rout.rout_param.iSubsetLength, rout.rout_param.iOutlets, 0, 1);

    printf("\nThe Ring, c-shifted...\n");
    print_array(rout.ring, rout.rout_param.iSubsetLength, rout.rout_param.iOutlets);

    printf("\nDoing convolution test!\n");
    int offset, outlet; /*2d indicies*/
    int rind, uhind; /*1d indicies*/

    /*Loop through all sources*/
    for (s = 0; s < rout.rout_param.iSources; s++) {

        outlet = rout.rout_param.source2outlet_ind[s];
        offset = rout.rout_param.source_time_offset[s];

        /* Do the convolution */
        // i is the position in the unit hydrograph
        // j is the position in the ring
        for (i = 0; i < rout.rout_param.iSubsetLength; i++) {
            j = i + offset;

            //1d index locations
            rind = j * rout.rout_param.iOutlets + outlet;
            uhind = i * rout.rout_param.iSources + s;

            //            ring[rind] += unit_hydrograph[uhind] * aggrunin[xyind];
            //            ring[rind] += unit_hydrograph[uhind];
            //ring[rind] += out_data[rout.rout_param.source_VIC_index[s]][OUT_RUNOFF].data[0] + out_data[rout.rout_param.source_VIC_index[s]][OUT_BASEFLOW].data[0];
            rout.ring[rind] += (rout.rout_param.unit_hydrograph[uhind] *
                    (out_data[rout.rout_param.source_VIC_index[s]][OUT_RUNOFF].data[0] +
                    out_data[rout.rout_param.source_VIC_index[s]][OUT_BASEFLOW].data[0]));
        }
    }

    printf("\nThe Ring, after convolution...\n");
    print_array(rout.ring, rout.rout_param.iSubsetLength, rout.rout_param.iOutlets);

    // Write to output struct...
    for (s = 0; s < rout.rout_param.iOutlets; s++) {
        out_data[rout.rout_param.outlet_VIC_index[s]][OUT_DISCHARGE].aggdata[0] = rout.ring[s];
    }
}
