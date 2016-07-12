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
rout_run(void) {
    extern int mpi_rank;
    extern double ***out_data;
    extern global_param_struct global_param;
    extern domain_struct local_domain;
    extern domain_struct global_domain;
    double *dvar_global = NULL;
    double *dvar_local = NULL;
    double *dvar_global_total_runoff = NULL;
    double *dvar_local_total_runoff = NULL;
    size_t i;

    //////////////////
    // Gather runoff and baseflow from local nodes to the masternode
    // and add them together in total_runoff
    dvar_global_total_runoff = malloc(global_domain.ncells_total * sizeof (*dvar_global_total_runoff));
    check_alloc_status(dvar_global_total_runoff, "Memory allocation error.");
    dvar_local_total_runoff = malloc(local_domain.ncells_active * sizeof (*dvar_local_total_runoff));
    check_alloc_status(dvar_local_total_runoff, "Memory allocation error.");
    
    // Read from out_data...
    for (i = 0; i < local_domain.ncells_active; i++) {
        dvar_local_total_runoff[i] = out_data[i][OUT_RUNOFF][0] + out_data[i][OUT_BASEFLOW][0];
    }

    gather_put_var_double(dvar_global_total_runoff, dvar_local_total_runoff);

    free(dvar_local_total_runoff);
    //////////////////////

    // Do the routing on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        debug("In Routing Lohmann Model");

        extern rout_struct rout;
        extern domain_struct global_domain;

        dvar_global = malloc(global_domain.ncells_total * sizeof (*dvar_global));
        check_alloc_status(dvar_global, "Memory allocation error.");

        size_t iSource, iOutlet, iTimestep, jTimestep;
        int offset; /*2d indicies*/
        size_t iRing, iUH; /*1d indicies*/

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

                rout.ring[iRing] += rout.rout_param.unit_hydrograph[iUH] *
                        dvar_global_total_runoff[rout.rout_param.source_VIC_index[iSource]];
            }
        }

        // Write to dvar_global prior to scattering over local domains...
        for (iOutlet = 0; iOutlet < rout.rout_param.nOutlets; iOutlet++) {
            dvar_global[rout.rout_param.outlet_VIC_index[iOutlet]] =
                    rout.ring[iOutlet] *
                    global_domain.locations[rout.rout_param.outlet_VIC_index[iOutlet]].area
                    / (MM_PER_M * global_param.dt);
        }

        // cleanup
        free(dvar_global);
    }

    // Scatter dvar_global to local domains (dvar_local)...
    dvar_local = malloc(local_domain.ncells_active * sizeof (*dvar_local));
    check_alloc_status(dvar_local, "Memory allocation error.");

    get_scatter_var_double(dvar_global, dvar_local);

    // Write to output struct...
    for (i = 0; i < local_domain.ncells_active; i++) {
        out_data[i][OUT_DISCHARGE][0] = dvar_local[i];
    }

    // cleanup
    free(dvar_local);
    free(dvar_global_total_runoff);

}
