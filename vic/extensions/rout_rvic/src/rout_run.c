/******************************************************************************
 * @section DESCRIPTION
 *
 * Run routing over the domain.
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

#include <rout.h>

/******************************************************************************
* @brief        This subroutine controls the RVIC convolution.
******************************************************************************/
void
rout_run(void)
{
    extern int                 mpi_rank;
    extern double           ***out_data;
    extern global_param_struct global_param;
    extern domain_struct       local_domain;
    extern domain_struct       global_domain;
    extern rout_struct         rout;
    double                    *dvar_global = NULL;
    double                    *dvar_local = NULL;
    double                    *dvar_global_total_runoff = NULL;
    double                    *dvar_local_total_runoff = NULL;
    size_t                     i;
    size_t                     i_source;
    size_t                     i_outlet;
    size_t                     i_timestep;
    size_t                     j_timestep;
    int                        offset; /*2d indicies*/
    size_t                     i_ring;
    size_t                     i_uh; /*1d indicies*/

    // Gather runoff and baseflow from local nodes to the masternode
    // and add them together in total_runoff
    dvar_global_total_runoff =
        malloc(global_domain.ncells_total * sizeof(*dvar_global_total_runoff));
    check_alloc_status(dvar_global_total_runoff, "Memory allocation error.");
    dvar_local_total_runoff =
        malloc(local_domain.ncells_active * sizeof(*dvar_local_total_runoff));
    check_alloc_status(dvar_local_total_runoff, "Memory allocation error.");

    // Read from out_data...
    for (i = 0; i < local_domain.ncells_active; i++) {
        dvar_local_total_runoff[i] = out_data[i][OUT_RUNOFF][0] +
                                     out_data[i][OUT_BASEFLOW][0];
    }

    gather_put_var_double(dvar_global_total_runoff, dvar_local_total_runoff);

    free(dvar_local_total_runoff);

    // Do the routing on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        debug("RVIC");

        dvar_global = malloc(global_domain.ncells_total * sizeof(*dvar_global));
        check_alloc_status(dvar_global, "Memory allocation error.");

        // Zero out current ring
        // in python: (from variables.py) self.ring[tracer][0, :] = 0.
        for (i_outlet = 0; i_outlet < rout.rout_param.n_outlets; i_outlet++) {
            rout.ring[i_outlet] = 0.0;
        }

        // Equivalent to Fortran 90 cshift function, in python: (from variables.py)
        // self.ring[tracer] = np.roll(self.ring[tracer], -1, axis=0)
        cshift(rout.ring, rout.rout_param.full_time_length,
               rout.rout_param.n_outlets, 0,
               1);

        /*Loop through all sources*/
        for (i_source = 0; i_source < rout.rout_param.n_sources; i_source++) {
            i_outlet = rout.rout_param.source2outlet_ind[i_source];
            offset = rout.rout_param.source_time_offset[i_source];

            /* Do the convolution */
            // iTimestep is the position in the unit hydrograph
            // jTimestep is the position in the ring
            for (i_timestep = 0; i_timestep < rout.rout_param.n_timesteps;
                 i_timestep++) {
                j_timestep = i_timestep + offset;

                // index locations
                i_ring = (j_timestep * rout.rout_param.n_outlets) + i_outlet;
                i_uh = (i_timestep * rout.rout_param.n_sources) + i_source;

                rout.ring[i_ring] += rout.rout_param.unit_hydrograph[i_uh] *
                                     dvar_global_total_runoff[rout.rout_param.
                                                              source_VIC_index[
                                                                  i_source]];
            }
        }

        // Write to dvar_global prior to scattering over local domains...
        for (i_outlet = 0; i_outlet < rout.rout_param.n_outlets; i_outlet++) {
            dvar_global[rout.rout_param.outlet_VIC_index[i_outlet]] =
                rout.ring[i_outlet] *
                global_domain.locations[rout.rout_param.outlet_VIC_index[
                                            i_outlet]].area /
                (MM_PER_M * global_param.dt);
        }

        // cleanup
        free(dvar_global);
    }

    // Scatter dvar_global to local domains (dvar_local)...
    dvar_local = malloc(local_domain.ncells_active * sizeof(*dvar_local));
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
