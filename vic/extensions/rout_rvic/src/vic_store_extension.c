/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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

#include <vic_driver_shared_image.h>
#include <rout.h>

/******************************************************************************
 * @brief    Save model state.
 *****************************************************************************/
void
vic_store_rout_extension(nc_file_struct *nc_state_file)
{
    extern int         mpi_rank;
    extern rout_struct rout;

    int                status;
    size_t             d2start[2];
    nc_var_struct     *nc_var;

    // write state variables

    // routing ring
    if (mpi_rank == VIC_MPI_ROOT) {
        d2start[0] = 0;
        d2start[1] = 0;
        nc_var = &(nc_state_file->nc_vars[N_STATE_VARS + STATE_ROUT_RING]);

        status =
            nc_put_vara_double(nc_state_file->nc_id, nc_var->nc_varid, d2start,
                               nc_var->nc_counts,
                               rout.ring);
        check_nc_status(status, "Error writing values.");
    }
}

/******************************************************************************
 * @brief   Setup state file netcdf structure
 *****************************************************************************/
void
set_nc_state_file_info_rout_extension(nc_file_struct *nc_state_file)
{
    extern rout_struct rout;

    // set ids to MISSING
    nc_state_file->outlet_dimid = MISSING;
    nc_state_file->routing_timestep_dimid = MISSING;

    // set dimension sizes
    nc_state_file->outlet_size = rout.rout_param.n_outlets;
    nc_state_file->routing_timestep_size = rout.rout_param.full_time_length;
}

/******************************************************************************
 * @brief   Setup state variable dimensions, types, etc.
 *****************************************************************************/
void
set_nc_state_var_info_rout_extension(nc_file_struct *nc)
{
    size_t i;
    size_t j;

    for (i = N_STATE_VARS; i < (N_STATE_VARS + N_STATE_VARS_EXT); i++) {
        nc->nc_vars[i].nc_varid = i;
        for (j = 0; j < MAXDIMS; j++) {
            nc->nc_vars[i].nc_dimids[j] = -1;
            nc->nc_vars[i].nc_counts[j] = 0;
        }
        nc->nc_vars[i].nc_dims = 0;
        nc->nc_vars[i].nc_type = NC_DOUBLE;

        // Set the number of dimensions and dimids for each state variable
        switch (i) {
        case (N_STATE_VARS + STATE_ROUT_RING):
            // 2d vars [routing_timestep, outlet]
            nc->nc_vars[i].nc_dims = 2;
            nc->nc_vars[i].nc_dimids[0] = nc->routing_timestep_dimid;
            nc->nc_vars[i].nc_dimids[1] = nc->outlet_dimid;
            nc->nc_vars[i].nc_counts[0] = nc->routing_timestep_size;
            nc->nc_vars[i].nc_counts[1] = nc->outlet_size;
            break;
        default:
            log_err("state variable %zu not found when setting dimensions", i);
        }

        if (nc->nc_vars[i].nc_dims > MAXDIMS) {
            log_err("Too many dimensions specified in variable %zu", i);
        }
    }
}

/******************************************************************************
 * @brief   Initialize state file by creating dimensions, variables,
            and adding metadata.
 *****************************************************************************/
void
initialize_state_file_rout_extension(char           *filename,
                                     nc_file_struct *nc_state_file)
{
    int status;

    // Add routing dimensions
    status = nc_def_dim(nc_state_file->nc_id, "outlet",
                        nc_state_file->outlet_size,
                        &(nc_state_file->outlet_dimid));
    check_nc_status(status, "Error defining outlet in %s", filename);

    status = nc_def_dim(nc_state_file->nc_id, "routing_timestep",
                        nc_state_file->routing_timestep_size,
                        &(nc_state_file->routing_timestep_dimid));
    check_nc_status(status, "Error defining routing_timestep in %s", filename);
}
