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
vic_restore_rout_extension(nameid_struct   *init_state_file,
                           metadata_struct *state_metadata)
{
    extern int         mpi_rank;
    extern rout_struct rout;

    size_t             d2start[2];
    size_t             d2count[2];

    // write state variables

    // routing ring
    if (mpi_rank == VIC_MPI_ROOT) {
        d2start[0] = 0;
        d2start[1] = 0;
        d2count[0] = rout.rout_param.full_time_length;
        d2count[1] = rout.rout_param.n_outlets;

        get_nc_field_double(
            init_state_file,
            state_metadata[N_STATE_VARS + STATE_ROUT_RING].varname,
            d2start, d2count, rout.ring);
    }
}
