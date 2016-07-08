/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for Routing structures.
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
#include <vic_run.h>
#include <vic_driver_image.h>
#include <rout.h>

/******************************************************************************
 * @brief   Gather double precision variable
 * @details Values are gathered to the master node
 *****************************************************************************/
void gather_put_var_double(double *dvar, double *var_local) {
    extern MPI_Comm MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int mpi_rank;
    extern int *mpi_map_global_array_offsets;
    extern int *mpi_map_local_array_sizes;
    extern size_t *filter_active_cells;
    extern size_t *mpi_map_mapping_array;
    int status;
    double *dvar_gathered = NULL;
    double *dvar_remapped = NULL;
    size_t               grid_size;
    size_t               i;
    
    if (mpi_rank == VIC_MPI_ROOT) {
        grid_size = global_domain.n_nx * global_domain.n_ny;

        for (i = 0; i < grid_size; i++) {
            dvar[i] = 0;
        }

        dvar_gathered =
                malloc(global_domain.ncells_active * sizeof (*dvar_gathered));
        check_alloc_status(dvar_gathered, "Memory allocation error.");

        dvar_remapped =
                malloc(global_domain.ncells_active * sizeof (*dvar_remapped));
        check_alloc_status(dvar_remapped, "Memory allocation error.");
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var_local, local_domain.ncells_active, MPI_DOUBLE,
            dvar_gathered, mpi_map_local_array_sizes,
            mpi_map_global_array_offsets, MPI_DOUBLE,
            VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");
    if (mpi_rank == VIC_MPI_ROOT) {
        // remap the array
        map(sizeof (double), global_domain.ncells_active, NULL,
                mpi_map_mapping_array, dvar_gathered, dvar_remapped);
        // expand to full grid size
        map(sizeof (double), global_domain.ncells_active, NULL,
                filter_active_cells, dvar_remapped, dvar);

        // cleanup
        free(dvar_gathered);
        free(dvar_remapped);
    }
}

