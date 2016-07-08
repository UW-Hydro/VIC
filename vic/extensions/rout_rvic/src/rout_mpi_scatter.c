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
 * @brief   Scatter double precision variable
 * @details values from master node are scattered to the local nodes
 *****************************************************************************/
void get_scatter_var_double(double *dvar, double *var_local) {

    extern MPI_Comm MPI_COMM_VIC;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int mpi_rank;
    extern int *mpi_map_global_array_offsets;
    extern int *mpi_map_local_array_sizes;
    extern size_t *filter_active_cells;
    extern size_t *mpi_map_mapping_array;
    int status;
    double *dvar_filtered = NULL;
    double *dvar_mapped = NULL;

    if (mpi_rank == VIC_MPI_ROOT) {
        dvar_filtered =
                malloc(global_domain.ncells_active * sizeof (*dvar_filtered));
        check_alloc_status(dvar_filtered, "Memory allocation error.");

        dvar_mapped =
                malloc(global_domain.ncells_active * sizeof (*dvar_mapped));
        check_alloc_status(dvar_mapped, "Memory allocation error.");

        // filter the active cells only
        map(sizeof (double), global_domain.ncells_active, filter_active_cells,
                NULL, dvar, dvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof (double), global_domain.ncells_active, mpi_map_mapping_array,
                NULL, dvar_filtered, dvar_mapped);
        free(dvar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(dvar_mapped, mpi_map_local_array_sizes,
            mpi_map_global_array_offsets, MPI_DOUBLE,
            var_local, local_domain.ncells_active, MPI_DOUBLE,
            VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        free(dvar_mapped);
    }
}

