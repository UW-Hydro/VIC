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
rout_run(void) {
    extern int mpi_rank;
    extern int *mpi_map_global_array_offsets;
    extern int *mpi_map_local_array_sizes;
    extern size_t *filter_active_cells;
    extern size_t *mpi_map_mapping_array;
    extern MPI_Comm MPI_COMM_VIC;

    extern double ***out_data;
    extern domain_struct local_domain;
    extern domain_struct global_domain;
    int status;
    double *var_local_runoff = NULL;
    double *var2 = NULL;
    double *dvar = NULL;
    double *dvar_runoff = NULL;
    double *dvar_filtered = NULL;
    double *dvar_mapped = NULL;
    double *dvar_gathered = NULL;
    double *dvar_remapped = NULL;
    size_t i;

    // allocate memory for variables to be read
    var_local_runoff = malloc(local_domain.ncells_active * sizeof (*var_local_runoff));
    check_alloc_status(var_local_runoff, "Memory allocation error.");

    // Read from out_data...
    for (i = 0; i < local_domain.ncells_active; i++) {
        var_local_runoff[i] = out_data[i][OUT_RUNOFF][0] +
                out_data[i][OUT_BASEFLOW][0];
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        dvar_runoff = malloc(global_domain.ncells_total * sizeof (*dvar_runoff));
        check_alloc_status(dvar_runoff, "Memory allocation error.");
        for (i = 0; i < global_domain.ncells_total; i++) {
            dvar_runoff[i] = 0.0;
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
    status = MPI_Gatherv(var_local_runoff, local_domain.ncells_active, MPI_DOUBLE,
            dvar_gathered, mpi_map_local_array_sizes,
            mpi_map_global_array_offsets, MPI_DOUBLE,
            VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");
    free(var_local_runoff);

    if (mpi_rank == VIC_MPI_ROOT) {
        map(sizeof (double), global_domain.ncells_active, NULL,
                mpi_map_mapping_array, dvar_gathered, dvar_remapped);
        // expand to full grid size
        map(sizeof (double), global_domain.ncells_active, NULL,
                filter_active_cells, dvar_remapped, dvar_runoff);

        free(dvar_gathered);
        free(dvar_remapped);
    }

    // Do the routing on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        debug("RVIC");

        dvar = malloc(global_domain.ncells_total * sizeof (*dvar));
        check_alloc_status(dvar, "Memory allocation error.");
        for (i = 0; i < global_domain.ncells_total; i++) {
            dvar[i] = 0.0;
        }

        dvar_filtered =
                malloc(global_domain.ncells_active * sizeof (*dvar_filtered));
        check_alloc_status(dvar_filtered, "Memory allocation error.");

        dvar_mapped =
                malloc(global_domain.ncells_active * sizeof (*dvar_mapped));
        check_alloc_status(dvar_mapped, "Memory allocation error.");

        convolution(dvar_runoff, dvar);

        // filter the active cells only
        map(sizeof (double), global_domain.ncells_active, filter_active_cells,
                NULL, dvar, dvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof (double), global_domain.ncells_active, mpi_map_mapping_array,
                NULL, dvar_filtered, dvar_mapped);

        free(dvar);
        free(dvar_runoff);
        free(dvar_filtered);
    }
    var2 = malloc(local_domain.ncells_active * sizeof (*var2));
    check_alloc_status(var2, "Memory allocation error.");

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(dvar_mapped, mpi_map_local_array_sizes,
            mpi_map_global_array_offsets, MPI_DOUBLE,
            var2, local_domain.ncells_active, MPI_DOUBLE,
            VIC_MPI_ROOT, MPI_COMM_VIC);
    check_mpi_status(status, "MPI error.");

    if (mpi_rank == VIC_MPI_ROOT) {
        free(dvar_mapped);
    }

    // Write to output struct...
    for (i = 0; i < local_domain.ncells_active; i++) {
        out_data[i][OUT_DISCHARGE][0] = var2[i];
    }
    free(var2);
}
