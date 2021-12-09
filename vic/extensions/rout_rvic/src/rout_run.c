/******************************************************************************
 * @section DESCRIPTION
 *
 * Run routing over the domain.
 *****************************************************************************/

#include <rout.h>

/******************************************************************************
* @brief        This subroutine controls the RVIC convolution.
******************************************************************************/
void
rout_run(void)
{
    extern int           mpi_rank;
    extern double     ***out_data;
    extern domain_struct local_domain;
    extern domain_struct global_domain;
    double              *var_local_runoff = NULL;
    double              *var_local_discharge = NULL;
    double              *var_domain_runoff = NULL;
    double              *var_domain_discharge = NULL;
    size_t               i;

    debug("RVIC");

    // Allocate memory for the local_domain variables
    var_local_runoff =
        malloc(local_domain.ncells_active * sizeof(*var_local_runoff));
    check_alloc_status(var_local_runoff, "Memory allocation error.");

    var_local_discharge =
        malloc(local_domain.ncells_active * sizeof(*var_local_discharge));
    check_alloc_status(var_local_discharge, "Memory allocation error.");


    // Allocate memory for entire global_domain variables on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        var_domain_runoff =
            malloc(global_domain.ncells_total * sizeof(*var_domain_runoff));
        check_alloc_status(var_domain_runoff, "Memory allocation error.");

        var_domain_discharge =
            malloc(global_domain.ncells_total * sizeof(*var_domain_discharge));
        check_alloc_status(var_domain_discharge, "Memory allocation error.");

        // Initialize discharge to zero
        for (i = 0; i < global_domain.ncells_total; i++) {
            var_domain_discharge[i] = 0;
        }
    }

    // Read from runoff and baseflow from out_data and sum to runoff
    for (i = 0; i < local_domain.ncells_active; i++) {
        var_local_runoff[i] = out_data[i][OUT_RUNOFF][0] +
                              out_data[i][OUT_BASEFLOW][0];
    }

    // Gather the runoff for the local nodes
    gather_field_double(0.0, var_domain_runoff, var_local_runoff);

    // Run the convolution on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        convolution(var_domain_runoff, var_domain_discharge);
    }

    // Scatter the discharge back to the local nodes
    scatter_field_double(var_domain_discharge, var_local_discharge);

    // Write to output struct
    for (i = 0; i < local_domain.ncells_active; i++) {
        out_data[i][OUT_DISCHARGE][0] = var_local_discharge[i];
    }

    // Free variables on the local nodes
    free(var_local_runoff);
    free(var_local_discharge);

    // Free variables on the master node
    if (mpi_rank == VIC_MPI_ROOT) {
        free(var_domain_runoff);
    }
}
