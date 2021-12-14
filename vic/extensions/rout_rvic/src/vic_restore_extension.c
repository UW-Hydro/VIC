/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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
