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
state_metadata_rout_extension()
{
    extern metadata_struct state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];

    // STATE_ROUT_RING
    strcpy(state_metadata[N_STATE_VARS + STATE_ROUT_RING].varname,
           "STATE_ROUT_RING");
    strcpy(state_metadata[N_STATE_VARS + STATE_ROUT_RING].long_name,
           "routing_ring");
    strcpy(state_metadata[N_STATE_VARS + STATE_ROUT_RING].standard_name,
           "routing_ring");
    strcpy(state_metadata[N_STATE_VARS + STATE_ROUT_RING].units, "-");
    strcpy(state_metadata[N_STATE_VARS + STATE_ROUT_RING].description,
           "unit hydrographs in the routing ring");
}
