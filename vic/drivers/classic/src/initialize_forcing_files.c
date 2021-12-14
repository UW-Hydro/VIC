/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes the forcing file parameters
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize Forcing File Parameters
 *****************************************************************************/
void
initialize_forcing_files()
{
    extern param_set_struct param_set;

    int                     i, j;

    /** Initialize forcing file input controls **/

    for (j = 0; j < N_FORCING_TYPES; j++) {
        param_set.TYPE[j].SUPPLIED = false;
        param_set.TYPE[j].SIGNED = 1;
        param_set.TYPE[j].multiplier = 1;
    }
    for (i = 0; i < 2; i++) {
        param_set.FORCE_DT[i] = MISSING;
        param_set.force_steps_per_day[i] = 0;
        param_set.N_TYPES[i] = 0;
        param_set.FORCE_FORMAT[i] = MISSING;
        for (j = 0; j < N_FORCING_TYPES; j++) {
            param_set.FORCE_INDEX[i][j] = MISSING;
        }
    }
}
