/******************************************************************************
 * @section DESCRIPTION
 *
 * Top level wrapper for globals used in vic_run but defined at the global
 * level.
 *****************************************************************************/

#define _MAIN_
#include <vic_driver_python.h>

// global variables
int                 flag;
size_t              NR; /* array index for atmos struct that indicates
                           the model step avarage or sum */
size_t              NF; /* array index loop counter limit for atmos
                           struct that indicates the SNOW_STEP values */

global_param_struct global_param;
option_struct       options;
parameters_struct   param;
param_set_struct    param_set;
metadata_struct     out_metadata[N_OUTVAR_TYPES];
