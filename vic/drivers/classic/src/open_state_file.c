/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine opens the model state file for output.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Open state file to write to.
 *****************************************************************************/
FILE *
open_state_file(global_param_struct *global,
                filenames_struct     filenames,
                size_t               Nlayer,
                size_t               Nnodes)
{
    extern option_struct options;

    FILE                *statefile;
    char                 filename[MAXSTRING];

    /* open state file */
    sprintf(filename, "%s", filenames.statefile);
    if (options.STATE_FORMAT == BINARY) {
        statefile = open_file(filename, "wb");
    }
    else {
        statefile = open_file(filename, "w");
    }

    /* Write save state date information */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&global->stateyear, sizeof(int), 1, statefile);
        fwrite(&global->statemonth, sizeof(int), 1, statefile);
        fwrite(&global->stateday, sizeof(int), 1, statefile);
    }
    else {
        fprintf(statefile, "%i %i %i\n", global->stateyear,
                global->statemonth, global->stateday);
    }

    /* Write simulation flags */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&Nlayer, sizeof(size_t), 1, statefile);
        fwrite(&Nnodes, sizeof(size_t), 1, statefile);
    }
    else {
        fprintf(statefile, "%zu %zu\n", Nlayer, Nnodes);
    }

    return(statefile);
}
