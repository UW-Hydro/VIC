/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine opens a model state file and verifys that the starting date,
 * number of layers and number of thermal nodes in the file agrees with what
 * was defined in the model global control file.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This subroutine opens a model state file and verifys that the
             starting date, number of layers and number of thermal nodes in the
             file agrees with what was defined in the model global control file.
 *****************************************************************************/
FILE *
check_state_file(char  *init_state_name,
                 size_t Nlayer,
                 size_t Nnodes,
                 int   *startrec)
{
    extern option_struct options;

    FILE                *init_state;
    size_t               tmp_Nlayer;
    size_t               tmp_Nnodes;
    unsigned short int   startday, startmonth, startyear;

    /* open state file */
    if (options.STATE_FORMAT == BINARY) {
        init_state = open_file(init_state_name, "rb");
    }
    else {
        init_state = open_file(init_state_name, "r");
    }

    /* Initialize startrec */
    *startrec = 0;

    /* Check state date information */
    if (options.STATE_FORMAT == BINARY) {
        fread(&startyear, sizeof(int), 1, init_state);
        fread(&startmonth, sizeof(int), 1, init_state);
        fread(&startday, sizeof(int), 1, init_state);
    }
    else {
        fscanf(init_state, "%hu %hu %hu\n", &startyear, &startmonth, &startday);
    }

    /* Check simulation options */
    if (options.STATE_FORMAT == BINARY) {
        fread(&tmp_Nlayer, sizeof(size_t), 1, init_state);
        fread(&tmp_Nnodes, sizeof(size_t), 1, init_state);
    }
    else {
        fscanf(init_state, "%zu %zu\n", &tmp_Nlayer, &tmp_Nnodes);
    }
    if (tmp_Nlayer != Nlayer) {
        log_err("The number of soil moisture layers in the model state file "
                "(%zu) does not equal that defined in the global control file "
                "(%zu).  Check your input files.", tmp_Nlayer, Nlayer);
    }
    if (tmp_Nnodes != Nnodes) {
        log_err("The number of soil thermal nodes in the model state file "
                "(%zu) does not equal that defined in the global control file "
                "(%zu).  Check your input files.", tmp_Nnodes, Nnodes);
    }

    return(init_state);
}
