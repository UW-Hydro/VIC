/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes all filefilenames before they are called by
 * the model.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Initialize all filenames before they are called by the
 *           model.
 *****************************************************************************/
void
initialize_filenames()
{
    extern filenames_struct filenames;

    size_t                  i;

    strcpy(filenames.init_state.nc_filename, "MISSING");
    strcpy(filenames.statefile, "MISSING");
    strcpy(filenames.constants, "MISSING");
    strcpy(filenames.params.nc_filename, "MISSING");
    strcpy(filenames.rout_params.nc_filename, "MISSING");
    strcpy(filenames.domain.nc_filename, "MISSING");
    strcpy(filenames.result_dir, "MISSING");
    strcpy(filenames.log_path, "MISSING");
    for (i = 0; i < 2; i++) {
        strcpy(filenames.f_path_pfx[i], "MISSING");
    }
}

/******************************************************************************
 * @brief    Initialize all file pointers
 *****************************************************************************/
void
initialize_fileps()
{
    extern filep_struct filep;

    filep.globalparam = NULL;
    filep.constants = NULL;
    filep.logfile = NULL;
}
