/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine opens files for soil, vegetation, and global parameters.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine opens files for soil, vegetation, and global
 *           parameters.
 *****************************************************************************/
void
check_files(filep_struct     *filep,
            filenames_struct *fnames)
{
    extern option_struct options;
    extern FILE          *open_file(char string[], char type[]);

    filep->soilparam = open_file(fnames->soil, "r");
    filep->veglib = open_file(fnames->veglib, "r");
    filep->vegparam = open_file(fnames->veg, "r");
    if (options.SNOW_BAND > 1) {
        filep->snowband = open_file(fnames->snowband, "r");
    }
    if (options.LAKES) {
        filep->lakeparam = open_file(fnames->lakeparam, "r");
    }
}
