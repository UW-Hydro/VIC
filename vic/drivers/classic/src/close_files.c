/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine closes all forcing data files, and output files.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine closes all forcing data files, and output files.
 *****************************************************************************/
void
close_files(filep_struct   *filep,
            stream_struct **streams)
{
    extern option_struct options;

    size_t               streamnum;

    /**********************
       Close All Input Files
    **********************/

    fclose(filep->forcing[0]);

    if (filep->forcing[1] != NULL) {
        fclose(filep->forcing[1]);
    }

    /*******************
       Close Output Files
    *******************/
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        fclose((*streams)[streamnum].fh);
        if ((*streams)[streamnum].compress) {
            compress_files((*streams)[streamnum].filename,
                           (*streams)[streamnum].compress);
        }
    }
}
