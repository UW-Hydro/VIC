/******************************************************************************
 * @section DESCRIPTION
 *
 * Finalize VIC by freeing memory and closing open files.
 *****************************************************************************/

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/
void
vic_image_finalize(void)
{
    extern dmy_struct *dmy;

    // free data structures specific to to image driver
    free(dmy);

    vic_finalize();
}
