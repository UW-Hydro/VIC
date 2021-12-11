/******************************************************************************
 * @section DESCRIPTION
 *
 * Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/
void
vic_cesm_finalize(void)
{
    extern x2l_data_struct *x2l_vic;
    extern l2x_data_struct *l2x_vic;

    // free VIC/CESM data structures
    free(x2l_vic);
    free(l2x_vic);

    vic_finalize();
}
