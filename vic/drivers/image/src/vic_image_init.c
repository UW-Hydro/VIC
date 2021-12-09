/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize model parameters for image driver.
 *****************************************************************************/

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Initialize model parameters
 *****************************************************************************/
void
vic_image_init(void)
{
    extern dmy_struct         *dmy;
    extern global_param_struct global_param;

    // make_dmy()
    initialize_time();
    dmy = make_dmy(&global_param);

    vic_init();
}
