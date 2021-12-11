/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine frees all components of the veg_con structure.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This subroutine frees all components of the veg_con structure.
 *****************************************************************************/
void
free_vegcon(veg_con_struct **veg_con)
{
    extern option_struct options;
    size_t               i;

    for (i = 0; i < veg_con[0][0].vegetat_type_num; i++) {
        free((char *) veg_con[0][i].zone_depth);
        free((char *) veg_con[0][i].zone_fract);
        if (options.CARBON) {
            free((char *) veg_con[0][i].CanopLayerBnd);
        }
    }
    free((char *) veg_con[0]);
}
