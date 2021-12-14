/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of vegitation variables for each vegitation type.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of vegitation variables for each vegitation type.
 *****************************************************************************/
veg_var_struct **
make_veg_var(size_t veg_type_num)
{
    extern option_struct options;

    size_t               i, j;
    veg_var_struct     **temp = NULL;

    temp = calloc(veg_type_num, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    for (i = 0; i < veg_type_num; i++) {
        temp[i] = calloc(options.SNOW_BAND, sizeof(*(temp[i])));
        check_alloc_status(temp[i], "Memory allocation error.");

        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                temp[i][j].NscaleFactor = calloc(options.Ncanopy,
                                                 sizeof(*(temp[i][j].
                                                          NscaleFactor)));
                check_alloc_status(temp[i][j].NscaleFactor,
                                   "Memory allocation error.");
                temp[i][j].aPARLayer = calloc(options.Ncanopy,
                                              sizeof(*(temp[i][j].aPARLayer)));
                check_alloc_status(temp[i][j].aPARLayer,
                                   "Memory allocation error.");
                temp[i][j].CiLayer = calloc(options.Ncanopy,
                                            sizeof(*(temp[i][j].CiLayer)));
                check_alloc_status(temp[i][j].CiLayer,
                                   "Memory allocation error.");
                temp[i][j].rsLayer = calloc(options.Ncanopy,
                                            sizeof(*(temp[i][j].rsLayer)));
                check_alloc_status(temp[i][j].rsLayer,
                                   "Memory allocation error.");
            }
        }
    }

    return temp;
}
