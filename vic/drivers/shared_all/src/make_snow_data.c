/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of snow cover data structures, one for each
 * vegetation type plus bare soil.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of snow cover data structures, one for each
 *           vegetation type plus bare soil.
 *****************************************************************************/
snow_data_struct **
make_snow_data(size_t nveg)
{
    extern option_struct options;

    size_t               i;
    snow_data_struct   **temp = NULL;

    temp = calloc(nveg, sizeof(*temp));

    check_alloc_status(temp, "Memory allocation error.");

    for (i = 0; i < nveg; i++) {
        temp[i] = calloc(options.SNOW_BAND, sizeof(*(temp[i])));
        check_alloc_status(temp[i], "Memory allocation error.");
    }

    return temp;
}
