/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine makes an array of frozen soil data structures, one for each
 * vegetation type and bare soil.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine makes an array of frozen soil data structures, one
 *           for each vegetation type and bare soil.
 *****************************************************************************/
energy_bal_struct **
make_energy_bal(size_t nveg)
{
    extern option_struct options;

    size_t               i, j;
    energy_bal_struct  **temp = NULL;

    temp = calloc(nveg, sizeof(*temp));
    check_alloc_status(temp, "Memory allocation error.");

    /** Initialize all records to unfrozen conditions */
    for (i = 0; i < nveg; i++) {
        temp[i] = calloc(options.SNOW_BAND, sizeof(*(temp[i])));
        check_alloc_status(temp[i], "Memory allocation error.");

        for (j = 0; j < options.SNOW_BAND; j++) {
            temp[i][j].frozen = false;
        }
    }

    return temp;
}
