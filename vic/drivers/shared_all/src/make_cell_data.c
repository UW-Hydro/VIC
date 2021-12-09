/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine makes an array of type cell, which contains soil column
 * variables for a single grid cell.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Make an array of type cell, which contains soil column variables
 *           for a single grid cell.
 *****************************************************************************/
cell_data_struct **
make_cell_data(size_t veg_type_num)
{
    extern option_struct options;

    size_t               i;
    cell_data_struct   **temp;

    temp = calloc(veg_type_num, sizeof(*temp));
    for (i = 0; i < veg_type_num; i++) {
        temp[i] = calloc(options.SNOW_BAND, sizeof(*(temp[i])));
    }
    return temp;
}
