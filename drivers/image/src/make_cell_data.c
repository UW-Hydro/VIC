#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

cell_data_struct **
make_cell_data(int veg_type_num,
               int Nlayer)

/**********************************************************************
        make_cell_data	Keith Cherkauer		July 9, 1997

   This subroutine makes an array of type cell, which contains soil
   column variables for a single grid cell.

**********************************************************************/
{
    extern option_struct options;

    int                  i;
    cell_data_struct   **temp = NULL;

    temp = (cell_data_struct**) calloc(veg_type_num,
                                       sizeof(cell_data_struct*));
    if (temp == NULL) {
        nrerror("Memory allocation error in make_cell_data().");
    }

    for (i = 0; i < veg_type_num; i++) {
        temp[i] = (cell_data_struct*) calloc(options.SNOW_BAND,
                                             sizeof(cell_data_struct));
        if (temp[i] == NULL) {
            nrerror("Memory allocation error in make_cell_data().");
        }

    }
    return temp;
}
