#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

snow_data_struct **
make_snow_data(size_t nveg)

/**********************************************************************
        make_snow_data	Keith Cherkauer		January 22, 1997

   This routine makes an array of snow cover data structures, one
   for each vegetation type plus bare soil.

   modifications:
   07-09-98 modified to make te make a two dimensional array which
           also accounts for a variable number of snow elevation
           bands                                               KAC

**********************************************************************/
{
    extern option_struct options;

    size_t               i;
    snow_data_struct   **temp = NULL;

    temp = (snow_data_struct **) calloc(nveg,
                                        sizeof(snow_data_struct *));
    if (temp == NULL) {
        nrerror("Memory allocation error in make_snow_data().");
    }


    for (i = 0; i < nveg; i++) {
        temp[i] = (snow_data_struct *) calloc(options.SNOW_BAND,
                                              sizeof(snow_data_struct));
        if (temp[i] == NULL) {
            nrerror("Memory allocation error in make_snow_data().");
        }
    }

    return temp;
}
