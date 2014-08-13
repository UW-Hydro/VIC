#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

energy_bal_struct **
make_energy_bal(int nveg)

/**********************************************************************
        make_energy_bal	Keith Cherkauer		May 26, 1996

   This routine makes an array of frozen soil data structures, one
   for each vegetation type and bare soil.

   Modifications:
   01-Nov-04 Removed modification of Nnodes, as this was preventing
            correct reading/writing of state files for QUICK_FLUX
            =TRUE.						TJB

**********************************************************************/
{
    extern option_struct options;

    int                  i, j;
    energy_bal_struct  **temp = NULL;

    temp = (energy_bal_struct**) calloc(nveg,
                                        sizeof(energy_bal_struct*));
    if (temp == NULL) {
        nrerror("Memory allocation error in make_energy_bal().");
    }


    /** Initialize all records to unfrozen conditions */
    for (i = 0; i < nveg; i++) {
        temp[i] = (energy_bal_struct*) calloc(options.SNOW_BAND,
                                              sizeof(energy_bal_struct));
        if (temp[i] == NULL) {
            nrerror("Memory allocation error in make_energy_bal().");
        }
        for (j = 0; j < options.SNOW_BAND; j++) {
            temp[i][j].frozen = FALSE;
        }
    }

    return temp;
}
