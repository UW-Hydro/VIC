#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

veg_var_struct **
make_veg_var(int veg_type_num)

/**********************************************************************
        make_veg_var	Dag Lohman		January 1996

   This routine makes an array of vegitation variables for each
   vegetation type.

**********************************************************************/
{
    extern option_struct options;

    int                  i, j;
    veg_var_struct     **temp = NULL;

    temp = (veg_var_struct **) calloc(veg_type_num, sizeof(veg_var_struct *));
    if (temp == NULL) {
        nrerror("Memory allocation error in make_veg_var().");
    }

    for (i = 0; i < veg_type_num; i++) {
        temp[i] =
            (veg_var_struct *) calloc(options.SNOW_BAND,
                                      sizeof(veg_var_struct));
        if (temp[i] == NULL) {
            nrerror("Memory allocation error in make_veg_var().");
        }

        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                temp[i][j].NscaleFactor = (double *) calloc(options.Ncanopy,
                                                            sizeof(double));
                if (temp[i][j].NscaleFactor == NULL) {
                    nrerror("Memory allocation error in make_veg_var().");
                }
                temp[i][j].aPARLayer =
                    (double *) calloc(options.Ncanopy, sizeof(double));
                if (temp[i][j].aPARLayer == NULL) {
                    nrerror("Memory allocation error in make_veg_var().");
                }
                temp[i][j].CiLayer =
                    (double *) calloc(options.Ncanopy, sizeof(double));
                if (temp[i][j].CiLayer == NULL) {
                    nrerror("Memory allocation error in make_veg_var().");
                }
                temp[i][j].rsLayer =
                    (double *) calloc(options.Ncanopy, sizeof(double));
                if (temp[i][j].rsLayer == NULL) {
                    nrerror("Memory allocation error in make_veg_var().");
                }
            }
        }
    }

    return temp;
}
