#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
free_all_vars(all_vars_struct *all_vars,
              int              Nveg)

/**********************************************************************
        free_all_vars	Keith Cherkauer		March 1998

   This routine frees all memory allocated down the all_vars
   data structure.  This include all grid cell specific
   variables (soil, vegetation, energy, snow).

**********************************************************************/
{
    extern option_struct options;

    int                  i, j, k, Nitems;

    Nitems = Nveg + 1;

    for (j = 0; j < Nitems; j++) {
        free((char *) all_vars[0].cell[j]);
    }
    free((char *) all_vars[0].cell);
    for (j = 0; j < Nitems; j++) {
        if (options.CARBON) {
            for (k = 0; k < options.SNOW_BAND; k++) {
                free((char *) all_vars[0].veg_var[j][k].NscaleFactor);
                free((char *) all_vars[0].veg_var[j][k].aPARLayer);
                free((char *) all_vars[0].veg_var[j][k].CiLayer);
                free((char *) all_vars[0].veg_var[j][k].rsLayer);
            }
        }
        free((char *) (*all_vars).veg_var[j]);
    }
    free((char *) (*all_vars).veg_var);
    for (j = 0; j < Nitems; j++) {
        free((char *) all_vars[0].energy[j]);
    }
    free((char *) all_vars[0].energy);
    for (i = 0; i < Nitems; i++) {
        free((char *) all_vars[0].snow[i]);
    }
    free((char *) all_vars[0].snow);
}
