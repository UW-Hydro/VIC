#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

all_vars_struct
make_all_vars(size_t nveg)

/**********************************************************************
        read_all_vars	Keith Cherkauer		May 21, 1996

   This routine creates an array of structures that contain information
   about a cell's states and fluxes.

**********************************************************************/
{
    all_vars_struct      temp;
    size_t               Nitems;

    Nitems = nveg + 1;

    temp.snow = make_snow_data(Nitems);
    temp.energy = make_energy_bal(Nitems);
    temp.veg_var = make_veg_var(Nitems);
    temp.cell = make_cell_data(Nitems);

    return (temp);
}
