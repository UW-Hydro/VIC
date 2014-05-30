#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_init_output()
{
    extern size_t              current;
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern lake_con_struct     lake_con;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;

    size_t                     i;

    for (i = 0; i < global_domain.ncells_global; i++) {
        put_data(&(all_vars[i]), &(atmos[i]), &(soil_con[i]), veg_con[i],
                 &lake_con, out_data[i], &(save_data[i]), &dmy[current],
                 -1);
    }
}
