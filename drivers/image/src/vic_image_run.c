#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_image_run(void)
{
    extern size_t              current;
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern global_param_struct global_param;
    extern lake_con_struct     lake_con;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;

    size_t                     i;

    for (i = 0; i < global_domain.ncells_global; i++) {
        vic_run(current, &(atmos[i]), &(all_vars[i]), dmy, &global_param,
                &lake_con, &(soil_con[i]), veg_con[i], veg_lib[i], veg_hist[i]);
        put_data(&(all_vars[i]), &(atmos[i]), &(soil_con[i]), veg_con[i],
                 veg_lib[i], &lake_con, out_data[i], &(save_data[i]),
                 current);
    }
}
