#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_restore()
{
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    double surf_temp;
    size_t                     i;
    size_t        nveg;

    // read the model state from the netcdf file if there is one
    if (options.INIT_STATE) {
    }
    
    // run through the remaining VIC initialization routines
    for (i = 0; i < global_domain.ncells_global; i++) {
        // TBD: do something sensible for surf_temp
        surf_temp = 0.;
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_model_state(&(all_vars[i]), nveg, options.Nnode, 
                                surf_temp, &(soil_con[i]), veg_con[i],
                                veg_lib[i]);
    }
}
