#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_finalize()
{
    extern all_vars_struct    *all_vars;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    size_t                     i;
    size_t                     j;

    for (i = 0; i < global_domain.ncells_global; i++) {
        free(soil_con[i].AreaFract);
        free(soil_con[i].BandElev);
        free(soil_con[i].Tfactor);
        free(soil_con[i].Pfactor);
        free(soil_con[i].AboveTreeLine);
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            free(veg_con[i][j].zone_depth);
            free(veg_con[i][j].zone_fract);
            if (options.CARBON) {
                free(veg_con[i][j].CanopLayerBnd);
            }
        }
        free(veg_con_map[i].vidx);
        free(veg_con_map[i].Cv);
        free(veg_con[i]);
        free(veg_lib[i]);
        free_all_vars(&(all_vars[i]), veg_con[i][0].vegetat_type_num);
    }
    free(soil_con);
    free(veg_con_map);
    free(veg_con);
    free(veg_lib);
    free(all_vars);
    free(global_domain.locations);
    free(dmy);
}