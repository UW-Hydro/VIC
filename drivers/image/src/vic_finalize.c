#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_finalize()
{
    extern domain_struct    global_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern veg_con_struct **veg_con;
    extern veg_lib_struct **veg_lib;

    size_t                  i;
    size_t                  j;

    for (i = 0; i < global_domain.ncells_global; i++) {
        printf("%zd\n", i);
        free(soil_con[i].AreaFract);
        free(soil_con[i].BandElev);
        free(soil_con[i].Tfactor);
        free(soil_con[i].Pfactor);
        free(soil_con[i].AboveTreeLine);
        for (j = 0; j < options.NVEGTYPES; j++) {
            free(veg_con[i][j].zone_depth);
            free(veg_con[i][j].zone_fract);
            if (options.CARBON) {
                free(veg_con[i][j].CanopLayerBnd);
            }
        }
        free(veg_con[i]);
        free(veg_lib[i]);
    }
    free(soil_con);
    free(veg_con);
    free(veg_lib);
    free(global_domain.locations);
}
