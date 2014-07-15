#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(vic_finalize): %s\n", nc_strerror(e)); }

void
vic_finalize(void)
{
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern filep_struct        filep;
    extern nc_file_struct      nc_hist_file;    
    extern option_struct       options;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    size_t                     i;
    size_t                     j;
    int                        status;

    // close the global parameter file
    fclose(filep.globalparam);
    
    // close the netcdf history file if it is still open
    if (nc_hist_file.open == true) {
        status = nc_close(nc_hist_file.nc_id);
        if (status != NC_NOERR) {
            ERR(status);
        }
    }

    for (i = 0; i < global_domain.ncells_global; i++) {
        free_atmos(&(atmos[i]));
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
        free_out_data(&(out_data[i]));
    }
    free(atmos);
    free(soil_con);
    free(veg_con_map);
    free(veg_con);
    free(veg_lib);
    free(all_vars);
    free(out_data);
    free(save_data);
    free(global_domain.locations);
    free(dmy);
}
