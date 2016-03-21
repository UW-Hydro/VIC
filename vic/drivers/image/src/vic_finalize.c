/******************************************************************************
 * @section DESCRIPTION
 *
 * Finalize VIC run by freeing memory and closing open files.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/
void
vic_finalize(void)
{
    extern size_t             *filter_active_cells;
    extern size_t             *mpi_map_mapping_array;
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern dmy_struct         *dmy;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filep_struct        filep;
    extern int                *mpi_map_local_array_sizes;
    extern int                *mpi_map_global_array_offsets;
    extern int                 mpi_rank;
    extern nc_file_struct      nc_hist_file;
    extern option_struct       options;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;

    size_t                     i;
    size_t                     j;
    int                        status;

    if (mpi_rank == 0) {
        // close the global parameter file
        fclose(filep.globalparam);

        // close the netcdf history file if it is still open
        if (nc_hist_file.open == true) {
            status = nc_close(nc_hist_file.nc_id);
            if (status != NC_NOERR) {
                log_err("Error history file");
            }
        }
    }

    for (i = 0; i < local_domain.ncells_active; i++) {
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
            free_veg_hist(&(veg_hist[i][j]));
        }
        free_all_vars(&(all_vars[i]), veg_con[i][0].vegetat_type_num);
        free_out_data(&(out_data[i]));
        free(veg_con_map[i].vidx);
        free(veg_con_map[i].Cv);
        free(veg_con[i]);
        free(veg_hist[i]);
        free(veg_lib[i]);
    }
    free(atmos);
    free(soil_con);
    free(veg_con_map);
    free(veg_con);
    free(veg_hist);
    free(veg_lib);
    free(all_vars);
    free(out_data);
    free(save_data);
    free(local_domain.locations);
    free(dmy);
    if (mpi_rank == 0) {
        free(filter_active_cells);
        free(global_domain.locations);
        free(mpi_map_local_array_sizes);
        free(mpi_map_global_array_offsets);
        free(mpi_map_mapping_array);
    }
    finalize_logging();
}
