/******************************************************************************
 * @section DESCRIPTION
 *
 * Finalize VIC run by freeing memory and closing open files.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Finalize VIC run by freeing memory and closing open files.
 *****************************************************************************/
void
vic_finalize(void)
{
    extern size_t             *filter_active_cells;
    extern size_t             *mpi_map_mapping_array;
    extern all_vars_struct    *all_vars;
    extern force_data_struct  *force;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filep_struct        filep;
    extern int                *mpi_map_local_array_sizes;
    extern int                *mpi_map_global_array_offsets;
    extern int                 mpi_rank;
    extern nc_file_struct     *nc_hist_files;
    extern option_struct       options;
    extern double           ***out_data;
    extern stream_struct      *output_streams;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;
    extern MPI_Datatype        mpi_global_struct_type;
    extern MPI_Datatype        mpi_filenames_struct_type;
    extern MPI_Datatype        mpi_location_struct_type;
    extern MPI_Datatype        mpi_alarm_struct_type;
    extern MPI_Datatype        mpi_option_struct_type;
    extern MPI_Datatype        mpi_param_struct_type;

    size_t                     i;
    size_t                     j;
    int                        status;


    if (mpi_rank == VIC_MPI_ROOT) {
        // close the global parameter file
        fclose(filep.globalparam);

        // close the netcdf history file if it is still open
        for (i = 0; i < options.Noutstreams; i++) {
            if (nc_hist_files[i].open == true) {
                status = nc_close(nc_hist_files[i].nc_id);
                check_nc_status(status, "Error closing history file");
            }
            free(nc_hist_files[i].nc_vars);
        }
        free(nc_hist_files);
    }

    for (i = 0; i < local_domain.ncells_active; i++) {
        free_force(&(force[i]));
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
        free_all_vars(&(all_vars[i]), veg_con_map[i].nv_active - 1);
        free(veg_con_map[i].vidx);
        free(veg_con_map[i].Cv);
        free(veg_con[i]);
        free(veg_hist[i]);
        free(veg_lib[i]);
    }

    free_streams(&output_streams);
    free_out_data(local_domain.ncells_active, out_data);
    free(force);
    free(soil_con);
    free(veg_con_map);
    free(veg_con);
    free(veg_hist);
    free(veg_lib);
    free(all_vars);
    free(save_data);
    free(local_domain.locations);
    if (mpi_rank == VIC_MPI_ROOT) {
        free(filter_active_cells);
        free(global_domain.locations);
        free(mpi_map_local_array_sizes);
        free(mpi_map_global_array_offsets);
        free(mpi_map_mapping_array);
    }

    MPI_Type_free(&mpi_global_struct_type);
    MPI_Type_free(&mpi_filenames_struct_type);
    MPI_Type_free(&mpi_location_struct_type);
    MPI_Type_free(&mpi_alarm_struct_type);
    MPI_Type_free(&mpi_option_struct_type);
    MPI_Type_free(&mpi_param_struct_type);
    finalize_logging();
}
