/******************************************************************************
 * @section DESCRIPTION
 *
 * C interface for CESM driver.
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

#include <vic_driver_cesm.h>

size_t              NF, NR;
size_t              current;
size_t             *filter_active_cells = NULL;
size_t             *mpi_map_mapping_array = NULL;
all_vars_struct    *all_vars = NULL;
atmos_data_struct  *atmos = NULL;
x2l_data_struct    *x2l_vic = NULL;
l2x_data_struct    *l2x_vic = NULL;
dmy_struct          dmy;
filenames_struct    filenames;
filep_struct        filep;
domain_struct       global_domain;
domain_struct       local_domain;
global_param_struct global_param;
lake_con_struct     lake_con;
MPI_Comm            MPI_COMM_VIC;
MPI_Datatype        mpi_domain_struct_type;
MPI_Datatype        mpi_global_struct_type;
MPI_Datatype        mpi_filenames_struct_type;
MPI_Datatype        mpi_location_struct_type;
MPI_Datatype        mpi_nc_file_struct_type;
MPI_Datatype        mpi_option_struct_type;
MPI_Datatype        mpi_param_struct_type;
int                *mpi_map_local_array_sizes = NULL;
int                *mpi_map_global_array_offsets = NULL;
int                 mpi_rank;
int                 mpi_size;
nc_file_struct      nc_hist_file;
nc_var_struct       nc_vars[N_OUTVAR_TYPES];
option_struct       options;
parameters_struct   param;
out_data_struct   **out_data;
save_data_struct   *save_data;
param_set_struct    param_set;
soil_con_struct    *soil_con = NULL;
veg_con_map_struct *veg_con_map = NULL;
veg_con_struct    **veg_con = NULL;
veg_hist_struct   **veg_hist = NULL;
veg_lib_struct    **veg_lib = NULL;

/******************************************************************************
 * @brief    Initialization function for CESM driver
 *****************************************************************************/
int
vic_cesm_init(vic_clock     *vclock,
              case_metadata *cmeta)
{
    // read global parameters
    vic_cesm_start(vclock, cmeta);

    // Initialize time
    initialize_cesm_time();

    // allocate memory
    vic_cesm_alloc();

    // initialize model parameters from parameter files
    vic_init();

    // restore model state, either using a cold start or from a restart file
    vic_restore(trim(cmeta->starttype));

    // initialize output structures
    vic_init_output();

    return 0;
}

/******************************************************************************
 * @brief    Run function for CESM driver
 *****************************************************************************/
int
vic_cesm_run(vic_clock *vclock)
{

    // reset l2x fields
    initialize_l2x_data();

    // advance the clock
    advance_time();
    assert_time_insync(vclock, &dmy);

    // read forcing data
    vic_force();

    // run vic over the domain
    vic_image_run();

    // return fields to coupler
    vic_cesm_put_data();

    // if output:
    if (check_write_flag(current)) {
        vic_write();
    }

    // if save:
    if (vclock->state_flag) {
        vic_store();
    }

    // reset x2l fields
    initialize_x2l_data();

    return 0;
}

/******************************************************************************
 * @brief    Finalize function for CESM driver
 *****************************************************************************/
int
vic_cesm_final()
{
    // clean up
    vic_cesm_finalize();

    return 0;
}
