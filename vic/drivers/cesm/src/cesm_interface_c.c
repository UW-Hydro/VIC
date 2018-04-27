/******************************************************************************
 * @section DESCRIPTION
 *
 * C interface for CESM driver.
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

#include <vic_driver_cesm.h>

size_t              NF, NR;
size_t              current;
size_t             *filter_active_cells = NULL;
size_t             *mpi_map_mapping_array = NULL;
all_vars_struct    *all_vars = NULL;
force_data_struct  *force = NULL;
x2l_data_struct    *x2l_vic = NULL;
l2x_data_struct    *l2x_vic = NULL;
dmy_struct          dmy_current;
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
MPI_Datatype        mpi_alarm_struct_type;
MPI_Datatype        mpi_option_struct_type;
MPI_Datatype        mpi_param_struct_type;
int                *mpi_map_local_array_sizes = NULL;
int                *mpi_map_global_array_offsets = NULL;
int                 mpi_rank;
int                 mpi_size;
option_struct       options;
parameters_struct   param;
param_set_struct    param_set;
soil_con_struct    *soil_con = NULL;
veg_con_map_struct *veg_con_map = NULL;
veg_con_struct    **veg_con = NULL;
veg_hist_struct   **veg_hist = NULL;
veg_lib_struct    **veg_lib = NULL;
metadata_struct     state_metadata[N_STATE_VARS];
metadata_struct     out_metadata[N_OUTVAR_TYPES];
save_data_struct   *save_data;  // [ncells]
double           ***out_data = NULL;  // [ncells, nvars, nelem]
stream_struct      *output_streams = NULL;  // [nstreams]
nc_file_struct     *nc_hist_files = NULL;  // [nstreams]
timer_struct        global_timers[N_TIMERS];

/******************************************************************************
 * @brief    Initialization function for CESM driver
 *****************************************************************************/
int
vic_cesm_init(vic_clock     *vclock,
              case_metadata *cmeta)
{
    // start vic all timer
    timer_start(&(global_timers[TIMER_VIC_ALL]));
    // start vic init timer
    timer_start(&(global_timers[TIMER_VIC_INIT]));

    // read global parameters
    vic_cesm_start(vclock, cmeta);

    // Initialize time
    initialize_cesm_time();

    // allocate memory
    vic_cesm_alloc();

    // initialize model parameters from parameter files
    vic_init();

    // populate model state, either using a cold start or from a restart file
    vic_populate_model_state(trimstr(cmeta->starttype), &dmy_current);

    // initialize forcings
    timer_start(&(global_timers[TIMER_VIC_FORCE]));
    vic_force();
    timer_stop(&(global_timers[TIMER_VIC_FORCE]));

    // initialize output structures
    vic_init_output(&dmy_current);

    // initialize albedo
    vic_initialize_albedo();

    // initialize temperature
    vic_initialize_temperature();

    // initialize upwelling longwave
    vic_initialize_lwup();

    // initialization is complete, print settings
    log_info(
        "Initialization is complete, print global param, parameters and options structures");
    print_global_param(&global_param);
    print_option(&options);
    print_parameters(&param);

    // stop init timer
    timer_stop(&(global_timers[TIMER_VIC_INIT]));
    // stop vic all timer
    timer_stop(&(global_timers[TIMER_VIC_ALL]));
    // init vic run timer
    timer_init(&(global_timers[TIMER_VIC_RUN]));

    return EXIT_SUCCESS;
}

/******************************************************************************
 * @brief    Run function for CESM driver
 *****************************************************************************/
int
vic_cesm_run(vic_clock *vclock)
{
    char state_filename[MAXSTRING];

    // continue vic all timer
    timer_continue(&(global_timers[TIMER_VIC_ALL]));
    // continue vic run timer
    timer_continue(&(global_timers[TIMER_VIC_RUN]));

    // reset l2x fields
    initialize_l2x_data();

    // read forcing data
    timer_continue(&(global_timers[TIMER_VIC_FORCE]));
    vic_force();
    timer_stop(&(global_timers[TIMER_VIC_FORCE]));

    // run vic over the domain
    vic_image_run(&dmy_current);

    // return fields to coupler
    vic_cesm_put_data();

    // Write history files
    timer_continue(&(global_timers[TIMER_VIC_WRITE]));
    vic_write_output(&dmy_current);
    timer_stop(&(global_timers[TIMER_VIC_WRITE]));

    // advance the clock
    advance_vic_time();
    assert_time_insync(vclock, &dmy_current);

    // if save:
    if (vclock->state_flag) {
        // write state file
        debug("writing state file for timestep %zu", current);
        vic_store(&dmy_current, state_filename);
        write_rpointer_file(state_filename);
        debug("finished storing state file: %s", state_filename)
    }

    // reset x2l fields
    initialize_x2l_data();

    // stop vic run timer
    timer_stop(&(global_timers[TIMER_VIC_RUN]));
    // stop vic all timer
    timer_stop(&(global_timers[TIMER_VIC_ALL]));

    return EXIT_SUCCESS;
}

/******************************************************************************
 * @brief    Finalize function for CESM driver
 *****************************************************************************/
int
vic_cesm_final(vic_clock *vclock)
{
    // continue vic all timer
    timer_continue(&(global_timers[TIMER_VIC_ALL]));
    // start vic final timer
    timer_start(&(global_timers[TIMER_VIC_FINAL]));

    // finalize time
    finalize_cesm_time(vclock);

    // clean up
    vic_cesm_finalize();

    // stop vic final timer
    timer_stop(&(global_timers[TIMER_VIC_FINAL]));
    // stop vic all timer
    timer_stop(&(global_timers[TIMER_VIC_ALL]));

    if (mpi_rank == VIC_MPI_ROOT) {
        // write timing info
        write_vic_timing_table(global_timers, VIC_DRIVER);
    }
    return EXIT_SUCCESS;
}
