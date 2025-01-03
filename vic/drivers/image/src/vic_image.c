/******************************************************************************
 * @section DESCRIPTION
 *
 * Stand-alone image mode driver of the VIC model
 *****************************************************************************/

#define _MAIN_
#include <vic_driver_image.h>
#include <rout.h>   // Routing routine (extension)

size_t              NF, NR;
size_t              current;
size_t             *filter_active_cells = NULL;
size_t             *mpi_map_mapping_array = NULL;
all_vars_struct    *all_vars = NULL;
force_data_struct  *force = NULL;
dmy_struct         *dmy = NULL;
dmy_struct          dmy_state;
filenames_struct    filenames;
filep_struct        filep;
domain_struct       global_domain;
global_param_struct global_param;
lake_con_struct    *lake_con = NULL;
domain_struct       local_domain;
MPI_Comm            MPI_COMM_VIC = MPI_COMM_WORLD;
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
metadata_struct     state_metadata[N_STATE_VARS + N_STATE_VARS_EXT];
metadata_struct     out_metadata[N_OUTVAR_TYPES];
save_data_struct   *save_data;  // [ncells]
double           ***out_data = NULL;  // [ncells, nvars, nelem]
stream_struct      *output_streams = NULL;  // [nstreams]
nc_file_struct     *nc_hist_files = NULL;  // [nstreams]

// Extensions
rout_struct         rout; // Routing routine (extension)

/******************************************************************************
 * @brief   Stand-alone image mode driver of the VIC model
 * @details The image mode driver runs VIC for a single timestep for all grid
 *          cells before moving on to the next timestep.
 *
 * @param argc Argument count
 * @param argv Argument vector
 *****************************************************************************/
int
main(int    argc,
     char **argv)
{
    int          status;
    timer_struct global_timers[N_TIMERS];
    char         state_filename[MAXSTRING];

    // start vic all timer
    timer_start(&(global_timers[TIMER_VIC_ALL]));
    // start vic init timer
    timer_start(&(global_timers[TIMER_VIC_INIT]));

    // Initialize MPI - note: logging not yet initialized
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }

    // Initialize Log Destination
    initialize_log();

    // initialize mpi
    initialize_mpi();

    // process command line arguments
    if (mpi_rank == VIC_MPI_ROOT) {
        cmd_proc(argc, argv, filenames.global);
    }

    // read global parameters
    vic_image_start();

    // allocate memory
    vic_alloc();

    // allocate memory for routing
    rout_alloc();   // Routing routine (extension)

    // initialize model parameters from parameter files
    vic_image_init();

    // initialize routing parameters from parameter files
    rout_init();    // Routing routine (extension)

    // populate model state, either using a cold start or from a restart file
    vic_populate_model_state(&(dmy[0]));

    // initialize output structures
    vic_init_output(&(dmy[0]));

    // Initialization is complete, print settings
    log_info(
        "Initialization is complete, print global param, parameters and options structures");
    print_global_param(&global_param);
    print_option(&options);
    print_parameters(&param);

    // stop init timer
    timer_stop(&(global_timers[TIMER_VIC_INIT]));
    // start vic run timer
    timer_start(&(global_timers[TIMER_VIC_RUN]));

    // loop over all timesteps
    for (current = 0; current < global_param.nrecs; current++) {
        // read forcing data
        timer_continue(&(global_timers[TIMER_VIC_FORCE]));
        vic_force();
        timer_stop(&(global_timers[TIMER_VIC_FORCE]));

        // run vic over the domain
        vic_image_run(&(dmy[current]));

        // Write history files
        timer_continue(&(global_timers[TIMER_VIC_WRITE]));
        vic_write_output(&(dmy[current]));
        timer_stop(&(global_timers[TIMER_VIC_WRITE]));

        // Write state file
        if (check_save_state_flag(current, &dmy_state)) {
            debug("writing state file for timestep %zu", current);
            vic_store(&dmy_state, state_filename);
            debug("finished storing state file: %s", state_filename)
        }
    }
    // stop vic run timer
    timer_stop(&(global_timers[TIMER_VIC_RUN]));
    // start vic final timer
    timer_start(&(global_timers[TIMER_VIC_FINAL]));
    // clean up
    vic_image_finalize();

    // clean up routing
    rout_finalize();    // Routing routine (extension)

    // finalize MPI
    status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
        log_err("MPI error: %d", status);
    }

    log_info("Completed running VIC %s", VIC_DRIVER);

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
