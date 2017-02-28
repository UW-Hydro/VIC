/******************************************************************************
 * @section DESCRIPTION
 *
 * Classic driver of the VIC model
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

#include <vic_driver_classic.h>
#include <cuda.h>

// global variables
int                 flag;
size_t              NR; /* array index for atmos struct that indicates
                           the model step avarage or sum */
size_t              NF; /* array index loop counter limit for atmos
                           struct that indicates the SNOW_STEP values */

global_param_struct global_param;
veg_lib_struct     *veg_lib;
option_struct       options;
Error_struct        Error;
param_set_struct    param_set;
parameters_struct   param;
filenames_struct    filenames;
filep_struct        filep;
metadata_struct     out_metadata[N_OUTVAR_TYPES];

/******************************************************************************
 * @brief   Classic driver of the VIC model
 * @details The classic driver runs VIC for a single grid cell for all
 *          timesteps before moving on to the next grid cell.
 *
 * @param argc Argument count
 * @param argv Argument vector
 *****************************************************************************/
int
main(int   argc,
     char *argv[])
{
    /** Variable Declarations **/
    extern FILE       *LOG_DEST;

    bool               MODEL_DONE;
    bool               RUN_MODEL;
    char               dmy_str[MAXSTRING];
    size_t             rec;
    size_t             Nveg_type;
    int                cellnum;
    int                startrec;
    int                ErrorFlag;
    int                n;
    size_t             streamnum;
    dmy_struct        *dmy;
    force_data_struct *force;
    veg_hist_struct  **veg_hist;
    veg_con_struct    *veg_con;
    soil_con_struct    soil_con;
    all_vars_struct    all_vars;
    lake_con_struct    lake_con;
    stream_struct     *streams = NULL;
    double          ***out_data;   // [1, nvars, nelem]
    save_data_struct   save_data;
    timer_struct       global_timers[N_TIMERS];
    timer_struct       cell_timer;

    // start vic all timer
    timer_start(&(global_timers[TIMER_VIC_ALL]));
    // start vic init timer
    timer_start(&(global_timers[TIMER_VIC_INIT]));

    // Initialize Log Destination
    initialize_log();

    /** Read Model Options **/
    cmd_proc(argc, argv, filenames.global);

    // Initialize global structures
    initialize_options();
    initialize_global();
    initialize_parameters();
    initialize_filenames();

    /* Initilize forcing file param structure */
    initialize_forcing_files();

    /** Read Global Control File **/
    filep.globalparam = open_file(filenames.global, "r");
    get_global_param(filep.globalparam);
    fclose(filep.globalparam);

    // Set Log Destination
    setup_logging(MISSING, filenames.log_path, &(filep.logfile));

    /** Set model constants **/
    if (strcmp(filenames.constants, "MISSING") != 0) {
        filep.constants = open_file(filenames.constants, "r");
        get_parameters(filep.constants);
    }
    // Check that model parameters are valid
    validate_parameters();

    /** Make Date Data Structure **/
    initialize_time();
    dmy = make_dmy(&global_param);

    /** Set up output data structures **/
    set_output_met_data_info();
    // out_data is shape [ngridcells (1), N_OUTVAR_TYPES]
    alloc_out_data(1, &out_data);
    filep.globalparam = open_file(filenames.global, "r");
    parse_output_info(filep.globalparam, &streams, &(dmy[0]));
    validate_streams(&streams);

    /** Check and Open Files **/
    check_files(&filep, &filenames);

    /** Read Vegetation Library File **/
    veg_lib = read_veglib(filep.veglib, &Nveg_type);

    /** Initialize Parameters **/
    cellnum = -1;

    /** allocate memory for the force_data_struct **/
    alloc_atmos(global_param.nrecs, &force);

    /** Initial state **/
    startrec = 0;
    if (options.INIT_STATE) {
        filep.init_state = check_state_file(filenames.init_state,
                                            options.Nlayer, options.Nnode,
                                            &startrec);
    }

    /** open state file if model state is to be saved **/
    if (options.SAVE_STATE && strcmp(filenames.statefile, "NONE") != 0) {
        filep.statefile = open_state_file(&global_param, filenames,
                                          options.Nlayer,
                                          options.Nnode);
    }
    else {
        filep.statefile = NULL;
    }

    /************************************
       Run Model for all Active Grid Cells
    ************************************/
    MODEL_DONE = false;

    // stop init timer
    timer_stop(&(global_timers[TIMER_VIC_INIT]));
    // start vic run timer
    timer_start(&(global_timers[TIMER_VIC_RUN]));

    while (!MODEL_DONE) {
        read_soilparam(filep.soilparam, &soil_con, &RUN_MODEL, &MODEL_DONE);

        if (RUN_MODEL) {
            cellnum++;

            /** Read Grid Cell Vegetation Parameters **/
            veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
                                    Nveg_type);
            calc_root_fractions(veg_con, &soil_con);

            if (options.LAKES) {
                lake_con =
                    read_lakeparam(filep.lakeparam, soil_con, veg_con);
            }

            /** Build Gridded Filenames, and Open **/
            make_in_and_outfiles(&filep, &filenames, &soil_con,
                                 &streams, dmy);

            /** Reset agg_alarm for Each Stream **/
            for (streamnum = 0;
                 streamnum < (size_t) options.Noutstreams;
                 streamnum++) {
                n = streams[streamnum].agg_alarm.n;
                set_alarm(&(dmy[0]), streams[streamnum].agg_alarm.freq,
                          &n,
                          &(streams[streamnum].agg_alarm));
            }

            /** Read Elevation Band Data if Used **/
            read_snowband(filep.snowband, &soil_con);

            /** Make Top-level Control Structure **/
            all_vars = make_all_vars(veg_con[0].vegetat_type_num);

            /** allocate memory for the veg_hist_struct **/
            alloc_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num,
                           &veg_hist);

            /**************************************************
               Initialize Meteological Forcing Values That
               Have not Been Specifically Set
            **************************************************/

            vic_force(force, dmy, filep.forcing, veg_con, veg_hist, &soil_con);

            /**************************************************
               Initialize Energy Balance and Snow Variables
            **************************************************/

            vic_populate_model_state(&all_vars, filep, soil_con.gridcel,
                                     &soil_con, veg_con, lake_con);

            /** Initialize the storage terms in the water and energy balances **/
            initialize_save_data(&all_vars, &force[0], &soil_con, veg_con,
                                 veg_lib, &lake_con, out_data[0], &save_data,
                                 &cell_timer);

            /******************************************
               Run Model in Grid Cell for all Time Steps
            ******************************************/

            for (rec = startrec; rec < global_param.nrecs; rec++) {
                // Set global reference string (for debugging inside vic_run)
                sprint_dmy(dmy_str, &(dmy[rec]));
                sprintf(vic_run_ref_str,
                        "Gridcell cellnum: %i, timestep info: %s",
                        cellnum, dmy_str);

                /**************************************************
                   Update data structures for current time step
                **************************************************/
                ErrorFlag = update_step_vars(&all_vars, veg_con,
                                             veg_hist[rec]);

                /**************************************************
                   Compute cell physics for 1 timestep
                **************************************************/
                timer_start(&cell_timer);
                ErrorFlag = vic_run(&force[rec], &all_vars,
                                    &(dmy[rec]), &global_param, &lake_con,
                                    &soil_con, veg_con, veg_lib);
                timer_stop(&cell_timer);

                /**************************************************
                   Calculate cell average values for current time step
                **************************************************/
                put_data(&all_vars, &force[rec], &soil_con, veg_con, veg_lib,
                         &lake_con, out_data[0], &save_data, &cell_timer);

                for (streamnum = 0;
                     streamnum < options.Noutstreams;
                     streamnum++) {
                    agg_stream_data(&(streams[streamnum]), &(dmy[rec]),
                                    out_data);
                }

                // Write cell average values for current time step
                write_output(&streams, &dmy[rec]);

                /************************************
                   Save model state at assigned date
                   (after the final time step of the assigned date)
                ************************************/
                if (filep.statefile != NULL &&
                    check_save_state_flag(dmy, rec)) {
                    write_model_state(&all_vars, veg_con->vegetat_type_num,
                                      soil_con.gridcel, &filep, &soil_con);
                }


                if (ErrorFlag == ERROR) {
                    if (options.CONTINUEONERROR) {
                        // Handle grid cell solution error
                        log_warn("ERROR: Grid cell %i failed in record %zu "
                                 "so the simulation has not finished.  An "
                                 "incomplete output file has been "
                                 "generated, check your inputs before "
                                 "rerunning the simulation.",
                                 soil_con.gridcel, rec);
                        break;
                    }
                    else {
                        // Else exit program on cell solution error as in previous versions
                        log_err("ERROR: Grid cell %i failed in record %zu "
                                "so the simulation has ended. Check your "
                                "inputs before rerunning the simulation.",
                                soil_con.gridcel, rec);
                    }
                }
            } /* End Rec Loop */

            close_files(&filep, &streams);

            free_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num,
                          &veg_hist);
            free_all_vars(&all_vars, veg_con[0].vegetat_type_num);
            free_vegcon(&veg_con);
            free((char *) soil_con.AreaFract);
            free((char *) soil_con.BandElev);
            free((char *) soil_con.Tfactor);
            free((char *) soil_con.Pfactor);
            free((char *) soil_con.AboveTreeLine);
        } /* End Run Model Condition */
    }   /* End Grid Loop */

    // stop vic run timer
    timer_stop(&(global_timers[TIMER_VIC_RUN]));
    // start vic final timer
    timer_start(&(global_timers[TIMER_VIC_FINAL]));

    /** cleanup **/
    free_atmos(global_param.nrecs, &force);
    free_dmy(&dmy);
    free_streams(&streams);
    free_out_data(1, out_data);  // 1 is for the number of gridcells, 1 in classic driver
    fclose(filep.soilparam);
    free_veglib(&veg_lib);
    fclose(filep.vegparam);
    fclose(filep.veglib);
    if (options.SNOW_BAND > 1) {
        fclose(filep.snowband);
    }
    if (options.LAKES) {
        fclose(filep.lakeparam);
    }
    if (options.INIT_STATE) {
        fclose(filep.init_state);
    }
    if (options.SAVE_STATE && strcmp(filenames.statefile, "NONE") != 0) {
        fclose(filep.statefile);
    }
    finalize_logging();

    log_info("Completed running VIC %s", VIC_DRIVER);

    // stop vic final timer
    timer_stop(&(global_timers[TIMER_VIC_FINAL]));
    // stop vic all timer
    timer_stop(&(global_timers[TIMER_VIC_ALL]));
    // write timing info
    write_vic_timing_table(global_timers);

    return EXIT_SUCCESS;
}       /* End Main Program */
