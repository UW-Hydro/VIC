/******************************************************************************
 * @section DESCRIPTION
 *
 * Classic driver of the VIC model
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

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
    extern FILE          *LOG_DEST;

    char                  MODEL_DONE;
    char                  RUN_MODEL;
    char                  write_flag;
    size_t                rec;
    size_t                Nveg_type;
    int                   cellnum;
    int                   startrec;
    int                   ErrorFlag;
    dmy_struct           *dmy;
    atmos_data_struct    *atmos;
    veg_hist_struct     **veg_hist;
    veg_con_struct       *veg_con;
    soil_con_struct       soil_con;
    all_vars_struct       all_vars;
    lake_con_struct       lake_con;
    out_data_file_struct *out_data_files;
    out_data_struct      *out_data;
    save_data_struct      save_data;

    /** Read Model Options **/
    cmd_proc(argc, argv, filenames.global);

    // Initialize Log Destination
    initialize_log();

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

    // Set Log Destination
    setup_logging();

    /** Set model constants **/
    if (strcmp(filenames.constants, "MISSING") != 0) {
        filep.constants = open_file(filenames.constants, "r");
        get_parameters(filep.constants);
    }

    /** Set up output data structures **/
    out_data = create_output_list();
    out_data_files = set_output_defaults(out_data);
    fclose(filep.globalparam);
    filep.globalparam = open_file(filenames.global, "r");
    parse_output_info(filep.globalparam, &out_data_files, out_data);

    /** Check and Open Files **/
    check_files(&filep, &filenames);

    if (!options.OUTPUT_FORCE) {
        /** Read Vegetation Library File **/
        veg_lib = read_veglib(filep.veglib, &Nveg_type);
    } /* !OUTPUT_FORCE */

    /** Initialize Parameters **/
    cellnum = -1;

    /** Make Date Data Structure **/
    dmy = make_dmy(&global_param);

    /** allocate memory for the atmos_data_struct **/
    alloc_atmos(global_param.nrecs, &atmos);

    /** Initial state **/
    startrec = 0;
    if (!options.OUTPUT_FORCE) {
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
    } /* !OUTPUT_FORCE */

    /************************************
       Run Model for all Active Grid Cells
    ************************************/
    MODEL_DONE = false;
    while (!MODEL_DONE) {
        soil_con = read_soilparam(filep.soilparam, &RUN_MODEL, &MODEL_DONE);

        if (RUN_MODEL) {
            cellnum++;

            if (!options.OUTPUT_FORCE) {
                /** Read Grid Cell Vegetation Parameters **/
                veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
                                        Nveg_type);
                calc_root_fractions(veg_con, &soil_con);

                if (options.LAKES) {
                    lake_con =
                        read_lakeparam(filep.lakeparam, soil_con, veg_con);
                }
            } /* !OUTPUT_FORCE */

            /** Build Gridded Filenames, and Open **/
            make_in_and_outfiles(&filep, &filenames, &soil_con, out_data_files);

            if (options.PRT_HEADER) {
                /** Write output file headers **/
                write_header(out_data_files, out_data, dmy, global_param);
            }

            if (!options.OUTPUT_FORCE) {
                /** Read Elevation Band Data if Used **/
                read_snowband(filep.snowband, &soil_con);

                /** Make Top-level Control Structure **/
                all_vars = make_all_vars(veg_con[0].vegetat_type_num);

                /** allocate memory for the veg_hist_struct **/
                alloc_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num,
                               &veg_hist);
            } /* !OUTPUT_FORCE */

            /**************************************************
               Initialize Meteological Forcing Values That
               Have not Been Specifically Set
            **************************************************/

            initialize_atmos(atmos, dmy, filep.forcing, veg_lib, veg_con,
                             veg_hist,
                             &soil_con, out_data_files, out_data);

            if (!options.OUTPUT_FORCE) {
                /**************************************************
                   Initialize Energy Balance and Snow Variables
                **************************************************/

                ErrorFlag = initialize_model_state(&all_vars, &global_param,
                                                   filep, soil_con.gridcel,
                                                   veg_con[0].vegetat_type_num,
                                                   options.Nnode,
                                                   atmos[0].air_temp[NR],
                                                   &soil_con, veg_con,
                                                   lake_con);
                if (ErrorFlag == ERROR) {
                    if (options.CONTINUEONERROR) {
                        // Handle grid cell solution error
                        log_warn("ERROR: Grid cell %i failed in record %zu so "
                                 "the simulation has not finished.  An "
                                 "incomplete output file has been generated, "
                                 "check your inputs before rerunning the "
                                 "simulation.\n", soil_con.gridcel, rec);
                        break;
                    }
                    else {
                        // Else exit program on cell solution error as in previous versions
                        log_err("ERROR: Grid cell %i failed in record %zu so "
                                "the simulation has ended. Check your inputs "
                                "before rerunning the simulation.\n",
                                soil_con.gridcel, rec);
                    }
                }

                /** Update Error Handling Structure **/
                Error.filep = filep;
                Error.out_data_files = out_data_files;

                /** Initialize the storage terms in the water and energy balances **/

                /** Sending a negative record number (-global_param.nrecs) to
                    put_data() will accomplish this **/
                ErrorFlag = put_data(&all_vars, &atmos[0], &soil_con, veg_con,
                                     veg_lib, &lake_con, out_data, &save_data,
                                     -global_param.nrecs);

                /******************************************
                   Run Model in Grid Cell for all Time Steps
                ******************************************/

                for (rec = startrec; rec < global_param.nrecs; rec++) {
                    /**************************************************
                       Compute cell physics for 1 timestep
                    **************************************************/
                    ErrorFlag = vic_run(rec, &atmos[rec], &all_vars,
                                        dmy, &global_param, &lake_con,
                                        &soil_con, veg_con, veg_lib,
                                        veg_hist[rec]);

                    /**************************************************
                       Calculate cell average values for current time step
                    **************************************************/
                    write_flag = put_data(&all_vars, &atmos[rec], &soil_con,
                                          veg_con, veg_lib, &lake_con, out_data,
                                          &save_data, rec);

                    // Write cell average values for current time step
                    if (write_flag) {
                        write_output(out_data, out_data_files, &dmy[rec], rec);
                    }

                    /************************************
                       Save model state at assigned date
                       (after the final time step of the assigned date)
                    ************************************/
                    if (filep.statefile != NULL &&
                        (dmy[rec].year == global_param.stateyear &&
                         dmy[rec].month == global_param.statemonth &&
                         dmy[rec].day == global_param.stateday &&
                         (rec + 1 == global_param.nrecs ||
                          dmy[rec + 1].day != global_param.stateday))) {
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
                                     "rerunning the simulation.\n",
                                     soil_con.gridcel, rec);
                            break;
                        }
                        else {
                            // Else exit program on cell solution error as in previous versions
                            log_err("ERROR: Grid cell %i failed in record %zu "
                                    "so the simulation has ended. Check your "
                                    "inputs before rerunning the simulation.\n",
                                    soil_con.gridcel, rec);
                        }
                    }
                } /* End Rec Loop */
            } /* !OUTPUT_FORCE */

            close_files(&filep, out_data_files, &filenames);

            if (!options.OUTPUT_FORCE) {
                free_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num,
                              &veg_hist);
                free_all_vars(&all_vars, veg_con[0].vegetat_type_num);
                free_vegcon(&veg_con);
                free((char *)soil_con.AreaFract);
                free((char *)soil_con.BandElev);
                free((char *)soil_con.Tfactor);
                free((char *)soil_con.Pfactor);
                free((char *)soil_con.AboveTreeLine);
            } /* !OUTPUT_FORCE */
        } /* End Run Model Condition */
    }   /* End Grid Loop */

    /** cleanup **/
    free_atmos(global_param.nrecs, &atmos);
    free_dmy(&dmy);
    free_out_data_files(&out_data_files);
    free_out_data(&out_data);
    fclose(filep.soilparam);
    if (!options.OUTPUT_FORCE) {
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
    } /* !OUTPUT_FORCE */

    return EXIT_SUCCESS;
}       /* End Main Program */
