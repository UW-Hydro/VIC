/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting values for
 * global parameters, model options, and debugging controls.
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

/******************************************************************************
 * @brief    Read the VIC model global control file, getting values for
 *           global parameters, model options, and debugging controls.
 *****************************************************************************/
void
get_global_param(FILE *gp)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern param_set_struct    param_set;
    extern filenames_struct    filenames;

    char                       cmdstr[MAXSTRING];
    char                       optstr[MAXSTRING];
    char                       flgstr[MAXSTRING];
    char                       flgstr2[MAXSTRING];

    /** Read through global control file to find parameters **/

    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);

    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#') {
                fgets(cmdstr, MAXSTRING, gp);
                continue;
            }

            /*************************************
               Get Model Global Parameters
            *************************************/
            if (strcasecmp("NLAYER", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.Nlayer);
            }
            else if (strcasecmp("NODES", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.Nnode);
            }
            else if (strcasecmp("OUT_TIME_UNITS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.time_units = str_to_timeunits(flgstr);
            }
            else if (strcasecmp("FULL_ENERGY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.FULL_ENERGY = str_to_bool(flgstr);
                if (options.FULL_ENERGY == false) {
                    log_warn("FULL_ENERGY is set to FALSE. Please double check "
                             "that this is the setting you intended.");
                }
            }
            else if (strcasecmp("FROZEN_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.FROZEN_SOIL = str_to_bool(flgstr);
                // TODO: move these steps to a option validation
                if (options.FROZEN_SOIL) {
                    options.QUICK_FLUX = false;
                }
                else {
                    options.IMPLICIT = false;
                    options.EXP_TRANS = false;
                }
            }
            else if (strcasecmp("QUICK_FLUX", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.QUICK_FLUX = str_to_bool(flgstr);
            }
            else if (strcasecmp("QUICK_SOLVE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.QUICK_SOLVE = str_to_bool(flgstr);
            }
            else if ((strcasecmp("NOFLUX",
                                 optstr) == 0) ||
                     (strcasecmp("NO_FLUX", optstr) == 0)) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.NOFLUX = str_to_bool(flgstr);
                }
            }
            else if (strcasecmp("IMPLICIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.IMPLICIT = str_to_bool(flgstr);
            }
            else if (strcasecmp("EXP_TRANS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.EXP_TRANS = str_to_bool(flgstr);
            }
            else if (strcasecmp("SNOW_DENSITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("DENS_SNTHRM", flgstr) == 0) {
                    options.SNOW_DENSITY = DENS_SNTHRM;
                }
                else if (strcasecmp("DENS_BRAS", flgstr) == 0) {
                    options.SNOW_DENSITY = DENS_BRAS;
                }
                else {
                    log_err("Unknown SNOW_DENSITY option: %s", flgstr);
                }
            }
            else if (strcasecmp("BLOWING", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING = str_to_bool(flgstr);
            }
            else if (strcasecmp("BLOWING_VAR_THRESHOLD", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING_VAR_THRESHOLD = str_to_bool(flgstr);
            }
            else if (strcasecmp("BLOWING_CALC_PROB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING_CALC_PROB = str_to_bool(flgstr);
            }
            else if (strcasecmp("BLOWING_SIMPLE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING_SIMPLE = str_to_bool(flgstr);
            }
            else if (strcasecmp("BLOWING_FETCH", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING_FETCH = str_to_bool(flgstr);
            }
            else if (strcasecmp("BLOWING_SPATIAL_WIND", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.BLOWING_SPATIAL_WIND = str_to_bool(flgstr);
            }
            else if (strcasecmp("CORRPREC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CORRPREC = str_to_bool(flgstr);
            }
            else if (strcasecmp("CLOSE_ENERGY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CLOSE_ENERGY = str_to_bool(flgstr);
            }
            else if (strcasecmp("COMPUTE_TREELINE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.COMPUTE_TREELINE = false;
                }
                else {
                    options.COMPUTE_TREELINE = true;
                    options.AboveTreelineVeg = atoi(flgstr);
                }
            }
            else if (strcasecmp("AERO_RESIST_CANSNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("AR_406", flgstr) == 0) {
                    options.AERO_RESIST_CANSNOW = AR_406;
                }
                else if (strcasecmp("AR_406_LS", flgstr) == 0) {
                    options.AERO_RESIST_CANSNOW = AR_406_LS;
                }
                else if (strcasecmp("AR_406_FULL", flgstr) == 0) {
                    options.AERO_RESIST_CANSNOW = AR_406_FULL;
                }
                else if (strcasecmp("AR_410", flgstr) == 0) {
                    options.AERO_RESIST_CANSNOW = AR_410;
                }
                else {
                    log_err("Unknown AERO_RESIST_CANSNOW option: %s", flgstr);
                }
            }
            else if (strcasecmp("GRND_FLUX_TYPE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("GF_406", flgstr) == 0) {
                    options.GRND_FLUX_TYPE = GF_406;
                }
                else if (strcasecmp("GF_410", flgstr) == 0) {
                    options.GRND_FLUX_TYPE = GF_410;
                }
                else {
                    log_err("Unknown GRND_FLUX_TYPE option: %s", flgstr);
                }
            }
            else if (strcasecmp("SPATIAL_FROST", optstr) == 0) {
                sscanf(cmdstr, "%*s %s %s", flgstr, flgstr2);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.SPATIAL_FROST = true;
                    options.Nfrost = atoi(flgstr2);
                }
                else {
                    options.SPATIAL_FROST = false;
                }
            }
            else if (strcasecmp("SPATIAL_SNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.SPATIAL_SNOW = str_to_bool(flgstr);
            }
            else if (strcasecmp("TFALLBACK", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.TFALLBACK = str_to_bool(flgstr);
            }
            else if (strcasecmp("SHARE_LAYER_MOIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.SHARE_LAYER_MOIST = str_to_bool(flgstr);
            }
            else if (strcasecmp("CANOPY_LAYERS", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.Ncanopy);
            }
            else if (strcasecmp("CARBON", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CARBON = str_to_bool(flgstr);
            }
            else if (strcasecmp("RC_MODE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("RC_PHOTO", flgstr) == 0) {
                    options.RC_MODE = RC_PHOTO;
                }
                else if (strcasecmp("RC_JARVIS", flgstr) == 0) {
                    options.RC_MODE = RC_JARVIS;
                }
                else {
                    log_err("Unknown RC_MODE option: %s", flgstr);
                }
            }

            /*************************************
               Define log directory
            *************************************/
            else if (strcasecmp("LOG_DIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.log_path);
            }
            else if (strcasecmp("INIT_STATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.INIT_STATE = false;
                }
                else {
                    options.INIT_STATE = true;
                    strcpy(filenames.init_state.nc_filename, flgstr);
                }
            }
            // Define state file format
            else if (strcasecmp("STATE_FORMAT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("NETCDF3_CLASSIC", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF3_CLASSIC;
                }
                else if (strcasecmp("NETCDF3_64BIT_OFFSET", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF3_64BIT_OFFSET;
                }
                else if (strcasecmp("NETCDF4_CLASSIC", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF4_CLASSIC;
                }
                else if (strcasecmp("NETCDF4", flgstr) == 0) {
                    options.STATE_FORMAT = NETCDF4;
                }
                else {
                    log_err("STATE_FORMAT must be either NETCDF3_CLASSIC, "
                            "NETCDF3_64BIT_OFFSET, NETCDF4_CLASSIC, or NETCDF4.");
                }
            }

            /*************************************
               Define parameter files
            *************************************/

            else if (strcasecmp("CONSTANTS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.constants);
            }
            else if (strcasecmp("DOMAIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.domain.nc_filename);
            }
            else if (strcasecmp("DOMAIN_TYPE", optstr) == 0) {
                get_domain_type(cmdstr);
            }
            else if (strcasecmp("PARAMETERS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.params.nc_filename);
            }
            else if (strcasecmp("ARNO_PARAMS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err("Please change \"ARNO_PARAMS TRUE\" to \"BASEFLOW "
                            "NIJSSEN2001\" in your global parameter file.");
                }
                else {
                    log_err("Please change \"ARNO_PARAMS FALSE\" to \"BASEFLOW "
                            "ARNO\" in your global parameter file.");
                }
            }
            else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err("Please change \"NIJSSEN2001_BASEFLOW TRUE\" to "
                            "\"BASEFLOW NIJSSEN2001\" in your global "
                            "parameter file.");
                }
                else {
                    log_err("Please change \"NIJSSEN2001_BASEFLOW FALSE\" to "
                            "\"BASEFLOW ARNO\" in your global parameter file.");
                }
            }
            else if (strcasecmp("BASEFLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("NIJSSEN2001", flgstr) == 0) {
                    options.BASEFLOW = NIJSSEN2001;
                }
                else if (strcasecmp("ARNO", flgstr) == 0) {
                    options.BASEFLOW = ARNO;
                }
                else {
                    log_err("Unknown BASEFLOW option: %s", flgstr);
                }
            }
            else if (strcasecmp("JULY_TAVG_SUPPLIED", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.JULY_TAVG_SUPPLIED = str_to_bool(flgstr);
            }
            else if (strcasecmp("ORGANIC_FRACT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.ORGANIC_FRACT = str_to_bool(flgstr);
            }
            else if (strcasecmp("VEGLIB_PHOTO", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGLIB_PHOTO = str_to_bool(flgstr);
            }
            else if (strcasecmp("LAI_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGHIST", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGHIST;
                }
                else if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGPARAM;
                }
                else if (strcasecmp("FROM_VEGLIB", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGLIB;
                }
                else {
                    log_err("Unrecognized value of LAI_SRC in the global "
                            "control file.");
                }
            }
            else if (strcasecmp("FCAN_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGHIST", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGHIST;
                }
                else if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGPARAM;
                }
                else if (strcasecmp("FROM_VEGLIB", flgstr) == 0) {
                    options.FCAN_SRC = FROM_VEGLIB;
                }
                else if (strcasecmp("FROM_DEFAULT", flgstr) == 0) {
                    options.FCAN_SRC = FROM_DEFAULT;
                }
                else {
                    log_err("Unrecognized value of FCAN_SRC in the global "
                            "control file.");
                }
            }
            else if (strcasecmp("ALB_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGHIST", flgstr) == 0) {
                    options.ALB_SRC = FROM_VEGHIST;
                }
                else if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.ALB_SRC = FROM_VEGPARAM;
                }
                else if (strcasecmp("FROM_VEGLIB", flgstr) == 0) {
                    options.ALB_SRC = FROM_VEGLIB;
                }
                else {
                    log_err("Unrecognized value of ALB_SRC in the global "
                            "control file.");
                }
            }
            else if (strcasecmp("SNOW_BAND", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (str_to_bool(flgstr)) {
                    options.SNOW_BAND = SNOW_BAND_TRUE_BUT_UNSET;
                }
            }
            else if (strcasecmp("LAKES", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.LAKES = false;
                }
                else {
                    options.LAKES = true;
                }
            }
            else if (strcasecmp("LAKE_PROFILE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.LAKE_PROFILE = str_to_bool(flgstr);
            }

            /*************************************
               Define output files
            *************************************/
            else if (strcasecmp("RESULT_DIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.result_dir);
            }

            /*************************************
               Define output file contents
            *************************************/
            else if (strcasecmp("OUTFILE", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUTVAR", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("AGGFREQ", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUTPUT_STEPS_PER_DAY", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("SKIPYEAR", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                ; // do nothing
            }
            // vegetation history not yet implemented in image mode
            // TBD: feature in VIC 4.2 that has been ported to classic
            // mode, but that does not exist in image mode (yet)
            else if (strcasecmp("ALBEDO", optstr) == 0 ||
                     strcasecmp("LAI", optstr) == 0 ||
                     strcasecmp("FCANOPY", optstr) == 0) {
                log_err("Time-varying vegetation parameters not implemented "
                        "in CESM driver");
            }

            /***********************************
               Unrecognized Global Parameter Flag
            ***********************************/
            else {
                log_warn("Unrecognized option in the global parameter file: %s"
                         "\n - check your spelling", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    // Output major options
    display_current_settings(DISP_VERSION);
}

/******************************************************************************
 * @brief    Validate the filenames_struct
 *****************************************************************************/
void
validate_filenames(filenames_struct *filenames)
{
    extern option_struct options;

    // Validate log directory
    if (strcmp(filenames->log_path, "MISSING") == 0) {
        log_err("Log directory must be specified in CESM driver");
    }

    // Validate result directory
    if (strcmp(filenames->result_dir, "MISSING") == 0) {
        log_err("No results directory has been defined.  Make sure that the "
                "global file defines the result directory on the line that "
                "begins with \"RESULT_DIR\".");
    }

    // Validate parameter file information
    if (strcmp(filenames->params.nc_filename, "MISSING") == 0) {
        log_err("A parameters file has not been defined.  Make sure that the "
                "global file defines the parameters parameter file on the line "
                "that begins with \"PARAMETERS\".");
    }

    // Validate snow parameter file
    if (options.SNOW_BAND > 1) {
        log_err("Snowbands not implemented in CESM driver");
    }

    // Validate lake parameter information
    if (options.LAKES) {
        log_err("Lakes are not implemented in CESM driver");
    }
}

/******************************************************************************
 * @brief    Validate the global_param_struct
 *****************************************************************************/
void
validate_global_param(global_param_struct *gp)
{
    if (gp->model_steps_per_day != gp->snow_steps_per_day) {
        log_err("snow_steps_per_day must match model_steps_per_day");
    }
    if (gp->model_steps_per_day != gp->runoff_steps_per_day) {
        log_err("runoff_steps_per_day must match model_steps_per_day");
    }
    if (gp->model_steps_per_day != gp->atmos_steps_per_day) {
        log_err("atmos_steps_per_day must match model_steps_per_day");
    }
}

/******************************************************************************
 * @brief    Validate the option_struct
 *****************************************************************************/
void
validate_options(option_struct *options)
{
    // Validate SPATIAL_FROST information
    if (options->SPATIAL_FROST) {
        if (options->Nfrost > MAX_FROST_AREAS) {
            log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                    "subareas, which is greater than the maximum of %d.",
                    options->Nfrost, MAX_FROST_AREAS);
        }
        if (options->Nfrost < 1) {
            log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                    "subareas, which is less than the mainmum of 1.",
                    options->Nfrost);
        }
    }

    // Carbon-cycling options
    if (!options->CARBON) {
        if (options->RC_MODE == RC_PHOTO) {
            log_warn("If CARBON==FALSE, RC_MODE must be set to "
                     "RC_JARVIS.  Setting RC_MODE to set to RC_JARVIS.");
            options->RC_MODE = RC_JARVIS;
        }
    }
    else {
        if (!options->VEGLIB_PHOTO) {
            log_err("Currently, CARBON==TRUE and VEGLIB_PHOTO==FALSE.  "
                    "If CARBON==TRUE, VEGLIB_PHOTO must be set to TRUE and "
                    "carbon-specific veg parameters must be listed in your "
                    "veg library file.");
        }
    }

    // Validate the elevation band file information
    if (options->SNOW_BAND <= 0) {
        log_err("Invalid number of elevation bands specified in global "
                "file (%zu).  Number of bands must be >= 1.",
                options->SNOW_BAND);
    }

    // Validate soil parameter/simulation mode combinations
    if (options->QUICK_FLUX) {
        if (options->Nnode != 3) {
            log_warn("To run the model QUICK_FLUX=TRUE, you must "
                     "define exactly 3 soil thermal nodes.  Currently "
                     "Nnodes is set to %zu.  Setting Nnodes to 3.",
                     options->Nnode);
            options->Nnode = 3;
        }
        if (options->IMPLICIT || options->EXP_TRANS) {
            log_err("To run the model with QUICK_FLUX=TRUE, you cannot "
                    "have IMPLICIT=TRUE or EXP_TRANS=TRUE.");
        }
    }
    if (options->QUICK_SOLVE && !options->QUICK_FLUX) {
        if (options->NOFLUX) {
            log_err("NOFLUX must be set to FALSE when QUICK_SOLVE=TRUE "
                    "and QUICK_FLUX=FALSE");
        }
        if (options->EXP_TRANS) {
            log_err("EXP_TRANS must be set to FALSE when QUICK_SOLVE=TRUE "
                    "and QUICK_FLUX=FALSE");
        }
    }
    if (options->Nlayer < 3) {
        log_err("You must define at least 3 soil moisture layers to run "
                "the model in FULL_ENERGY or FROZEN_SOIL modes.  "
                "Currently Nlayers is set to  %zu.", options->Nlayer);
    }
    if (options->Nlayer > MAX_LAYERS) {
        log_err("Global file wants more soil moisture layers (%zu) than "
                "are defined by MAX_LAYERS (%d).  Edit vicNl_def.h and "
                "recompile.", options->Nlayer, MAX_LAYERS);
    }
    if (options->Nnode > MAX_NODES) {
        log_err("Global file wants more soil thermal nodes (%zu) than "
                "are defined by MAX_NODES (%d).  Edit vicNl_def.h and "
                "recompile.", options->Nnode, MAX_NODES);
    }

    // Validate treeline option
    if (options->COMPUTE_TREELINE) {
        log_err("COMPUTE_TREELINE not implemented in cesm driver");
    }
    if (options->COMPUTE_TREELINE && !options->JULY_TAVG_SUPPLIED) {
        log_err("COMPUTE_TREELINE is TRUE but JULY_TAVG_SUPPLIED is "
                "FALSE.\n You must supply July average temperature if"
                "you want to use the treeline option.");
    }

    // Validate lake parameter information
    if (options->LAKES) {
        if (options->COMPUTE_TREELINE) {
            log_err("LAKES = TRUE and COMPUTE_TREELINE = TRUE are "
                    "incompatible options.");
        }
    }

    // Default file formats (if unset)
    if (options->STATE_FORMAT == UNSET_FILE_FORMAT) {
        options->STATE_FORMAT = NETCDF3_CLASSIC;
    }
}
