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

#include <vic_driver_classic.h>

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
    extern size_t              NF, NR;

    char                       cmdstr[MAXSTRING];
    char                       optstr[MAXSTRING];
    char                       flgstr[MAXSTRING];
    char                       flgstr2[MAXSTRING];
    size_t                     file_num;
    int                        field;
    int                        i;
    unsigned int               tmpstartdate;
    unsigned int               tmpenddate;
    unsigned short int         lastday[MONTHS_PER_YEAR];

    file_num = 0;

    /** Read through global control file to find parameters **/

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
            else if (strcasecmp("MODEL_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.model_steps_per_day);
            }
            else if (strcasecmp("SNOW_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.snow_steps_per_day);
            }
            else if (strcasecmp("RUNOFF_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.runoff_steps_per_day);
            }
            else if (strcasecmp("ATMOS_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.atmos_steps_per_day);
            }
            else if (strcasecmp("STARTYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startyear);
            }
            else if (strcasecmp("STARTMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startmonth);
            }
            else if (strcasecmp("STARTDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.startday);
            }
            else if (strcasecmp("STARTSEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.startsec);
            }
            else if (strcasecmp("NRECS", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.nrecs);
            }
            else if (strcasecmp("ENDYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endyear);
            }
            else if (strcasecmp("ENDMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endmonth);
            }
            else if (strcasecmp("ENDDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.endday);
            }
            else if (strcasecmp("CALENDAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.calendar = str_to_calendar(flgstr);
            }
            else if (strcasecmp("OUT_TIME_UNITS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.time_units = str_to_timeunits(flgstr);
            }
            else if (strcasecmp("FULL_ENERGY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.FULL_ENERGY = str_to_bool(flgstr);
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
                options.NOFLUX = str_to_bool(flgstr);
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
            else if (strcasecmp("CONTINUEONERROR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.CONTINUEONERROR = str_to_bool(flgstr);
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
            else if (strcasecmp("EQUAL_AREA", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.EQUAL_AREA = str_to_bool(flgstr);
            }
            else if (strcasecmp("RESOLUTION", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &global_param.resolution);
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
                options.SPATIAL_FROST = str_to_bool(flgstr);
                if (options.SPATIAL_FROST) {
                    options.Nfrost = atoi(flgstr2);
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

            /*************************************
               Define state files
            *************************************/
            else if (strcasecmp("INIT_STATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.INIT_STATE = false;
                }
                else {
                    options.INIT_STATE = true;
                    strcpy(filenames.init_state, flgstr);
                }
            }
            else if (strcasecmp("STATENAME", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.statefile);
                options.SAVE_STATE = true;
            }
            else if (strcasecmp("STATEYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.stateyear);
            }
            else if (strcasecmp("STATEMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.statemonth);
            }
            else if (strcasecmp("STATEDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.stateday);
            }
            else if (strcasecmp("STATESEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.statesec);
            }
            else if (strcasecmp("STATE_FORMAT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("BINARY", flgstr) == 0) {
                    options.STATE_FORMAT = BINARY;
                }
                else if (strcasecmp("ASCII", flgstr) == 0) {
                    options.STATE_FORMAT = ASCII;
                }
                else {
                    log_err("STATE_FORMAT must be either ASCII or BINARY.");
                }
            }

            /*************************************
               Define forcing files
            *************************************/
            else if (strcasecmp("FORCING1", optstr) == 0) {
                if (strcmp(filenames.f_path_pfx[0], "MISSING") != 0) {
                    log_err(
                        "Tried to define FORCING1 twice, if you want to use "
                        "two forcing files, the second must be defined as "
                        "FORCING2");
                }
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[0]);
                file_num = 0;
                field = 0;
                // count the number of forcing variables in this file
                param_set.N_TYPES[file_num] = count_force_vars(gp);
            }
            else if (strcasecmp("FORCING2", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[1]);
                if (strcasecmp("FALSE", filenames.f_path_pfx[1]) == 0) {
                    strcpy(filenames.f_path_pfx[1], "MISSING");
                }
                file_num = 1;
                field = 0;
                // count the number of forcing variables in this file
                param_set.N_TYPES[file_num] = count_force_vars(gp);
            }
            else if (strcasecmp("FORCE_FORMAT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp(flgstr, "BINARY") == 0) {
                    param_set.FORCE_FORMAT[file_num] = BINARY;
                }
                else if (strcasecmp(flgstr, "ASCII") == 0) {
                    param_set.FORCE_FORMAT[file_num] = ASCII;
                }
                else {
                    log_err("FORCE_FORMAT must be either ASCII or BINARY.");
                }
            }
            else if (strcasecmp("FORCE_ENDIAN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp(flgstr, "LITTLE") == 0) {
                    param_set.FORCE_ENDIAN[file_num] = LITTLE;
                }
                else if (strcasecmp(flgstr, "BIG") == 0) {
                    param_set.FORCE_ENDIAN[file_num] = BIG;
                }
                else {
                    log_err("FORCE_ENDIAN must be either BIG or LITTLE.");
                }
            }
            else if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                get_force_type(cmdstr, file_num, &field);
            }
            else if (strcasecmp("FORCE_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu",
                       &param_set.force_steps_per_day[file_num]);
            }
            else if (strcasecmp("FORCEYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu",
                       &global_param.forceyear[file_num]);
            }
            else if (strcasecmp("FORCEMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu",
                       &global_param.forcemonth[file_num]);
            }
            else if (strcasecmp("FORCEDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forceday[file_num]);
            }
            else if (strcasecmp("FORCESEC", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.forcesec[file_num]);
            }
            else if (strcasecmp("GRID_DECIMAL", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &options.GRID_DECIMAL);
            }
            else if (strcasecmp("WIND_H", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &global_param.wind_h);
            }

            /*************************************
               Define parameter files
            *************************************/

            else if (strcasecmp("CONSTANTS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.constants);
            }
            else if (strcasecmp("SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.soil);
            }
            else if (strcasecmp("BASEFLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("NIJSSEN2001", flgstr) == 0) {
                    options.BASEFLOW = NIJSSEN2001;
                }
                else {
                    options.BASEFLOW = ARNO;
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
            else if (strcasecmp("VEGLIB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veglib);
            }
            else if (strcasecmp("VEGLIB_PHOTO", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGLIB_PHOTO = str_to_bool(flgstr);
            }
            else if (strcasecmp("VEGLIB_FCAN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGLIB_FCAN = str_to_bool(flgstr);
            }
            else if (strcasecmp("VEGPARAM", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veg);
            }
            else if (strcasecmp("VEGPARAM_LAI", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGPARAM_LAI = str_to_bool(flgstr);
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
            else if (strcasecmp("VEGPARAM_FCAN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGPARAM_FCAN = str_to_bool(flgstr);
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
            else if (strcasecmp("VEGPARAM_ALB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                options.VEGPARAM_ALB = str_to_bool(flgstr);
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
            else if (strcasecmp("ROOT_ZONES", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.ROOT_ZONES);
            }
            else if (strcasecmp("SNOW_BAND", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu %s", &options.SNOW_BAND,
                       filenames.snowband);
            }
            else if (strcasecmp("LAKES", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.LAKES = false;
                }
                else {
                    options.LAKES = true;
                    strcpy(filenames.lakeparam, flgstr);
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
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                ; // do nothing
            }

            /*************************************
               Fail when deprecated options are used.
            *************************************/
            else if (strcasecmp("ARC_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err("\"ARC_SOIL\" is no longer a supported option.  "
                            "Please convert your soil parameter file and "
                            "remove this option from your global file");
                }
            }
            else if (strcasecmp("ARNO_PARAMS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err(
                        "Please change \"ARNO_PARAMS  TRUE\" to \"BASEFLOW  "
                        "NIJSSEN2001\" in your global parameter file.");
                }
                else {
                    log_err(
                        "Please change \"ARNO_PARAMS  FALSE\" to \"BASEFLOW  "
                        "ARNO\" in your global parameter file.");
                }
            }
            else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err(
                        "Please change \"NIJSSEN2001_BASEFLOW  TRUE\" to "
                        "\"BASEFLOW  NIJSSEN2001\" in your global parameter "
                        "file.");
                }
                else {
                    log_err(
                        "Please change \"NIJSSEN2001_BASEFLOW  FALSE\" to "
                        "\"BASEFLOW  ARNO\" in your global parameter file.");
                }
            }
            else if (strcasecmp("VEGLIB_VEGCOVER", optstr) == 0) {
                log_err("The option VEGLIB_VEGCOVER has been replaced by "
                        "VEGLIB_FCAN.  Please edit your global parameter "
                        "file and re-run.");
            }
            else if (strcasecmp("VEGPARAM_VEGCOVER", optstr) == 0) {
                log_err("The option VEGPARAM_VEGCOVER has been replaced by "
                        "VEGPARAM_FCAN.  Please edit your global parameter "
                        "file and re-run.");
            }
            else if (strcasecmp("VEGCOVER_SRC", optstr) == 0) {
                log_err("The option VEGCOVER_SRC has been replaced by "
                        "FCAN_SRC.  Please edit your global parameter "
                        "file and re-run.");
            }
            else if (strcasecmp("TIME_STEP", optstr) == 0) {
                log_err("TIME_STEP has been replaced with MODEL_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("SNOW_STEP", optstr) == 0) {
                log_err("SNOW_STEP has been replaced with SNOW_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("OUT_DT", optstr) == 0) {
                log_err("OUT_DT has been replaced with OUTPUT_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("FORCE_DT", optstr) == 0) {
                log_err("FORCE_DT has been replaced with FORCE_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("BINARY_OUTPUT", optstr) == 0) {
                log_err("BINARY_OUTPUT has been replaced with OUT_FORMAT, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("BINARY_STATE_FILE", optstr) == 0) {
                log_err("BINARY_STATE_FILE has been replaced with STATE_FORMAT, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("ALMA_OUTPUT", optstr) == 0) {
                log_err("ALMA_OUTPUT has been deprecated, update your global "
                        "parameter file accordingly");
            }
            else if (strcasecmp("MOISTFRACT", optstr) == 0) {
                log_err("MOISTFRACT has been deprecated and has been replaced "
                        "with two new output variables OUT_SOIL_ICE_FRAC and "
                        "OUT_SOIL_LIQ_FRAC, update your global parameter file "
                        "accordingly");
            }
            else if (strcasecmp("PRT_HEADER", optstr) == 0) {
                log_err("PRT_HEADER has been deprecated. All output files "
                        "include a header including pertinent metadata.");
            }
            else if (strcasecmp("PRT_SNOW_BAND", optstr) == 0) {
                log_err("PRT_SNOW_BAND has been deprecated. To output band "
                        "specific variables, directly specify them in the "
                        "global parameter file");
            }
            else if (strcasecmp("SKIPYEAR", optstr) == 0) {
                log_err("SKIPYEAR has been deprecated. To avoid writing output"
                        "to history files, set AGGFREQ == FREQ_NEVER");
            }
            else if (strcasecmp("MAX_SNOW_TEMP", optstr) == 0) {
                log_err("MAX_SNOW_TEMP has been deprecated. To"
                        "specify a maximum snow temperature, use the option"
                        "SNOW_MAX_SNOW_TEMP in the vic constants file.")
            }
            else if (strcasecmp("MIN_RAIN_TEMP", optstr) == 0) {
                log_err("MIN_RAIN_TEMP has been deprecated. To"
                        "specify a minimum rain temperature, use the option"
                        "SNOW_MIN_RAIN_TEMP in the vic constants file.")
            }

            /***********************************
               Unrecognized Global Parameter Flag
            ***********************************/
            else {
                log_warn("Unrecognized option in the global parameter file: "
                         "%s is unknown - check your spelling", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    /******************************************
       Check for undefined required parameters
    ******************************************/

    // Validate model time step
    if (global_param.model_steps_per_day == 0) {
        log_err("Model time steps per day has not been defined.  Make sure "
                "that the global file defines MODEL_STEPS_PER_DAY.");
    }
    else if (global_param.model_steps_per_day != 1 &&
             global_param.model_steps_per_day <
             MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of model steps per day (%zu) > 1 and < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines MODEL_STEPS_PER_DAY of at "
                "least (%d).", global_param.model_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.model_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of model steps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines MODEL_STEPS_PER_DAY of at "
                "most (%d).", global_param.model_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else if ((global_param.model_steps_per_day > HOURS_PER_DAY) &&
             (global_param.model_steps_per_day % HOURS_PER_DAY) != 0) {
        log_err("The specified number of model steps per day (%zu) is > 24 "
                "and is not evenly divided by 24.",
                global_param.model_steps_per_day);
    }
    else {
        global_param.dt = SEC_PER_DAY /
                          (double) global_param.model_steps_per_day;
    }

    // Validate snow model time step
    if (global_param.snow_steps_per_day == 0) {
        log_err("Snow model time steps per day has not been defined.  Make "
                "sure that the global file defines SNOW_STEPS_PER_DAY.");
    }
    else if (global_param.model_steps_per_day != 1 &&
             global_param.snow_steps_per_day !=
             global_param.model_steps_per_day) {
        log_err("If the model step is smaller than daily, the snow model "
                "should run at the same time step as the rest of the model.");
    }
    else if (global_param.snow_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of snow model steps per day (%zu) < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines SNOW_STEPS_PER_DAY of at "
                "least (%d).", global_param.snow_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.snow_steps_per_day > MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of snow steps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines SNOW_STEPS_PER_DAY of at "
                "most (%d).", global_param.snow_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.snow_steps_per_day > HOURS_PER_DAY &&
             global_param.snow_steps_per_day % HOURS_PER_DAY != 0) {
        log_err("The specified number of snow model steps per day (%zu) is > "
                "24 and is not evenly divided by 24.",
                global_param.snow_steps_per_day);
    }
    else if (global_param.snow_steps_per_day %
             global_param.model_steps_per_day != 0) {
        log_err("The specified number of snow model timesteps (%zu) must be "
                "evenly divisible by the number of model timesteps per day "
                "(%zu)", global_param.snow_steps_per_day,
                global_param.model_steps_per_day);
    }
    else {
        global_param.snow_dt = SEC_PER_DAY /
                               (double) global_param.snow_steps_per_day;
    }

    // Validate runoff time step
    if (global_param.runoff_steps_per_day == 0) {
        log_err("Runoff time steps per day has not been defined.  Make "
                "sure that the global file defines RUNOFF_STEPS_PER_DAY.");
    }
    else if (global_param.runoff_steps_per_day <
             MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of runoff steps per day (%zu) < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines RUNOFF_STEPS_PER_DAY of at "
                "least (%d).", global_param.runoff_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.runoff_steps_per_day > HOURS_PER_DAY &&
             global_param.runoff_steps_per_day % HOURS_PER_DAY != 0) {
        log_err("The specified number of runoff steps per day (%zu) is > "
                "24 and is not evenly divided by 24.",
                global_param.runoff_steps_per_day);
    }
    else if (global_param.runoff_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of runoff steps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines RUNOFF_STEPS_PER_DAY of at "
                "most (%d).", global_param.runoff_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.runoff_steps_per_day %
             global_param.model_steps_per_day != 0) {
        log_err("The specified number of runoff timesteps (%zu) must be "
                "evenly divisible by the number of model timesteps per day "
                "(%zu)", global_param.runoff_steps_per_day,
                global_param.model_steps_per_day);
    }
    else {
        global_param.runoff_dt = SEC_PER_DAY /
                                 (double) global_param.runoff_steps_per_day;
    }
    // Validate atmos time step
    if (global_param.atmos_steps_per_day == 0) {
        // For classic driver default to hourly atmos timestep
        global_param.atmos_steps_per_day = HOURS_PER_DAY;
    }
    if (global_param.atmos_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of atmos steps per day (%zu) < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines ATMOS_STEPS_PER_DAY of at "
                "least (%d).", global_param.atmos_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.atmos_steps_per_day > HOURS_PER_DAY &&
             global_param.atmos_steps_per_day % HOURS_PER_DAY != 0) {
        log_err("The specified number of atmos steps per day (%zu) is > "
                "24 and is not evenly divided by 24.",
                global_param.atmos_steps_per_day);
    }
    else if (global_param.atmos_steps_per_day > MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of atmos timesteps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines ATMOS_STEPS_PER_DAY of at "
                "most (%d).", global_param.atmos_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.atmos_steps_per_day %
             global_param.model_steps_per_day != 0) {
        log_err("The specified number of atmos timesteps (%zu) must be "
                "evenly divisible by the number of model timesteps per day "
                "(%zu)", global_param.atmos_steps_per_day,
                global_param.model_steps_per_day);
    }
    else if (global_param.atmos_steps_per_day %
             global_param.snow_steps_per_day != 0) {
        log_err("The specified number of atmos timesteps (%zu) must be evenly "
                "divisible by the number of snow model timesteps per day (%zu)",
                global_param.atmos_steps_per_day,
                global_param.model_steps_per_day);
    }
    else {
        global_param.atmos_dt = SEC_PER_DAY /
                                (double) global_param.atmos_steps_per_day;
    }

    // set NR and NF
    NF = global_param.snow_steps_per_day / global_param.model_steps_per_day;
    if (NF == 1) {
        NR = 0;
    }
    else {
        NR = NF;
    }

    // Validate simulation start date
    if (global_param.startyear == 0) {
        log_err("Simulation start year has not been defined.  Make sure that "
                "the global file defines STARTYEAR.");
    }
    if (global_param.startmonth == 0) {
        log_err("Simulation start month has not been defined.  Make sure that "
                "the global file defines STARTMONTH.");
    }
    else if (global_param.startmonth > MONTHS_PER_YEAR) {
        log_err("The specified simulation start month (%hu) > 12. Make "
                "sure that the global file defines a positive integer for "
                "STARTMONTH.", global_param.startmonth);
    }
    if (global_param.startday == 0) {
        log_err("Simulation start day has not been defined.  Make sure that "
                "the global file defines STARTDAY.");
    }
    if (global_param.model_steps_per_day == 1) {
        global_param.startsec = 0;
    }
    else if (global_param.startsec > SEC_PER_DAY) {
        log_err("The specified simulation start second (%u) > 86400.  Make sure "
                "that the global file defines time between 0 and 86400.",
                global_param.startsec);
    }


    // Validate simulation end date and/or number of timesteps
    make_lastday(global_param.calendar, global_param.endyear, lastday);

    if (global_param.nrecs == 0 && global_param.endyear == 0 &&
        global_param.endmonth == 0 && global_param.endday == 0) {
        log_err("The model global file MUST define EITHER the number of "
                "records to simulate (NRECS), or the year (ENDYEAR), month "
                "(ENDMONTH), and day (ENDDAY) of the last full simulation day");
    }
    else if (global_param.nrecs == 0) {
        if (global_param.endyear == 0) {
            log_err("Simulation end year has not been defined.  Make sure "
                    "that the global file defines ENDYEAR.");
        }
        if (global_param.endmonth == 0) {
            log_err("Simulation end month has not been defined.  Make sure "
                    "that the global file defines ENDMONTH.");
        }
        else if (global_param.endmonth > MONTHS_PER_YEAR) {
            log_err("The specified simulation end month (%hu) < 0.  Make sure "
                    "that the global file defines a positive integer for "
                    "ENDMONTH.",
                    global_param.endmonth);
        }
        if (global_param.endday == 0) {
            log_err("Simulation end day has not been defined.  Make sure "
                    "that the global file defines ENDDAY.");
        }
        else if (global_param.endday > lastday[global_param.endmonth - 1]) {
            log_err("The specified simulation end day (%hu) > the number of "
                    "days in the ENDMONTH (%hu).  Make sure that the global "
                    "file defines a positive integer for ENDDAY.",
                    global_param.endday,
                    global_param.endmonth);
        }
        tmpstartdate = global_param.startyear * 10000 +
                       global_param.startmonth * 100 +
                       global_param.startday;
        tmpenddate = global_param.endyear * 10000 +
                     global_param.endmonth * 100 +
                     global_param.endday;
        if (tmpenddate < tmpstartdate) {
            log_err("The specified simulation end date (%04d-%02d-%02d) is "
                    "EARLIER than the specified start date (%04d-%02d-%02d).",
                    global_param.endyear, global_param.endmonth,
                    global_param.endday,
                    global_param.startyear, global_param.startmonth,
                    global_param.startday);
        }
    }
    else if (global_param.nrecs < 1) {
        log_err("The specified duration of simulation (%zu) < 1 time step. "
                "Make sure that the global file defines a positive integer "
                "for NRECS.", global_param.nrecs);
    }

    // Validate forcing files and variables
    if (strcmp(filenames.f_path_pfx[0], "MISSING") == 0) {
        log_err("No forcing file has been defined.  Make sure that the global "
                "file defines FORCING1.");
    }
    for (i = 0; i < 2; i++) {
        if (i == 0 || (i == 1 && param_set.N_TYPES[i] != 0)) {
            if (param_set.N_TYPES[i] == MISSING) {
                log_err("Need to specify the number forcing variables types "
                        "in forcing file %d.",
                        i);
            }
            if (param_set.FORCE_FORMAT[i] == MISSING) {
                log_err("FORCE_FORMAT%d: %d. Need to specify the FORCE_FORMAT "
                        "(ASCII or BINARY) for forcing file %d.", i,
                        param_set.FORCE_FORMAT[i], i);
            }
            if (param_set.FORCE_INDEX[i][param_set.N_TYPES[i] - 1] == MISSING) {
                log_err("Did not define enough forcing variables in forcing "
                        "file %d.",
                        i);
            }
            if (param_set.force_steps_per_day[i] == 0) {
                log_err("Forcing file %d time steps per day has not been "
                        "defined.  Make sure that the global file defines "
                        "FORCE_STEPS_PER_DAY.", i);
            }
            if (param_set.force_steps_per_day[i] !=
                global_param.snow_steps_per_day) {
                log_err("FORCE_STEPS_PER_DAY must match SNOW_STEPS_PER_DAY");
            }
            else {
                param_set.FORCE_DT[i] = SEC_PER_DAY /
                                        (double) param_set.force_steps_per_day[i
                                        ];
            }
        }
    }
    if (param_set.N_TYPES[1] != MISSING && global_param.forceyear[1] == 0) {
        global_param.forceyear[1] = global_param.forceyear[0];
        global_param.forcemonth[1] = global_param.forcemonth[0];
        global_param.forceday[1] = global_param.forceday[0];
        global_param.forcesec[1] = global_param.forcesec[0];
        global_param.forceskip[1] = 0;
        global_param.forceoffset[1] = global_param.forceskip[1];
    }

    // Validate result directory
    if (strcmp(filenames.result_dir, "MISSING") == 0) {
        log_err("No results directory has been defined.  Make sure that the "
                "global file defines the result directory on the line that "
                "begins with \"RESULT_DIR\".");
    }

    // Validate soil parameter file information
    if (strcmp(filenames.soil, "MISSING") == 0) {
        log_err("No soil parameter file has been defined.  Make sure that the "
                "global file defines the soil parameter file on the line that "
                "begins with \"SOIL\".");
    }

    /*******************************************************************************
       Validate parameters required for normal simulations
    *******************************************************************************/

    // Validate veg parameter information
    if (strcmp(filenames.veg, "MISSING") == 0) {
        log_err("No vegetation parameter file has been defined.  Make "
                "sure that the global file defines the vegetation "
                "parameter file on the line that begins with "
                "\"VEGPARAM\".");
    }
    if (strcmp(filenames.veglib, "MISSING") == 0) {
        log_err(
            "No vegetation library file has been defined.  Make sure that "
            "the global file defines the vegetation library file on the "
            "line that begins with \"VEGLIB\".");
    }
    if (options.ROOT_ZONES == 0) {
        log_err("ROOT_ZONES must be defined to a positive integer greater "
                "than 0, in the global control file.");
    }
    if (options.LAI_SRC == FROM_VEGHIST && !param_set.TYPE[LAI].SUPPLIED) {
        log_err("\"LAI_SRC\" was specified as \"FROM_VEGHIST\", but "
                "\"LAI\" was not specified as an input forcing in the "
                "global parameter file.  If you want VIC to read LAI "
                "values from the veg_hist file, you MUST make sure the veg "
                "hist file contains Nveg columns of LAI values, 1 for "
                "each veg tile in the grid cell, AND specify LAI as a "
                "forcing variable in the veg_hist forcing file in the "
                "global parameter file.");
    }
    if (options.LAI_SRC == FROM_VEGPARAM && !options.VEGPARAM_LAI) {
        log_err("\"LAI_SRC\" was specified as \"FROM_VEGPARAM\", but "
                "\"VEGPARAM_LAI\" was set to \"FALSE\" in the global "
                "parameter file.  If you want VIC to read LAI values from "
                "the vegparam file, you MUST make sure the veg param file "
                "contains 1 line of 12 monthly LAI values for EACH veg "
                "tile in EACH grid cell, and you MUST specify "
                "\"VEGPARAM_LAI\" as \"TRUE\" in the global parameter "
                "file.  Alternatively, if you want VIC to read LAI values "
                "from the veg library file, set \"LAI_SRC\" to "
                "\"FROM_VEGLIB\" in the global parameter file.  In "
                "either case, the setting of \"VEGPARAM_LAI\" must be "
                "consistent with the contents of the veg param file (i.e. "
                "whether or not it contains LAI values).");
    }
    if (options.ALB_SRC == FROM_VEGHIST && !param_set.TYPE[ALBEDO].SUPPLIED) {
        log_err("\"ALB_SRC\" was specified as \"FROM_VEGHIST\", but "
                "\"ALBEDO\" was not specified as an input forcing in the "
                "global parameter file.  If you want VIC to read ALBEDO "
                "values from the veg_hist file, you MUST make sure the veg "
                "hist file contains Nveg columns of ALBEDO values, 1 for "
                "each veg tile in the grid cell, AND specify ALBEDO as a "
                "forcing variable in the veg_hist forcing file in the "
                "global parameter file.");
    }
    if (options.ALB_SRC == FROM_VEGPARAM && !options.VEGPARAM_ALB) {
        log_err("\"ALB_SRC\" was specified as \"FROM_VEGPARAM\", but "
                "\"VEGPARAM_ALB\" was set to \"FALSE\" in the global "
                "parameter file.  If you want VIC to read albedo values from "
                "the vegparam file, you MUST make sure the veg param file "
                "contains 1 line of 12 monthly albedo values for EACH veg "
                "tile in EACH grid cell, and you MUST specify "
                "\"VEGPARAM_ALB\" as \"TRUE\" in the global parameter "
                "file.  Alternatively, if you want VIC to read albedo values "
                "from the veg library file, set \"ALB_SRC\" to "
                "\"FROM_VEGLIB\" in the global parameter file.  In "
                "either case, the setting of \"VEGPARAM_ALB\" must be "
                "consistent with the contents of the veg param file (i.e. "
                "whether or not it contains albedo values).");
    }
    if (options.FCAN_SRC == FROM_VEGHIST &&
        !param_set.TYPE[FCANOPY].SUPPLIED) {
        log_err("\"FCAN_SRC\" was specified as \"FROM_VEGHIST\", but "
                "\"FCANOPY\" was not specified as an input forcing in the "
                "global parameter file.  If you want VIC to read FCANOPY "
                "values from the veg_hist file, you MUST make sure the veg "
                "hist file contains Nveg columns of FCANOPY values, 1 for "
                "each veg tile in the grid cell, AND specify FCANOPY as a "
                "forcing variable in the veg_hist forcing file in the "
                "global parameter file.");
    }
    if (options.FCAN_SRC == FROM_VEGPARAM && !options.VEGPARAM_FCAN) {
        log_err("\"FCAN_SRC\" was specified as \"FROM_VEGPARAM\", but "
                "\"VEGPARAM_FCAN\" was set to \"FALSE\" in the global "
                "parameter file.  If you want VIC to read fcanopy values from "
                "the vegparam file, you MUST make sure the veg param file "
                "contains 1 line of 12 monthly fcanopy values for EACH veg "
                "tile in EACH grid cell, and you MUST specify "
                "\"VEGPARAM_FCAN\" as \"TRUE\" in the global parameter "
                "file.  Alternatively, if you want VIC to read fcanopy values "
                "from the veg library file, set \"FCAN_SRC\" to "
                "\"FROM_VEGLIB\" in the global parameter file.  In "
                "either case, the setting of \"VEGPARAM_FCAN\" must be "
                "consistent with the contents of the veg param file (i.e. "
                "whether or not it contains fcanopy values).");
    }
    if (options.FCAN_SRC == FROM_VEGLIB && !options.VEGLIB_FCAN) {
        log_err("\"FCAN_SRC\" was specified as \"FROM_VEGLIB\", but "
                "\"VEGLIB_FCAN\" was set to \"FALSE\" in the global "
                "parameter file.  If you want VIC to read fcanopy values from "
                "the veglib file, you MUST make sure the veg lib file "
                "contains 1 line of 12 monthly fcanopy values for EACH veg "
                "class, and you MUST specify "
                "\"VEGLIB_FCAN\" as \"TRUE\" in the global parameter "
                "file.  Alternatively, if you want VIC to read fcanopy values "
                "from the veg param file, set \"FCAN_SRC\" to "
                "\"FROM_VEGPARAM\" in the global parameter file.  In "
                "either case, the setting of \"VEGLIB_FCAN\" must be "
                "consistent with the contents of the veg lib file (i.e. "
                "whether or not it contains fcanopy values).");
    }

    // Validate SPATIAL_FROST information
    if (options.SPATIAL_FROST) {
        if (options.Nfrost > MAX_FROST_AREAS) {
            log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                    "subareas, which is greater than the maximum of %d.",
                    options.Nfrost, MAX_FROST_AREAS);
        }
        if (options.Nfrost < 1) {
            log_err("\"SPATIAL_FROST\" was specified with %zu frost "
                    "subareas, which is less than the mainmum of 1.",
                    options.Nfrost);
        }
    }

    // Carbon-cycling options
    if (!options.CARBON) {
        if (options.RC_MODE == RC_PHOTO) {
            log_warn("If CARBON==FALSE, RC_MODE must be set to "
                     "RC_JARVIS.  Setting RC_MODE to set to RC_JARVIS.");
            options.RC_MODE = RC_JARVIS;
        }
    }
    else {
        if (!options.VEGLIB_PHOTO) {
            log_err("Currently, CARBON==TRUE and VEGLIB_PHOTO==FALSE.  If "
                    "CARBON==TRUE, VEGLIB_PHOTO must be set to TRUE and "
                    "carbon-specific veg parameters must be listed in "
                    "your veg library file.");
        }
    }

    // Validate the elevation band file information
    if (options.SNOW_BAND > 1) {
        if (strcmp(filenames.snowband, "MISSING") == 0) {
            log_err("\"SNOW_BAND\" was specified with %zu elevation bands, "
                    "but no elevation band file has been defined.  Make "
                    "sure that the global file defines the elevation band "
                    "file on the line that begins with \"SNOW_BAND\" "
                    "(after the number of bands).", options.SNOW_BAND);
        }
    }
    else if (options.SNOW_BAND <= 0) {
        log_err("Invalid number of elevation bands specified in global "
                "file (%zu).  Number of bands must be >= 1.",
                options.SNOW_BAND);
    }

    // Validate the input state file information
    if (options.INIT_STATE) {
        if (strcmp(filenames.init_state, "MISSING") == 0) {
            log_err(
                "\"INIT_STATE\" was specified, but no input state file "
                "has been defined.  Make sure that the global file "
                "defines the inputstate file on the line that begins with "
                "\"INIT_STATE\".");
        }
    }

    // Validate the output state file information
    if (options.SAVE_STATE) {
        if (strcmp(filenames.statefile, "MISSING") == 0) {
            log_err(
                "\"SAVE_STATE\" was specified, but no output state file "
                "has been defined.  Make sure that the global file "
                "defines the output state file on the line that begins "
                "with \"SAVE_STATE\".");
        }
        if (global_param.stateyear == 0 || global_param.statemonth == 0 ||
            global_param.stateday == 0) {
            log_err("Incomplete specification of the date to save state "
                    "for state file (%s).\nSpecified date (yyyy-mm-dd-sssss): "
                    "%04d-%02d-%02d-%05u\nMake sure STATEYEAR, STATEMONTH, "
                    "and STATEDAY are set correctly in your global parameter "
                    "file.", filenames.statefile, global_param.stateyear,
                    global_param.statemonth, global_param.stateday,
                    global_param.statesec);
        }
        // Check for month, day in range
        make_lastday(global_param.calendar, global_param.stateyear,
                     lastday);
        if (global_param.stateday > lastday[global_param.statemonth - 1] ||
            global_param.statemonth < 1 ||
            global_param.statemonth > MONTHS_PER_YEAR ||
            global_param.stateday < 1 || global_param.stateday > 31 ||
            global_param.statesec > SEC_PER_DAY) {
            log_err("Unusual specification of the date to save state "
                    "for state file (%s).\nSpecified date (yyyy-mm-dd-sssss): "
                    "%04d-%02d-%02d-%05u\nMake sure STATEYEAR, STATEMONTH, "
                    "STATEDAY and STATESEC are set correctly in your global "
                    "parameter file.", filenames.statefile,
                    global_param.stateyear, global_param.statemonth,
                    global_param.stateday, global_param.statesec);
        }
    }
    // Set the statename here to be able to compare with INIT_STATE name
    if (options.SAVE_STATE) {
        sprintf(filenames.statefile, "%s_%04i%02i%02i_%05u",
                filenames.statefile, global_param.stateyear,
                global_param.statemonth, global_param.stateday,
                global_param.statesec);
    }
    if (options.INIT_STATE && options.SAVE_STATE &&
        (strcmp(filenames.init_state, filenames.statefile) == 0)) {
        log_err("The save state file (%s) has the same name as the "
                "initialize state file (%s).  The initialize state file "
                "will be destroyed when the save state file is opened.",
                filenames.statefile, filenames.init_state);
    }

    // Default file formats (if unset)
    if (options.SAVE_STATE && options.STATE_FORMAT == UNSET_FILE_FORMAT) {
        options.STATE_FORMAT = ASCII;
    }

    // Validate soil parameter/simulation mode combinations
    if (options.QUICK_FLUX) {
        if (options.Nnode != 3) {
            log_warn("To run the model QUICK_FLUX=TRUE, you must "
                     "define exactly 3 soil thermal nodes.  Currently "
                     "Nnodes is set to %zu.  Setting Nnodes to 3.",
                     options.Nnode);
            options.Nnode = 3;
        }
        if (options.IMPLICIT || options.EXP_TRANS) {
            log_err("To run the model with QUICK_FLUX=TRUE, you cannot "
                    "have IMPLICIT=TRUE or EXP_TRANS=TRUE.");
        }
    }
    else {
        if (!options.FULL_ENERGY && !options.FROZEN_SOIL) {
            log_err("To run the model in water balance mode (both "
                    "FULL_ENERGY and FROZEN_SOIL are FALSE) you MUST set "
                    "QUICK_FLUX to TRUE (or leave QUICK_FLUX out of your "
                    "global parameter file).");
        }
    }
    if (options.QUICK_SOLVE && !options.QUICK_FLUX) {
        if (options.NOFLUX) {
            log_err("NOFLUX must be set to FALSE when QUICK_SOLVE=TRUE "
                    "and QUICK_FLUX=FALSE");
        }
        if (options.EXP_TRANS) {
            log_err("EXP_TRANS must be set to FALSE when QUICK_SOLVE=TRUE "
                    "and QUICK_FLUX=FALSE");
        }
    }
    if ((options.FULL_ENERGY ||
         options.FROZEN_SOIL) && options.Nlayer < 3) {
        log_err("You must define at least 3 soil moisture layers to run "
                "the model in FULL_ENERGY or FROZEN_SOIL modes.  "
                "Currently Nlayers is set to %zu.", options.Nlayer);
    }
    if ((!options.FULL_ENERGY &&
         !options.FROZEN_SOIL) && options.Nlayer < 1) {
        log_err("You must define at least 1 soil moisture layer to run "
                "the model.  Currently Nlayers is set to  %zu.",
                options.Nlayer);
    }
    if (options.Nlayer > MAX_LAYERS) {
        log_err("Global file wants more soil moisture layers (%zu) than "
                "are defined by MAX_LAYERS (%d).  Edit vic_run/include/vic_def.h "
                "and recompile.", options.Nlayer,
                MAX_LAYERS);
    }
    if (options.Nnode > MAX_NODES) {
        log_err("Global file wants more soil thermal nodes (%zu) than are "
                "defined by MAX_NODES (%d).  Edit vic_run/include/vic_def.h and "
                "recompile.", options.Nnode,
                MAX_NODES);
    }
    if (!options.FULL_ENERGY && options.CLOSE_ENERGY) {
        log_err("CLOSE_ENERGY is TRUE but FULL_ENERGY is FALSE. Set "
                "FULL_ENERGY to TRUE to run CLOSE_ENERGY, or set "
                "CLOSE_ENERGY to FALSE.");
    }
    if (options.COMPUTE_TREELINE && !options.JULY_TAVG_SUPPLIED) {
        log_err("COMPUTE_TREELINE is TRUE but JULY_TAVG_SUPPLIED is "
                "FALSE. Set JULY_TAVG_SUPPLIED to TRUE and include July "
                "Average Temperature in the soil parameter file to run "
                "COMPUTE_TREELINE, or set both to FALSE.");
    }
    // Validate lake parameter information
    if (options.LAKES) {
        if (!options.FULL_ENERGY) {
            log_err("FULL_ENERGY must be TRUE if the lake model is to be "
                    "run.");
        }
        if (strcmp(filenames.lakeparam, "MISSING") == 0) {
            log_err("\"LAKES\" was specified, but no lake parameter file "
                    "has been defined.  Make sure that the global file defines "
                    "the lake parameter file on the line that begins with "
                    "\"LAKES\".");
        }
        if (global_param.resolution == 0) {
            log_err("The model grid cell resolution (RESOLUTION) must be "
                    "defined in the global control file when the lake "
                    "model is active.");
        }
        if (global_param.resolution > 360 && !options.EQUAL_AREA) {
            log_err("For EQUAL_AREA=FALSE, the model grid cell resolution "
                    "(RESOLUTION) must be set to the number of lat or lon "
                    "degrees per grid cell.  This cannot exceed 360.");
        }
        if (options.COMPUTE_TREELINE) {
            log_err("LAKES = TRUE and COMPUTE_TREELINE = TRUE are "
                    "incompatible options.");
        }
    }

    /*********************************
       Output major options
    *********************************/
    display_current_settings(DISP_VERSION);
}
