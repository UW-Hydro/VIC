/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting values for
 * global parameters, model options, and debugging controls.
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
    unsigned int               tmpstartdate;
    unsigned int               tmpenddate;
    unsigned short int         lastday[MONTHS_PER_YEAR];


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
            if (strcasecmp("NODES", optstr) == 0) {
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
                global_param.calendar = calendar_from_chars(flgstr);
            }
            else if (strcasecmp("OUT_TIME_UNITS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                global_param.time_units = timeunits_from_chars(flgstr);
            }
            else if (strcasecmp("FULL_ENERGY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.FULL_ENERGY = true;
                }
                else {
                    options.FULL_ENERGY = false;
                }
            }
            else if (strcasecmp("FROZEN_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.FROZEN_SOIL = true;
                    options.QUICK_FLUX = false;
                }
                else {
                    options.FROZEN_SOIL = false;
                    options.IMPLICIT = false;
                    options.EXP_TRANS = false;
                }
            }
            else if (strcasecmp("QUICK_FLUX", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.QUICK_FLUX = true;
                }
                else {
                    options.QUICK_FLUX = false;
                }
            }
            else if (strcasecmp("QUICK_SOLVE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.QUICK_SOLVE = true;
                }
                else {
                    options.QUICK_SOLVE = false;
                }
            }
            else if ((strcasecmp("NOFLUX",
                                 optstr) == 0) ||
                     (strcasecmp("NO_FLUX", optstr) == 0)) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.NOFLUX = true;
                }
                else {
                    options.NOFLUX = false;
                }
            }
            else if (strcasecmp("IMPLICIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.IMPLICIT = true;
                }
                else {
                    options.IMPLICIT = false;
                }
            }
            else if (strcasecmp("EXP_TRANS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.EXP_TRANS = true;
                }
                else {
                    options.EXP_TRANS = false;
                }
            }
            else if (strcasecmp("SNOW_DENSITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("DENS_SNTHRM", flgstr) == 0) {
                    options.SNOW_DENSITY = DENS_SNTHRM;
                }
                else {
                    options.SNOW_DENSITY = DENS_BRAS;
                }
            }
            else if (strcasecmp("BLOWING", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING = true;
                }
                else {
                    options.BLOWING = false;
                }
            }
            else if (strcasecmp("BLOWING_VAR_THRESHOLD", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING_VAR_THRESHOLD = true;
                }
                else {
                    options.BLOWING_VAR_THRESHOLD = false;
                }
            }
            else if (strcasecmp("BLOWING_CALC_PROB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING_CALC_PROB = true;
                }
                else {
                    options.BLOWING_CALC_PROB = false;
                }
            }
            else if (strcasecmp("BLOWING_SIMPLE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING_SIMPLE = true;
                }
                else {
                    options.BLOWING_SIMPLE = false;
                }
            }
            else if (strcasecmp("BLOWING_FETCH", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING_FETCH = true;
                }
                else {
                    options.BLOWING_FETCH = false;
                }
            }
            else if (strcasecmp("BLOWING_SPATIAL_WIND", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BLOWING_SPATIAL_WIND = true;
                }
                else {
                    options.BLOWING_SPATIAL_WIND = false;
                }
            }
            else if (strcasecmp("CORRPREC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.CORRPREC = true;
                }
                else {
                    options.CORRPREC = false;
                }
            }
            else if (strcasecmp("CLOSE_ENERGY", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.CLOSE_ENERGY = true;
                }
                else {
                    options.CLOSE_ENERGY = false;
                }
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
            }
            else if (strcasecmp("GRND_FLUX_TYPE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("GF_406", flgstr) == 0) {
                    options.GRND_FLUX_TYPE = GF_406;
                }
                else if (strcasecmp("GF_410", flgstr) == 0) {
                    options.GRND_FLUX_TYPE = GF_410;
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
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.SPATIAL_SNOW = true;
                }
                else {
                    options.SPATIAL_SNOW = false;
                }
            }
            else if (strcasecmp("TFALLBACK", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.TFALLBACK = true;
                }
                else {
                    options.TFALLBACK = false;
                }
            }
            else if (strcasecmp("SHARE_LAYER_MOIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.SHARE_LAYER_MOIST = true;
                }
                else {
                    options.SHARE_LAYER_MOIST = false;
                }
            }
            else if (strcasecmp("CANOPY_LAYERS", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &options.Ncanopy);
            }
            else if (strcasecmp("CARBON", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.CARBON = true;
                }
                else {
                    options.CARBON = false;
                }
            }
            else if (strcasecmp("RC_MODE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("RC_PHOTO", flgstr) == 0) {
                    options.RC_MODE = RC_PHOTO;
                }
                else {
                    options.RC_MODE = RC_JARVIS;
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

            /*************************************
               Define forcing files
            *************************************/
            else if (strcasecmp("FORCING1", optstr) == 0) {
                if (strcmp(filenames.f_path_pfx[0], "MISSING") != 0) {
                    log_err("Tried to define FORCING1 twice, if you want to "
                            "use two forcing files, the second must be "
                            "defined as FORCING2");
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
            else if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                set_force_type(cmdstr, file_num, &field);
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
            else if (strcasecmp("DOMAIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.domain);
            }
            else if (strcasecmp("DOMAIN_TYPE", optstr) == 0) {
                get_domain_type(cmdstr);
            }
            else if (strcasecmp("SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.soil);
            }
            else if (strcasecmp("ARC_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    log_err("\"ARC_SOIL\" is no longer a supported option.\n"
                            "Please convert your soil parameter file and "
                            "remove this option from your global file.");
                }
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
                else {
                    options.BASEFLOW = ARNO;
                }
            }
            else if (strcasecmp("JULY_TAVG_SUPPLIED", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.JULY_TAVG_SUPPLIED = false;
                }
                else {
                    options.JULY_TAVG_SUPPLIED = true;
                }
            }
            else if (strcasecmp("ORGANIC_FRACT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.ORGANIC_FRACT = false;
                }
                else {
                    options.ORGANIC_FRACT = true;
                }
            }
            else if (strcasecmp("VEGLIB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veglib);
            }
            else if (strcasecmp("VEGLIB_PHOTO", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGLIB_PHOTO = true;
                }
                else {
                    options.VEGLIB_PHOTO = false;
                }
            }
            else if (strcasecmp("VEGPARAM", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veg);
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
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.LAKE_PROFILE = false;
                }
                else {
                    options.LAKE_PROFILE = true;
                }
            }

            /*************************************
               Define output files
            *************************************/
            else if (strcasecmp("RESULT_DIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.result_dir);
            }
            else if (strcasecmp("OUTPUT_STEPS_PER_DAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %zu", &global_param.output_steps_per_day);
            }
            else if (strcasecmp("SKIPYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.skipyear);
            }
            else if (strcasecmp("ALMA_OUTPUT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.ALMA_OUTPUT = true;
                }
                else {
                    options.ALMA_OUTPUT = false;
                }
            }
            else if (strcasecmp("MOISTFRACT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.MOISTFRACT = true;
                }
                else {
                    options.MOISTFRACT = false;
                }
            }
            else if (strcasecmp("PRT_SNOW_BAND", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.PRT_SNOW_BAND = true;
                }
                else {
                    options.PRT_SNOW_BAND = false;
                }
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
            // vegetation history not yet implemented in image mode
            // TBD: feature in VIC 4.2 that has been ported to classic
            // mode, but that does not exist in image mode (yet)
            else if (strcasecmp("ALBEDO", optstr) == 0 ||
                     strcasecmp("LAI_IN", optstr) == 0 ||
                     strcasecmp("FCANOPY", optstr) == 0) {
                log_err("Time-varying vegetation parameters not implemented "
                        "in image mode");
            }

            /*************************************
               Fail when deprecated options are used.
            *************************************/
            else if (strcasecmp("TIME_STEP", optstr) == 0) {
                log_err("TIME_STEP has been replaced with MODEL_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("SNOW_STEP", optstr) == 0) {
                log_err("SNOW_STEP has been replaced with SNOW_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("OUT_STEP", optstr) == 0) {
                log_err("OUT_STEP has been replaced with OUTPUT_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("FORCE_DT", optstr) == 0) {
                log_err("FORCE_DT has been replaced with FORCE_STEPS_PER_DAY, "
                        "update your global parameter file accordingly");
            }
            else if (strcasecmp("NLAYER", optstr) == 0) {
                log_err("NLAYER has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCE_FORMAT", optstr) == 0) {
                log_err("FORCE_FORMAT has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCE_ENDIAN", optstr) == 0) {
                log_err("FORCE_ENDIAN has been deprecated in the image driver");
            }
            else if (strcasecmp("GRID_DECIMAL", optstr) == 0) {
                log_err("GRID_DECIMAL has been deprecated in the image driver");
            }
            else if (strcasecmp("BINARY_STATE_FILE", optstr) == 0) {
                log_err(
                    "BINARY_STATE_FILE has been deprecated in the image driver");
            }
            else if (strcasecmp("RESOLUTION", optstr) == 0) {
                log_err("RESOLUTION has been deprecated in the image driver");
            }
            else if (strcasecmp("EQUAL_AREA", optstr) == 0) {
                log_err("EQUAL_AREA has been deprecated in the image driver");
            }
            else if (strcasecmp("CONTINUEONERROR", optstr) == 0) {
                log_err(
                    "CONTINUEONERROR has been deprecated in the image driver");
            }
            else if (strcasecmp("BINARY_OUTPUT", optstr) == 0) {
                log_err("BINARY_OUTPUT has been deprecated in the image driver");
            }
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                log_err("COMPRESS has been deprecated in the image driver");
            }
            else if (strcasecmp("PRT_HEADER", optstr) == 0) {
                log_err("PRT_HEADER has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCE_STEPS_PER_DAY", optstr) == 0) {
                log_err(
                    "FORCE_STEPS_PER_DAY has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCEYEAR", optstr) == 0) {
                log_err("FORCEYEAR has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCEMONTH", optstr) == 0) {
                log_err("FORCEMONTH has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCEDAY", optstr) == 0) {
                log_err("FORCEDAY has been deprecated in the image driver");
            }
            else if (strcasecmp("FORCESEC", optstr) == 0) {
                log_err("FORCESEC has been deprecated in the image driver");
            }

            /*************************************
               Fail when classic driver specific options are used
            *************************************/
            else if (strcasecmp("ATMOS_STEPS_PER_DAY", optstr) == 0) {
                log_err("ATMOS_STEPS_PER_DAY is not a valid option for this "
                        "driver.  Update your global parameter file accordingly.");
            }
            else if (strcasecmp("OUTPUT_FORCE", optstr) == 0) {
                log_err("OUTPUT_FORCE is not a valid option for this driver.  "
                        "Update your global parameter file accordingly.");
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
        // For image drivers, set to model timestep
        global_param.atmos_steps_per_day = global_param.model_steps_per_day;
    }
    global_param.atmos_dt = SEC_PER_DAY /
                            (double) global_param.atmos_steps_per_day;

    // Validate the output step
    if (global_param.output_steps_per_day == 0) {
        global_param.output_steps_per_day = global_param.model_steps_per_day;
    }
    if (global_param.output_steps_per_day > global_param.model_steps_per_day) {
        log_err("Invalid value for OUTPUT_STEPS_PER_DAY (%zu).  "
                "OUTPUT_STEPS_PER_DAY must be <= MODEL_STEPS_PER_DAY (%zu)",
                global_param.output_steps_per_day,
                global_param.model_steps_per_day);
    }
    else if (global_param.model_steps_per_day %
             global_param.output_steps_per_day != 0) {
        log_err("Invalid value for OUTPUT_STEPS_PER_DAY (%zu).  "
                "MODEL_STEPS_PER_DAY (%zu) must be a multiple of "
                "OUTPUT_STEPS_PER_DAY.",
                global_param.output_steps_per_day,
                global_param.model_steps_per_day);
    }
    else if (global_param.output_steps_per_day != 1 &&
             global_param.output_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of output steps per day (%zu) > 1 and < "
                "the minimum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines OUTPUT_STEPS_PER_DAY of at "
                "least (%d).", global_param.model_steps_per_day,
                MIN_SUBDAILY_STEPS_PER_DAY,
                MIN_SUBDAILY_STEPS_PER_DAY);
    }
    else if (global_param.output_steps_per_day >
             MAX_SUBDAILY_STEPS_PER_DAY) {
        log_err("The specified number of model steps per day (%zu) > the "
                "the maximum number of subdaily steps per day (%d).  Make "
                "sure that the global file defines MODEL_STEPS_PER_DAY of at "
                "most (%d).", global_param.model_steps_per_day,
                MAX_SUBDAILY_STEPS_PER_DAY,
                MAX_SUBDAILY_STEPS_PER_DAY);
    }
    else {
        global_param.out_dt = SEC_PER_DAY /
                              (double) global_param.output_steps_per_day;
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
        log_err("The specified simulation start second (%u) > 86400  Make sure "
                "that the global file defines a positive integer "
                "for STARTSEC.",
                global_param.startsec);
    }

    // Validate simulation end date and/or number of timesteps
    make_lastday(global_param.endyear, global_param.calendar, lastday);

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
                    "ENDMONTH.", global_param.endmonth);
        }
        if (global_param.endday == 0) {
            log_err("Simulation end day has not been defined.  Make sure "
                    "that the global file defines ENDDAY.");
        }
        else if (global_param.endday > lastday[global_param.endmonth - 1]) {
            log_err("The specified simulation end day (%hu) > the number of "
                    "days in the ENDMONTH (%hu).  Make sure that the global "
                    "file defines a positive integer for ENDDAY.",
                    global_param.endday, global_param.endmonth);
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

    // Get information from the forcing file(s)
    sprintf(filenames.forcing[0], "%s%4d.nc", filenames.f_path_pfx[0],
            global_param.startyear);
    get_forcing_file_info(&param_set, 0);
    if (param_set.N_TYPES[1] != MISSING) {
        sprintf(filenames.forcing[1], "%s%4d.nc", filenames.f_path_pfx[1],
                global_param.startyear);
        get_forcing_file_info(&param_set, 1);
    }

    if (param_set.N_TYPES[1] != MISSING && global_param.forceyear[1] == 0) {
        global_param.forceyear[1] = global_param.forceyear[0];
        global_param.forcemonth[1] = global_param.forcemonth[0];
        global_param.forceday[1] = global_param.forceday[0];
        global_param.forcesec[1] = global_param.forcesec[0];
        global_param.forceskip[1] = 0;
        global_param.forceoffset[1] = global_param.forceskip[1];
    }
    if (param_set.force_steps_per_day[0] == 0) {
        log_err("Forcing file time steps per day has not been "
                "defined.  Make sure that the global file defines "
                "FORCE_STEPS_PER_DAY.");
    }
    else {
        param_set.FORCE_DT[0] = SEC_PER_DAY /
                                (double) param_set.force_steps_per_day[0];
    }
    if (param_set.force_steps_per_day[1] > 0) {
        param_set.FORCE_DT[1] = SEC_PER_DAY /
                                (double) param_set.force_steps_per_day[1];
    }
    else {
        param_set.FORCE_DT[1] = param_set.FORCE_DT[0];
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

    // Validate veg parameter information
    if (strcmp(filenames.veg, "MISSING") == 0) {
        log_err("No vegetation parameter file has been defined.  Make sure "
                "that the global file defines the vegetation parameter "
                "file on the line that begins with \"VEGPARAM\".");
    }
    if (strcmp(filenames.veglib, "MISSING") == 0) {
        log_err("No vegetation library file has been defined.  Make sure "
                "that the global file defines the vegetation library file "
                "on the line that begins with \"VEGLIB\".");
    }
    if (options.ROOT_ZONES == 0) {
        log_err("ROOT_ZONES must be defined to a positive integer greater "
                "than 0, in the global control file.");
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
            log_err("Currently, CARBON==TRUE and VEGLIB_PHOTO==FALSE.  "
                    "If CARBON==TRUE, VEGLIB_PHOTO must be set to TRUE and "
                    "carbon-specific veg parameters must be listed in your "
                    "veg library file.");
        }
    }

    // Validate the elevation band file information
    if (options.SNOW_BAND > 1) {
        if (strcmp(filenames.snowband, "MISSING") == 0) {
            log_err("\"SNOW_BAND\" was specified with %zu elevation bands, "
                    "but no elevation band file has been defined.  "
                    "Make sure that the global file defines the elevation "
                    "band file on the line that begins with \"SNOW_BAND\" "
                    "(after the number of bands).", options.SNOW_BAND);
        }
        if (options.SNOW_BAND > MAX_BANDS) {
            log_err("Global file wants more snow bands (%zu) than are "
                    "defined by MAX_BANDS (%d).  Edit vic_driver_shared.h and "
                    "recompile.", options.SNOW_BAND, MAX_BANDS);
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
            log_err("\"INIT_STATE\" was specified, but no input state file "
                    "has been defined.  Make sure that the global file "
                    "defines the inputstate file on the line that begins "
                    "with \"INIT_STATE\".");
        }
    }

    // Validate the output state file information
    if (options.SAVE_STATE) {
        if (strcmp(filenames.statefile, "MISSING") == 0) {
            log_err("\"SAVE_STATE\" was specified, but no output state "
                    "file has been defined.  Make sure that the global "
                    "file defines the output state file on the line that "
                    "begins with \"SAVE_STATE\".");
        }
        if (global_param.stateyear == 0 || global_param.statemonth == 0 ||
            global_param.stateday == 0) {
            log_err("Incomplete specification of the date to save state "
                    "for state file (%s).\nSpecified date (yyyy-mm-dd): "
                    "%04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and "
                    "STATEDAY are set correctly in your global parameter "
                    "file.", filenames.statefile, global_param.stateyear,
                    global_param.statemonth, global_param.stateday);
        }
        // Check for month, day in range
        make_lastday(global_param.stateyear, global_param.calendar,
                     lastday);
        if (global_param.stateday > lastday[global_param.statemonth - 1] ||
            global_param.statemonth > MONTHS_PER_YEAR ||
            global_param.statemonth < 1 ||
            global_param.stateday < 1) {
            log_err("Unusual specification of the date to save state for "
                    "state file (%s).\nSpecified date (yyyy-mm-dd): "
                    "%04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and "
                    "STATEDAY are set correctly in your global parameter "
                    "file.", filenames.statefile, global_param.stateyear,
                    global_param.statemonth, global_param.stateday);
        }
    }
    // Set the statename here to be able to compare with INIT_STATE name
    if (options.SAVE_STATE) {
        sprintf(filenames.statefile, "%s_%04i%02i%02i", filenames.statefile,
                global_param.stateyear, global_param.statemonth,
                global_param.stateday);
    }
    if (options.INIT_STATE && options.SAVE_STATE &&
        (strcmp(filenames.init_state, filenames.statefile) == 0)) {
        log_err("The save state file (%s) has the same name as the "
                "initialize state file (%s).  The initialize state file "
                "will be destroyed when the save state file is opened.",
                filenames.statefile, filenames.init_state);
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
    if (options.Nnode > MAX_NODES) {
        log_err("Global file wants more soil thermal nodes (%zu) than "
                "are defined by MAX_NODES (%d).  Edit vic_driver_shared.h and "
                "recompile.", options.Nnode, MAX_NODES);
    }
    if (!options.FULL_ENERGY && options.CLOSE_ENERGY) {
        log_err("CLOSE_ENERGY is TRUE but FULL_ENERGY is FALSE. Set "
                "FULL_ENERGY to TRUE to run CLOSE_ENERGY, or set "
                "CLOSE_ENERGY to FALSE.");
    }

    // Validate treeline option
    if (options.COMPUTE_TREELINE && !options.JULY_TAVG_SUPPLIED) {
        log_err("COMPUTE_TREELINE is TRUE but JULY_TAVG_SUPPLIED is "
                "FALSE.\n You must supply July average temperature if"
                "you want to use the treeline option.");
    }

    // Validate lake parameter information
    if (options.LAKES) {
        if (!options.FULL_ENERGY) {
            log_err("FULL_ENERGY must be TRUE if the lake model is to "
                    "be run.");
        }
        if (strcmp(filenames.lakeparam, "MISSING") == 0) {
            log_err("\"LAKES\" was specified, but no lake parameter "
                    "file has been defined.  Make sure that the global "
                    "file defines the lake parameter file on the line that "
                    "begins with \"LAKES\".");
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
