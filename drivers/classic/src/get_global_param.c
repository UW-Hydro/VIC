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

#include <vic_def.h>
#include <vic_run.h>
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
    char                       ErrStr[MAXSTRING];
    int                        file_num;
    int                        field;
    int                        i;
    int                        tmpstartdate;
    int                        tmpenddate;
    int                        lastvalidday;
    int                        lastday[] = {
        31,     /* JANUARY */
        28,     /* FEBRUARY */
        31,     /* MARCH */
        30,     /* APRIL */
        31,     /* MAY */
        30,     /* JUNE */
        31,     /* JULY */
        31,     /* AUGUST */
        30,     /* SEPTEMBER */
        31,     /* OCTOBER */
        30,     /* NOVEMBER */
        31,     /* DECEMBER */
    };

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
            else if (strcasecmp("TIME_STEP", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.dt);
            }
            else if (strcasecmp("SNOW_STEP", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &options.SNOW_STEP);
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
            else if (strcasecmp("STARTHOUR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.starthour);
            }
            else if (strcasecmp("NRECS", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.nrecs);
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
            else if (strcasecmp("CONTINUEONERROR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.CONTINUEONERROR = true;
                }
                else {
                    options.CONTINUEONERROR = false;
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
            else if (strcasecmp("EQUAL_AREA", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.EQUAL_AREA = true;
                }
                else {
                    options.EQUAL_AREA = false;
                }
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
            else if (strcasecmp("LW_TYPE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("LW_TVA", flgstr) == 0) {
                    options.LW_TYPE = LW_TVA;
                }
                else if (strcasecmp("LW_ANDERSON", flgstr) == 0) {
                    options.LW_TYPE = LW_ANDERSON;
                }
                else if (strcasecmp("LW_BRUTSAERT", flgstr) == 0) {
                    options.LW_TYPE = LW_BRUTSAERT;
                }
                else if (strcasecmp("LW_SATTERLUND", flgstr) == 0) {
                    options.LW_TYPE = LW_SATTERLUND;
                }
                else if (strcasecmp("LW_IDSO", flgstr) == 0) {
                    options.LW_TYPE = LW_IDSO;
                }
                else if (strcasecmp("LW_PRATA", flgstr) == 0) {
                    options.LW_TYPE = LW_PRATA;
                }
            }
            else if (strcasecmp("LW_CLOUD", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("LW_CLOUD_DEARDORFF", flgstr) == 0) {
                    options.LW_CLOUD = LW_CLOUD_DEARDORFF;
                }
                else {
                    options.LW_CLOUD = LW_CLOUD_BRAS;
                }
            }
            else if (strcasecmp("MTCLIM_SWE_CORR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.MTCLIM_SWE_CORR = true;
                }
                else {
                    options.MTCLIM_SWE_CORR = false;
                }
            }
            else if (strcasecmp("PLAPSE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.PLAPSE = false;
                }
                else {
                    options.PLAPSE = true;
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
            else if (strcasecmp("VP_INTERP", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VP_INTERP = true;
                }
                else {
                    options.VP_INTERP = false;
                }
            }
            else if (strcasecmp("VP_ITER", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("VP_ITER_NONE", flgstr) == 0) {
                    options.VP_ITER = VP_ITER_NONE;
                }
                if (strcasecmp("VP_ITER_ALWAYS", flgstr) == 0) {
                    options.VP_ITER = VP_ITER_ALWAYS;
                }
                if (strcasecmp("VP_ITER_ANNUAL", flgstr) == 0) {
                    options.VP_ITER = VP_ITER_ANNUAL;
                }
                if (strcasecmp("VP_ITER_CONVERGE", flgstr) == 0) {
                    options.VP_ITER = VP_ITER_CONVERGE;
                }
                else {
                    options.VP_INTERP = VP_ITER_ALWAYS;
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
                sscanf(cmdstr, "%*s %s", filenames.log_path_pfx);
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
            else if (strcasecmp("BINARY_STATE_FILE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FALSE", flgstr) == 0) {
                    options.BINARY_STATE_FILE = false;
                }
                else {
                    options.BINARY_STATE_FILE = true;
                }
            }
            /*************************************
               Define forcing files
            *************************************/
            else if (strcasecmp("FORCING1", optstr) == 0) {
                if (strcmp(filenames.f_path_pfx[0], "MISSING") != 0) {
                    nrerror(
                        "Tried to define FORCING1 twice, if you want to use "
                        "two forcing files, the second must be defined as "
                        "FORCING2");
                }
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[0]);
                file_num = 0;
                field = 0;
            }
            else if (strcasecmp("FORCING2", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.f_path_pfx[1]);
                if (strcasecmp("FALSE", filenames.f_path_pfx[1]) == 0) {
                    strcpy(filenames.f_path_pfx[1], "MISSING");
                }
                file_num = 1;
                field = 0;
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
                    nrerror("FORCE_FORMAT must be either ASCII or BINARY.");
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
                    nrerror("FORCE_ENDIAN must be either BIG or LITTLE.");
                }
            }
            else if (strcasecmp("N_TYPES", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param_set.N_TYPES[file_num]);
            }
            else if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                get_force_type(cmdstr, file_num, &field);
            }
            else if (strcasecmp("FORCE_DT", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu ", &param_set.FORCE_DT[file_num]);
            }
            else if (strcasecmp("FORCEYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forceyear[file_num]);
            }
            else if (strcasecmp("FORCEMONTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forcemonth[file_num]);
            }
            else if (strcasecmp("FORCEDAY", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forceday[file_num]);
            }
            else if (strcasecmp("FORCEHOUR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.forcehour[file_num]);
            }
            else if (strcasecmp("GRID_DECIMAL", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &options.GRID_DECIMAL);
            }
            else if (strcasecmp("WIND_H", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &global_param.wind_h);
            }
            else if (strcasecmp("MEASURE_H", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &global_param.measure_h);
            }
            else if (strcasecmp("ALMA_INPUT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.ALMA_INPUT = true;
                }
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
            else if (strcasecmp("ARC_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    nrerror("\"ARC_SOIL\" is no longer a supported option.  "
                            "Please convert your soil parameter file and "
                            "remove this option from your global file");
                }
            }
            else if (strcasecmp("ARNO_PARAMS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    nrerror(
                        "Please change \"ARNO_PARAMS  TRUE\" to \"BASEFLOW  "
                        "NIJSSEN2001\" in your global parameter file.");
                }
                else {
                    nrerror(
                        "Please change \"ARNO_PARAMS  FALSE\" to \"BASEFLOW  "
                        "ARNO\" in your global parameter file.");
                }
            }
            else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    nrerror(
                        "Please change \"NIJSSEN2001_BASEFLOW  TRUE\" to "
                        "\"BASEFLOW  NIJSSEN2001\" in your global parameter "
                        "file.");
                }
                else {
                    nrerror(
                        "Please change \"NIJSSEN2001_BASEFLOW  FALSE\" to "
                        "\"BASEFLOW  ARNO\" in your global parameter file.");
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
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.JULY_TAVG_SUPPLIED = true;
                }
                else {
                    options.JULY_TAVG_SUPPLIED = false;
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
            else if (strcasecmp("VEGLIB_VEGCOVER", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGLIB_VEGCOVER = true;
                }
                else {
                    options.VEGLIB_VEGCOVER = false;
                }
            }
            else if (strcasecmp("VEGPARAM", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", filenames.veg);
            }
            else if (strcasecmp("GLOBAL_LAI", optstr) == 0) {
                fprintf(stderr,
                        "WARNING: GLOBAL_LAI has been replaced by 2 new "
                        "options: VEGPARAM_LAI (whether the vegparam file "
                        "contains LAI values) and LAI_SRC (where to get LAI "
                        "values).\n");
                fprintf(stderr,
                        "\"GLOBAL_LAI  TRUE\" should now be: \"VEGPARAM_LAI  "
                        "TRUE\" and \"LAI_SRC  FROM_VEGPARAM\".\n");
                fprintf(stderr,
                        "\"GLOBAL_LAI  FALSE\" should now be: \"LAI_SRC  "
                        "FROM_VEGLIB\".\n");
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGPARAM_LAI = true;
                    options.LAI_SRC = FROM_VEGPARAM;
                }
                else {
                    options.LAI_SRC = FROM_VEGLIB;
                }
            }
            else if (strcasecmp("VEGPARAM_LAI", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGPARAM_LAI = true;
                }
                else {
                    options.VEGPARAM_LAI = false;
                }
            }
            else if (strcasecmp("LAI_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGPARAM",
                               flgstr) == 0 ||
                    strcasecmp("LAI_FROM_VEGPARAM", flgstr) == 0) {
                    options.LAI_SRC = FROM_VEGPARAM;
                }
                else {
                    options.LAI_SRC = FROM_VEGLIB;
                }
            }
            else if (strcasecmp("VEGPARAM_VEGCOVER", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGPARAM_VEGCOVER = true;
                }
                else {
                    options.VEGPARAM_VEGCOVER = false;
                }
            }
            else if (strcasecmp("VEGCOVER_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.VEGCOVER_SRC = FROM_VEGPARAM;
                }
                else {
                    options.VEGCOVER_SRC = FROM_VEGLIB;
                }
            }
            else if (strcasecmp("VEGPARAM_ALB", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.VEGPARAM_ALB = true;
                }
                else {
                    options.VEGPARAM_ALB = false;
                }
            }
            else if (strcasecmp("ALB_SRC", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("FROM_VEGPARAM", flgstr) == 0) {
                    options.ALB_SRC = FROM_VEGPARAM;
                }
                else {
                    options.ALB_SRC = FROM_VEGLIB;
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
            else if (strcasecmp("OUT_STEP", optstr) == 0) {
                sscanf(cmdstr, "%*s %u", &global_param.out_dt);
            }
            else if (strcasecmp("SKIPYEAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %hu", &global_param.skipyear);
            }
            else if (strcasecmp("COMPRESS", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.COMPRESS = true;
                }
                else {
                    options.COMPRESS = false;
                }
            }
            else if (strcasecmp("BINARY_OUTPUT", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.BINARY_OUTPUT = true;
                }
                else {
                    options.BINARY_OUTPUT = false;
                }
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
            else if (strcasecmp("OUTPUT_FORCE", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.OUTPUT_FORCE = true;
                }
                else {
                    options.OUTPUT_FORCE = false;
                }
            }
            else if (strcasecmp("PRT_HEADER", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", flgstr);
                if (strcasecmp("TRUE", flgstr) == 0) {
                    options.PRT_HEADER = true;
                }
                else {
                    options.PRT_HEADER = false;
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
            else if (strcasecmp("N_OUTFILES", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUTFILE", optstr) == 0) {
                ; // do nothing
            }
            else if (strcasecmp("OUTVAR", optstr) == 0) {
                ; // do nothing
            }
            /***********************************
               Unrecognized Global Parameter Flag
            ***********************************/
            else {
                fprintf(stderr,
                        "WARNING: Unrecognized option in the global parameter "
                        "file:\n\t%s is unknown - check your spelling\n",
                        optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    /******************************************
       Check for undefined required parameters
    ******************************************/

    // Validate model time step
    if (global_param.dt == HOURS_PER_DAY + 1) {
        nrerror("Model time step has not been defined.  Make sure that the "
                "global file defines TIME_STEP.");
    }
    else if (global_param.dt < 1) {
        sprintf(ErrStr,
                "The specified model time step (%d) < 1 hour.  Make sure that "
                "the global file defines a positive number of hours for "
                "TIME_STEP.", global_param.dt);
        nrerror(ErrStr);
    }

    // Validate the output step
    if (global_param.out_dt == 0) {
        global_param.out_dt = global_param.dt;
    }
    else if (global_param.out_dt < global_param.dt || global_param.out_dt >
             HOURS_PER_DAY ||
             (double)global_param.out_dt / (double)global_param.dt !=
             (double)(global_param.out_dt / global_param.dt)) {
        sprintf(ErrStr,
                "Invalid output step specified. Output step must be an "
                "integer multiple of the model time step; >= model time step "
                "and <= 24");
        nrerror(ErrStr);
    }

    // Validate SNOW_STEP and set NR and NF
    if (global_param.dt < HOURS_PER_DAY && global_param.dt !=
        options.SNOW_STEP) {
        sprintf(ErrStr,
                "If the model step is smaller than daily, the snow model "
                "should run\nat the same time step as the rest of the model.");
        nrerror(ErrStr);
    }
    if (global_param.dt % options.SNOW_STEP != 0 || options.SNOW_STEP >
        global_param.dt) {
        sprintf(ErrStr,
                "SNOW_STEP should be <= TIME_STEP and divide TIME_STEP "
                "evenly.");
        nrerror(ErrStr);
    }
    NF = global_param.dt / options.SNOW_STEP;
    if (NF == 1) {
        NR = 0;
    }
    else {
        NR = NF;
    }

    // Validate simulation start date
    if (global_param.startyear == 0) {
        sprintf(ErrStr,
                "Simulation start year has not been defined.  Make sure that "
                "the global file defines STARTYEAR.");
        nrerror(ErrStr);
    }
    if (global_param.startmonth == 0) {
        sprintf(ErrStr,
                "Simulation start month has not been defined.  Make sure that "
                "the global file defines STARTMONTH.");
        nrerror(ErrStr);
    }
    else if (global_param.startmonth > MONTHS_PER_YEAR) {
        sprintf(ErrStr, "The specified simulation start month (%hu) > 12. Make "
                "sure that the global file defines a positive integer for "
                "STARTMONTH.", global_param.startmonth);
        nrerror(ErrStr);
    }
    if (global_param.startday == 0) {
        sprintf(ErrStr,
                "Simulation start day has not been defined.  Make sure that "
                "the global file defines STARTDAY.");
        nrerror(ErrStr);
    }
    if (global_param.dt == HOURS_PER_DAY) {
        global_param.starthour = 0;
    }
    else if (global_param.starthour == HOURS_PER_DAY + 1) {
        sprintf(ErrStr,
                "Simulation start hour has not been defined, yet model "
                "time step is less than 24 hours.  Make sure that the "
                "global file defines STARTHOUR.");
        nrerror(ErrStr);
    }
    else if (global_param.starthour > HOURS_PER_DAY) {
        sprintf(ErrStr,
                "The specified simulation start hour (%hu) > 24.  Make sure "
                "that the global file defines a positive integer "
                "for STARTHOUR.", global_param.starthour);
        nrerror(ErrStr);
    }


    // Validate simulation end date and/or number of timesteps
    if (global_param.nrecs == 0 && global_param.endyear == 0 &&
        global_param.endmonth == 0 && global_param.endday == 0) {
        sprintf(ErrStr,
                "The model global file MUST define EITHER the number of "
                "records to simulate (NRECS), or the year (ENDYEAR), month "
                "(ENDMONTH), and day (ENDDAY) of the last full simulation day");
        nrerror(ErrStr);
    }
    else if (global_param.nrecs == 0) {
        if (global_param.endyear == 0) {
            sprintf(ErrStr,
                    "Simulation end year has not been defined.  Make sure "
                    "that the global file defines ENDYEAR.");
            nrerror(ErrStr);
        }
        if (global_param.endmonth == 0) {
            sprintf(ErrStr,
                    "Simulation end month has not been defined.  Make sure "
                    "that the global file defines ENDMONTH.");
            nrerror(ErrStr);
        }
        else if (global_param.endmonth > MONTHS_PER_YEAR) {
            sprintf(ErrStr,
                    "The specified simulation end month (%hu) < 0.  Make sure "
                    "that the global file defines a positive integer for "
                    "ENDMONTH.", global_param.endmonth);
            nrerror(ErrStr);
        }
        if (global_param.endday == 0) {
            sprintf(ErrStr,
                    "Simulation end day has not been defined.  Make sure "
                    "that the global file defines ENDDAY.");
            nrerror(ErrStr);
        }
        else if (global_param.endday > lastday[global_param.endmonth]) {
            sprintf(ErrStr,
                    "The specified simulation end day (%hu) > the number of "
                    "days in the ENDMONTH (%hu).  Make sure that the global "
                    "file defines a positive integer for ENDDAY.",
                    global_param.endday, global_param.endmonth);
            nrerror(ErrStr);
        }
        tmpstartdate = global_param.startyear * 10000 +
                       global_param.startmonth * 100 +
                       global_param.startday;
        tmpenddate = global_param.endyear * 10000 + global_param.endmonth *
                     100 +
                     global_param.endday;
        if (tmpenddate < tmpstartdate) {
            sprintf(ErrStr,
                    "The specified simulation end date (%04d-%02d-%02d) is "
                    "EARLIER than the specified start date (%04d-%02d-%02d).",
                    global_param.endyear, global_param.endmonth,
                    global_param.endday,
                    global_param.startyear, global_param.startmonth,
                    global_param.startday);
            nrerror(ErrStr);
        }
    }
    else if (global_param.nrecs < 1) {
        sprintf(ErrStr,
                "The specified duration of simulation (%d) < 1 time step. "
                "Make sure that the global file defines a positive integer "
                "for NRECS.", global_param.nrecs);
        nrerror(ErrStr);
    }

    // Validate forcing files and variables
    if (strcmp(filenames.f_path_pfx[0], "MISSING") == 0) {
        nrerror("No forcing file has been defined.  Make sure that the global "
                "file defines FORCING1.");
    }
    for (i = 0; i < 2; i++) {
        if (i == 0 || (i == 1 && param_set.N_TYPES[i] != 0)) {
            if (param_set.N_TYPES[i] == MISSING) {
                sprintf(ErrStr,
                        "Need to specify the number forcing variables types "
                        "in forcing file %d.",
                        i);
                nrerror(ErrStr);
            }
            if (param_set.FORCE_FORMAT[i] == MISSING) {
                fprintf(stderr, "FORCE_FORMAT%d: %d\n", i,
                        param_set.FORCE_FORMAT[i]);
                sprintf(ErrStr, "Need to specify the FORCE_FORMAT "
                        "(ASCII or BINARY) for forcing file %d.", i);
                nrerror(ErrStr);
            }
            if (param_set.FORCE_INDEX[i][param_set.N_TYPES[i] - 1] == MISSING) {
                sprintf(ErrStr,
                        "Did not define enough forcing variables in forcing "
                        "file %d.",
                        i);
                nrerror(ErrStr);
            }
            if (param_set.FORCE_DT[i] == 0) {
                sprintf(ErrStr,
                        "Must define time steps (FORCE_DT <dt>) in control "
                        "file for focing file %d.",
                        file_num);
                nrerror(ErrStr);
            }
        }
    }
    if (param_set.N_TYPES[1] != MISSING && global_param.forceyear[1] == 0) {
        global_param.forceyear[1] = global_param.forceyear[0];
        global_param.forcemonth[1] = global_param.forcemonth[0];
        global_param.forceday[1] = global_param.forceday[0];
        global_param.forcehour[1] = global_param.forcehour[0];
        global_param.forceskip[1] = 0;
        global_param.forceoffset[1] = global_param.forceskip[1];
    }

    // Validate result directory
    if (strcmp(filenames.result_dir, "MISSING") == 0) {
        nrerror("No results directory has been defined.  Make sure that the "
                "global file defines the result directory on the line that "
                "begins with \"RESULT_DIR\".");
    }

    // Validate soil parameter file information
    if (strcmp(filenames.soil, "MISSING") == 0) {
        nrerror("No soil parameter file has been defined.  Make sure that the "
                "global file defines the soil parameter file on the line that "
                "begins with \"SOIL\".");
    }

    /*******************************************************************************
       Validate parameters required for normal simulations but NOT for OUTPUT_FORCE
    *******************************************************************************/

    if (!options.OUTPUT_FORCE) {
        // Validate veg parameter information
        if (strcmp(filenames.veg, "MISSING") == 0) {
            nrerror("No vegetation parameter file has been defined.  Make "
                    "sure that the global file defines the vegetation "
                    "parameter file on the line that begins with "
                    "\"VEGPARAM\".");
        }
        if (strcmp(filenames.veglib, "MISSING") == 0) {
            nrerror(
                "No vegetation library file has been defined.  Make sure that "
                "the global file defines the vegetation library file on the "
                "line that begins with \"VEGLIB\".");
        }
        if (options.ROOT_ZONES == 0) {
            nrerror("ROOT_ZONES must be defined to a positive integer greater "
                    "than 0, in the global control file.");
        }
        if (options.LAI_SRC == LAI_FROM_VEGPARAM && !options.VEGPARAM_LAI) {
            sprintf(ErrStr,
                    "\"LAI_SRC\" was specified as \"LAI_FROM_VEGPARAM\", but "
                    "\"VEGPARAM_LAI\" was set to \"FALSE\" in the global "
                    "parameter file.  If you want VIC to read LAI values from "
                    "the vegparam file, you MUST make sure the veg param file "
                    "contains 1 line of 12 monthly LAI values for EACH veg "
                    "tile in EACH grid cell, and you MUST specify "
                    "\"VEGPARAM_LAI\" as \"TRUE\" in the global parameter "
                    "file.  Alternatively, if you want VIC to read LAI values "
                    "from the veg library file, set \"LAI_SRC\" ro "
                    "\"LAI_FROM_VEGLIB\" in the global parameter file.  In "
                    "either case, the setting of \"VEGPARAM_LAI\" must be "
                    "consistent with the contents of the veg param file (i.e. "
                    "whether or not it contains LAI values).");
            nrerror(ErrStr);
        }

        // Validate SPATIAL_FROST information
        if (options.SPATIAL_FROST) {
            if (options.Nfrost > MAX_FROST_AREAS) {
                sprintf(ErrStr,
                        "\"SPATIAL_FROST\" was specified with %zu frost "
                        "subareas, which is greater than the maximum of %d.",
                        options.Nfrost, MAX_FROST_AREAS);
                nrerror(ErrStr);
            }
            if (options.Nfrost < 1) {
                sprintf(ErrStr,
                        "\"SPATIAL_FROST\" was specified with %zu frost "
                        "subareas, which is less than the mainmum of 1.",
                        options.Nfrost);
                nrerror(ErrStr);
            }
        }

        // Carbon-cycling options
        if (!options.CARBON) {
            if (options.RC_MODE == RC_PHOTO) {
                fprintf(stderr,
                        "WARNING: If CARBON==FALSE, RC_MODE must be set to "
                        "RC_JARVIS.  Setting RC_MODE to set to RC_JARVIS.\n");
                options.RC_MODE = RC_JARVIS;
            }
        }
        else {
            if (!options.VEGLIB_PHOTO) {
                sprintf(ErrStr,
                        "Currently, CARBON==TRUE and VEGLIB_PHOTO==FALSE.  If "
                        "CARBON==TRUE, VEGLIB_PHOTO must be set to TRUE and "
                        "carbon-specific veg parameters must be listed in "
                        "your veg library file.");
                nrerror(ErrStr);
            }
        }

        // Validate the elevation band file information
        if (options.SNOW_BAND > 1) {
            if (strcmp(filenames.snowband, "MISSING") == 0) {
                sprintf(ErrStr,
                        "\"SNOW_BAND\" was specified with %zu elevation bands, "
                        "but no elevation band file has been defined.  Make "
                        "sure that the global file defines the elevation band "
                        "file on the line that begins with \"SNOW_BAND\" "
                        "(after the number of bands).", options.SNOW_BAND);
                nrerror(ErrStr);
            }
            if (options.SNOW_BAND > MAX_BANDS) {
                sprintf(ErrStr,
                        "Global file wants more snow bands (%zu) than are "
                        "defined by MAX_BANDS (%d).  Edit vic_run/include/vic_def.h "
                        "and recompile.", options.SNOW_BAND,
                        MAX_BANDS);
                nrerror(ErrStr);
            }
        }
        else if (options.SNOW_BAND <= 0) {
            sprintf(ErrStr,
                    "Invalid number of elevation bands specified in global "
                    "file (%zu).  Number of bands must be >= 1.",
                    options.SNOW_BAND);
            nrerror(ErrStr);
        }

        // Validate the input state file information
        if (options.INIT_STATE) {
            if (strcmp(filenames.init_state, "MISSING") == 0) {
                nrerror(
                    "\"INIT_STATE\" was specified, but no input state file "
                    "has been defined.  Make sure that the global file "
                    "defines the inputstate file on the line that begins with "
                    "\"INIT_STATE\".");
            }
        }

        // Validate the output state file information
        if (options.SAVE_STATE) {
            if (strcmp(filenames.statefile, "MISSING") == 0) {
                nrerror(
                    "\"SAVE_STATE\" was specified, but no output state file "
                    "has been defined.  Make sure that the global file "
                    "defines the output state file on the line that begins "
                    "with \"SAVE_STATE\".");
            }
            if (global_param.stateyear == 0 || global_param.statemonth == 0 ||
                global_param.stateday == 0) {
                sprintf(ErrStr,
                        "Incomplete specification of the date to save state "
                        "for state file (%s).\nSpecified date (yyyy-mm-dd): "
                        "%04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and "
                        "STATEDAY are set correctly in your global parameter "
                        "file.\n", filenames.statefile, global_param.stateyear,
                        global_param.statemonth, global_param.stateday);
                nrerror(ErrStr);
            }
            // Check for month, day in range
            lastvalidday = lastday[global_param.statemonth - 1];
            if (global_param.statemonth == 2) {
                if ((global_param.stateyear % 4) == 0 &&
                    ((global_param.stateyear % 100) != 0 ||
                     (global_param.stateyear % 400) ==
                     0)) {
                    lastvalidday = 29;
                }
            }
            if (global_param.stateday > lastvalidday ||
                global_param.statemonth > MONTHS_PER_YEAR ||
                global_param.statemonth < 1 || global_param.stateday > 31 ||
                global_param.stateday <
                1) {
                sprintf(ErrStr,
                        "Unusual specification of the date to save state for "
                        "state file (%s).\nSpecified date (yyyy-mm-dd): "
                        "%04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and "
                        "STATEDAY are set correctly in your global parameter "
                        "file.\n", filenames.statefile, global_param.stateyear,
                        global_param.statemonth, global_param.stateday);
                nrerror(ErrStr);
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
            sprintf(ErrStr,
                    "The save state file (%s) has the same name as the "
                    "initialize state file (%s).  The initialize state file "
                    "will be destroyed when the save state file is opened.",
                    filenames.statefile, filenames.init_state);
            nrerror(ErrStr);
        }

        // Validate soil parameter/simulation mode combinations
        if (options.QUICK_FLUX) {
            if (options.Nnode != 3) {
                fprintf(stderr,
                        "WARNING: To run the model QUICK_FLUX=TRUE, you must "
                        "define exactly 3 soil thermal nodes.  Currently "
                        "Nnodes is set to %zu.  Setting Nnodes to 3.\n",
                        options.Nnode);
                options.Nnode = 3;
            }
            if (options.IMPLICIT || options.EXP_TRANS) {
                sprintf(ErrStr,
                        "To run the model with QUICK_FLUX=TRUE, you cannot "
                        "have IMPLICIT=TRUE or EXP_TRANS=TRUE.");
                nrerror(ErrStr);
            }
        }
        else {
            if (!options.FULL_ENERGY && !options.FROZEN_SOIL) {
                sprintf(ErrStr,
                        "To run the model in water balance mode (both "
                        "FULL_ENERGY and FROZEN_SOIL are FALSE) you MUST set "
                        "QUICK_FLUX to TRUE (or leave QUICK_FLUX out of your "
                        "global parameter file).");
                nrerror(ErrStr);
            }
        }
        if ((options.FULL_ENERGY ||
             options.FROZEN_SOIL) && options.Nlayer < 3) {
            sprintf(ErrStr,
                    "You must define at least 3 soil moisture layers to run "
                    "the model in FULL_ENERGY or FROZEN_SOIL modes.  "
                    "Currently Nlayers is set to %zu.", options.Nlayer);
            nrerror(ErrStr);
        }
        if ((!options.FULL_ENERGY &&
             !options.FROZEN_SOIL) && options.Nlayer < 1) {
            sprintf(ErrStr,
                    "You must define at least 1 soil moisture layer to run "
                    "the model.  Currently Nlayers is set to  %zu.",
                    options.Nlayer);
            nrerror(ErrStr);
        }
        if (options.Nlayer > MAX_LAYERS) {
            sprintf(ErrStr,
                    "Global file wants more soil moisture layers (%zu) than "
                    "are defined by MAX_LAYERS (%d).  Edit vic_run/include/vic_def.h "
                    "and recompile.", options.Nlayer,
                    MAX_LAYERS);
            nrerror(ErrStr);
        }
        if (options.Nnode > MAX_NODES) {
            sprintf(ErrStr,
                    "Global file wants more soil thermal nodes (%zu) than are "
                    "defined by MAX_NODES (%d).  Edit vic_run/include/vic_def.h and "
                    "recompile.", options.Nnode,
                    MAX_NODES);
            nrerror(ErrStr);
        }
        if (!options.FULL_ENERGY && options.CLOSE_ENERGY) {
            sprintf(ErrStr,
                    "CLOSE_ENERGY is TRUE but FULL_ENERGY is FALSE. Set "
                    "FULL_ENERGY to TRUE to run CLOSE_ENERGY, or set "
                    "CLOSE_ENERGY to FALSE.");
            nrerror(ErrStr);
        }
        if (options.COMPUTE_TREELINE && !options.JULY_TAVG_SUPPLIED) {
            sprintf(ErrStr,
                    "COMPUTE_TREELINE is TRUE but JULY_TAVG_SUPPLIED is "
                    "FALSE. \nSet JULY_TAVG_SUPPLIED to TRUE and include July "
                    "Average Temperature in the soil parameter file to run "
                    "COMPUTE_TREELINE, or set both to FALSE.");
            nrerror(ErrStr);
        }
        // Validate lake parameter information
        if (options.LAKES) {
            if (!options.FULL_ENERGY) {
                sprintf(ErrStr,
                        "FULL_ENERGY must be TRUE if the lake model is to be "
                        "run.");
                nrerror(ErrStr);
            }
            if (strcmp(filenames.lakeparam, "MISSING") == 0) {
                nrerror(
                    "\"LAKES\" was specified, but no lake parameter file has "
                    "been defined.  Make sure that the global file defines "
                    "the lake parameter file on the line that begins with "
                    "\"LAKES\".");
            }
            if (global_param.resolution == 0) {
                sprintf(ErrStr,
                        "The model grid cell resolution (RESOLUTION) must be "
                        "defined in the global control file when the lake "
                        "model is active.");
                nrerror(ErrStr);
            }
            if (global_param.resolution > 360 && !options.EQUAL_AREA) {
                sprintf(ErrStr,
                        "For EQUAL_AREA=FALSE, the model grid cell resolution "
                        "(RESOLUTION) must be set to the number of lat or lon "
                        "degrees per grid cell.  This cannot exceed 360.");
                nrerror(ErrStr);
            }
            if (options.COMPUTE_TREELINE) {
                sprintf(ErrStr,
                        "LAKES = TRUE and COMPUTE_TREELINE = TRUE are "
                        "incompatible options.");
                nrerror(ErrStr);
            }
        }

        /*********************************
           Output major options to stderr
        *********************************/
        display_current_settings(DISP_VERSION);
    } // !OUTPUT_FORCE
}
