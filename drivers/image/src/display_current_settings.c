/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine displays the current settings of options defined in the header
 * files and the global parameter file.
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
#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Display the current settings of options defined in the header
 *           files and the global parameter file.
 *****************************************************************************/
void
display_current_settings(int mode)
{
    extern option_struct       options;
    extern param_set_struct    param_set;
    extern global_param_struct global_param;
    extern filenames_struct    filenames;

    int                        file_num;

    if (mode == DISP_VERSION) {
        fprintf(stderr, "***** VIC Version %s *****\n", version);
        return;
    }
    else {
        fprintf(stderr,
                "\n***** VIC Version %s - Current Model Settings *****\n",
                version);
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "COMPILE-TIME OPTIONS (set in .h files)\n");
    fprintf(stderr, "----------------------------------------\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "Output to Screen:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Maximum Array Sizes:\n");
    fprintf(stderr, "MAX_BANDS\t\t%2d\n", MAX_BANDS);
    fprintf(stderr, "MAX_FRONTS\t\t%2d\n", MAX_FRONTS);
    fprintf(stderr, "MAX_FROST_AREAS\t\t\t%2d\n", MAX_FROST_AREAS);
    fprintf(stderr, "MAX_LAKE_NODES\t\t%2d\n", MAX_LAKE_NODES);
    fprintf(stderr, "MAX_LAYERS\t\t%2d\n", MAX_LAYERS);
    fprintf(stderr, "MAX_NODES\t\t%2d\n", MAX_NODES);
    fprintf(stderr, "MAX_VEG\t\t\t%2d\n", MAX_VEG);
    fprintf(stderr, "\n");

    if (mode == DISP_COMPILE_TIME) {
        return;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "RUN-TIME OPTIONS (set in global parameter file)\n");
    fprintf(stderr, "-----------------------------------------------\n");

    fprintf(stderr, "Simulation Dimensions:\n");
    fprintf(stderr, "NLAYER\t\t\t%zu\n", options.Nlayer);
    if (options.EQUAL_AREA) {
        fprintf(stderr, "EQUAL_AREA\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "EQUAL_AREA\t\tFALSE\n");
    }
    fprintf(stderr, "RESOLUTION\t\t%f\n", global_param.resolution);
    fprintf(stderr, "TIME_STEP\t\t%d\n", global_param.dt);
    fprintf(stderr, "SNOW_STEP\t\t%d\n", options.SNOW_STEP);
    fprintf(stderr, "STARTYEAR\t\t%d\n", global_param.startyear);
    fprintf(stderr, "STARTMONTH\t\t%d\n", global_param.startmonth);
    fprintf(stderr, "STARTDAY\t\t%d\n", global_param.startday);
    fprintf(stderr, "STARTHOUR\t\t%d\n", global_param.starthour);
    if (global_param.nrecs > 0) {
        fprintf(stderr, "NRECS\t\t%u\n", global_param.nrecs);
    }
    else {
        fprintf(stderr, "ENDYEAR\t\t\t%d\n", global_param.endyear);
        fprintf(stderr, "ENDMONTH\t\t%d\n", global_param.endmonth);
        fprintf(stderr, "ENDDAY\t\t\t%d\n", global_param.endday);
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Simulation Parameters:\n");
    if (options.AERO_RESIST_CANSNOW == AR_406) {
        fprintf(stderr, "AERO_RESIST_CANSNOW\t\tAR_406\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_406_LS) {
        fprintf(stderr, "AERO_RESIST_CANSNOW\t\tAR_406_LS\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_406_FULL) {
        fprintf(stderr, "AERO_RESIST_CANSNOW\t\tAR_406_FULL\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_410) {
        fprintf(stderr, "AERO_RESIST_CANSNOW\t\tAR_410\n");
    }
    if (options.BLOWING) {
        fprintf(stderr, "BLOWING\t\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "BLOWING\t\t\tFALSE\n");
    }
    if (options.CLOSE_ENERGY) {
        fprintf(stderr, "CLOSE_ENERGY\t\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "CLOSE_ENERGY\t\t\tFALSE\n");
    }
    if (options.COMPUTE_TREELINE) {
        fprintf(stderr, "COMPUTE_TREELINE\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "COMPUTE_TREELINE\t\tFALSE\n");
    }
    if (options.CONTINUEONERROR) {
        fprintf(stderr, "CONTINUEONERROR\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "CONTINUEONERROR\t\tFALSE\n");
    }
    if (options.CORRPREC) {
        fprintf(stderr, "CORRPREC\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "CORRPREC\t\tFALSE\n");
    }
    if (options.EXP_TRANS) {
        fprintf(stderr, "EXP_TRANS\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "EXP_TRANS\t\tFALSE\n");
    }
    if (options.FROZEN_SOIL) {
        fprintf(stderr, "FROZEN_SOIL\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "FROZEN_SOIL\t\tFALSE\n");
    }
    if (options.FULL_ENERGY) {
        fprintf(stderr, "FULL_ENERGY\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "FULL_ENERGY\t\tFALSE\n");
    }
    if (options.GRND_FLUX_TYPE == GF_406) {
        fprintf(stderr, "GRND_FLUX_TYPE\t\tGF_406\n");
    }
    else if (options.GRND_FLUX_TYPE == GF_410) {
        fprintf(stderr, "GRND_FLUX_TYPE\t\tGF_410\n");
    }
    if (options.LOG_MATRIC) {
        fprintf(stderr, "LOG_MATRIC\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "LOG_MATRIC\t\tFALSE\n");
    }
    if (options.LW_TYPE == LW_TVA) {
        fprintf(stderr, "LW_TYPE\t\tLW_TVA\n");
    }
    else if (options.LW_TYPE == LW_ANDERSON) {
        fprintf(stderr, "LW_TYPE\t\tLW_ANDERSON\n");
    }
    else if (options.LW_TYPE == LW_BRUTSAERT) {
        fprintf(stderr, "LW_TYPE\t\tLW_BRUTSAERT\n");
    }
    else if (options.LW_TYPE == LW_SATTERLUND) {
        fprintf(stderr, "LW_TYPE\t\tLW_SATTERLUND\n");
    }
    else if (options.LW_TYPE == LW_IDSO) {
        fprintf(stderr, "LW_TYPE\t\tLW_IDSO\n");
    }
    else if (options.LW_TYPE == LW_PRATA) {
        fprintf(stderr, "LW_TYPE\t\tLW_PRATA\n");
    }
    if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
        fprintf(stderr, "LW_CLOUD\t\tLW_CLOUD_DEARDORFF\n");
    }
    else {
        fprintf(stderr, "LW_CLOUD\t\tLW_CLOUD_BRAS\n");
    }
    if (options.IMPLICIT) {
        fprintf(stderr, "IMPLICIT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "IMPLICIT\t\tFALSE\n");
    }
    if (options.NOFLUX) {
        fprintf(stderr, "NOFLUX\t\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "NOFLUX\t\t\tFALSE\n");
    }
    if (options.MTCLIM_SWE_CORR) {
        fprintf(stderr, "MTCLIM_SWE_CORR\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "MTCLIM_SWE_CORR\t\tFALSE\n");
    }
    if (options.PLAPSE) {
        fprintf(stderr, "PLAPSE\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "PLAPSE\t\tFALSE\n");
    }
    if (options.QUICK_FLUX) {
        fprintf(stderr, "QUICK_FLUX\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "QUICK_FLUX\t\tFALSE\n");
    }
    if (options.QUICK_SOLVE) {
        fprintf(stderr, "QUICK_SOLVE\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "QUICK_SOLVE\t\tFALSE\n");
    }
    if (options.SPATIAL_FROST) {
        fprintf(stderr, "SPATIAL_FROST\t\tTRUE\n");
        fprintf(stderr, "Nfrost\t\t%zu\n", options.Nfrost);
    }
    else {
        fprintf(stderr, "SPATIAL_FROST\t\tFALSE\n");
    }
    if (options.SPATIAL_SNOW) {
        fprintf(stderr, "SPATIAL_SNOW\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "SPATIAL_SNOW\t\tFALSE\n");
    }
    if (options.SNOW_DENSITY == DENS_BRAS) {
        fprintf(stderr, "SNOW_DENSITY\t\tDENS_BRAS\n");
    }
    else if (options.SNOW_DENSITY == DENS_SNTHRM) {
        fprintf(stderr, "SNOW_DENSITY\t\tDENS_SNTHRM\n");
    }
    if (options.TFALLBACK) {
        fprintf(stderr, "TFALLBACK\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "TFALLBACK\t\tFALSE\n");
    }
    if (options.VP_INTERP) {
        fprintf(stderr, "VP_INTERP\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "VP_INTERP\t\tFALSE\n");
    }
    if (options.VP_ITER == VP_ITER_NONE) {
        fprintf(stderr, "VP_ITER\t\tVP_ITER_NONE\n");
    }
    else if (options.VP_ITER == VP_ITER_ALWAYS) {
        fprintf(stderr, "VP_ITER\t\tVP_ITER_ALWAYS\n");
    }
    else if (options.VP_ITER == VP_ITER_ANNUAL) {
        fprintf(stderr, "VP_ITER\t\tVP_ITER_ANNUAL\n");
    }
    else if (options.VP_ITER == VP_ITER_CONVERGE) {
        fprintf(stderr, "VP_ITER\t\tVP_ITER_CONVERGE\n");
    }
    fprintf(stderr, "WIND_H\t\t\t%f\n", global_param.wind_h);
    fprintf(stderr, "MEASURE_H\t\t%f\n", global_param.measure_h);
    fprintf(stderr, "NODES\t\t\t%zu\n", options.Nnode);
    if (options.CARBON) {
        fprintf(stderr, "CARBON\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "CARBON\t\tFALSE\n");
    }
    if (options.SHARE_LAYER_MOIST) {
        fprintf(stderr, "SHARE_LAYER_MOIST\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "SHARE_LAYER_MOIST\t\tFALSE\n");
    }
    fprintf(stderr, "Ncanopy\t\t%zu\n", options.Ncanopy);

    fprintf(stderr, "\n");
    fprintf(stderr, "Input Forcing Data:\n");
    for (file_num = 0; file_num < 2; file_num++) {
        if (global_param.forceyear[file_num] > 0) {
            fprintf(stderr, "Forcing File %d:\t\t%s*\n", file_num + 1,
                    filenames.f_path_pfx[file_num]);
            fprintf(stderr, "FORCEYEAR\t\t%d\n",
                    global_param.forceyear[file_num]);
            fprintf(stderr, "FORCEMONTH\t\t%d\n",
                    global_param.forcemonth[file_num]);
            fprintf(stderr, "FORCEDAY\t\t%d\n",
                    global_param.forceday[file_num]);
            fprintf(stderr, "FORCEHOUR\t\t%d\n",
                    global_param.forcehour[file_num]);
            fprintf(stderr, "N_TYPES\t\t\t%d\n", param_set.N_TYPES[file_num]);
            fprintf(stderr, "FORCE_DT\t\t%d\n", param_set.FORCE_DT[file_num]);
            if (param_set.FORCE_ENDIAN[file_num] == LITTLE) {
                fprintf(stderr, "FORCE_ENDIAN\t\tLITTLE\n");
            }
            else {
                fprintf(stderr, "FORCE_ENDIAN\t\tBIG\n");
            }
            if (param_set.FORCE_FORMAT[file_num] == BINARY) {
                fprintf(stderr, "FORCE_FORMAT\t\tBINARY\n");
            }
            else {
                fprintf(stderr, "FORCE_FORMAT\t\tASCII\n");
            }
        }
    }
    fprintf(stderr, "GRID_DECIMAL\t\t%hu\n", options.GRID_DECIMAL);
    if (options.ALMA_INPUT) {
        fprintf(stderr, "ALMA_INPUT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "ALMA_INPUT\t\tFALSE\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Input Domain Data:\n");
    fprintf(stderr, "Domain file\t\t%s\n", filenames.domain);

    fprintf(stderr, "\n");
    fprintf(stderr, "Constants File\t\t%s\n", filenames.constants);
    fprintf(stderr, "Input Soil Data:\n");
    fprintf(stderr, "Soil file\t\t%s\n", filenames.soil);
    if (options.BASEFLOW == ARNO) {
        fprintf(stderr, "BASEFLOW\t\tARNO\n");
    }
    else if (options.BASEFLOW == NIJSSEN2001) {
        fprintf(stderr, "BASEFLOW\t\tNIJSSEN2001\n");
    }
    if (options.JULY_TAVG_SUPPLIED) {
        fprintf(stderr, "JULY_TAVG_SUPPLIED\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "JULY_TAVG_SUPPLIED\t\tFALSE\n");
    }
    if (options.ORGANIC_FRACT) {
        fprintf(stderr, "ORGANIC_FRACT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "ORGANIC_FRACT\t\tFALSE\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Input Veg Data:\n");
    fprintf(stderr, "Veg library file\t%s\n", filenames.veglib);
    if (options.VEGLIB_PHOTO) {
        fprintf(stderr, "VEGLIB_PHOTO\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "VEGLIB_PHOTO\t\tFALSE\n");
    }
    fprintf(stderr, "Veg param file\t\t%s\n", filenames.veg);
    fprintf(stderr, "ROOT_ZONES\t\t%zu\n", options.ROOT_ZONES);
    if (options.VEGPARAM_LAI) {
        fprintf(stderr, "VEGPARAM_LAI\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "VEGPARAM_LAI\t\tFALSE\n");
    }
    if (options.LAI_SRC == LAI_FROM_VEGPARAM) {
        fprintf(stderr, "LAI_SRC\t\tLAI_FROM_VEGPARAM\n");
    }
    else if (options.LAI_SRC == LAI_FROM_VEGLIB) {
        fprintf(stderr, "LAI_SRC\t\tLAI_FROM_VEGLIB\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Input Elevation Data:\n");
    if (options.SNOW_BAND > 1) {
        fprintf(stderr, "SNOW_BAND\t\t%zu\t%s\n", options.SNOW_BAND,
                filenames.snowband);
    }
    else if (options.SNOW_BAND == 1) {
        fprintf(stderr,
                "SNOW_BAND\t\t%zu\t(no input file needed for SNOW_BAND=1)\n",
                options.SNOW_BAND);
    }
    else {
        fprintf(stderr, "SNOW_BAND\t\t%zu\n", options.SNOW_BAND);
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Input Lake Data:\n");
    if (options.LAKES) {
        fprintf(stderr, "LAKES\t\tTRUE\t%s\n", filenames.lakeparam);
    }
    else {
        fprintf(stderr, "LAKES\t\tFALSE\n");
    }
    if (options.LAKE_PROFILE) {
        fprintf(stderr, "LAKE_PROFILE\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "LAKE_PROFILE\t\tFALSE\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Input State File:\n");
    if (options.INIT_STATE) {
        fprintf(stderr, "INIT_STATE\t\tTRUE\t%s\n", filenames.init_state);
        if (options.BINARY_STATE_FILE) {
            fprintf(stderr, "BINARY_STATE_FILE\tTRUE\n");
        }
        else {
            fprintf(stderr, "BINARY_STATE_FILE\tFALSE\n");
        }
    }
    else {
        fprintf(stderr, "INIT_STATE\t\tFALSE\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Output State File:\n");
    if (options.SAVE_STATE) {
        fprintf(stderr, "SAVE_STATE\t\tTRUE\n");
        fprintf(stderr, "STATENAME\t\t%s\n", filenames.statefile);
        fprintf(stderr, "STATEYEAR\t\t%d\n", global_param.stateyear);
        fprintf(stderr, "STATEMONTH\t\t%d\n", global_param.statemonth);
        fprintf(stderr, "STATEDAY\t\t%d\n", global_param.stateday);
        if (options.BINARY_STATE_FILE) {
            fprintf(stderr, "BINARY_STATE_FILE\tTRUE\n");
        }
        else {
            fprintf(stderr, "BINARY_STATE_FILE\tFALSE\n");
        }
    }
    else {
        fprintf(stderr, "SAVE_STATE\t\tFALSE\n");
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Output Data:\n");
    fprintf(stderr, "Result dir:\t\t%s\n", filenames.result_dir);
    fprintf(stderr, "OUT_STEP\t\t%d\n", global_param.out_dt);
    if (options.ALMA_OUTPUT) {
        fprintf(stderr, "ALMA_OUTPUT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "ALMA_OUTPUT\t\tFALSE\n");
    }
    if (options.BINARY_OUTPUT) {
        fprintf(stderr, "BINARY_OUTPUT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "BINARY_OUTPUT\t\tFALSE\n");
    }
    if (options.COMPRESS) {
        fprintf(stderr, "COMPRESS\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "COMPRESS\t\tFALSE\n");
    }
    if (options.MOISTFRACT) {
        fprintf(stderr, "MOISTFRACT\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "MOISTFRACT\t\tFALSE\n");
    }
    if (options.OUTPUT_FORCE) {
        fprintf(stderr, "OUTPUT_FORCE\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "OUTPUT_FORCE\t\tFALSE\n");
    }
    if (options.PRT_HEADER) {
        fprintf(stderr, "PRT_HEADER\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "PRT_HEADER\t\tFALSE\n");
    }
    if (options.PRT_SNOW_BAND) {
        fprintf(stderr, "PRT_SNOW_BAND\t\tTRUE\n");
    }
    else {
        fprintf(stderr, "PRT_SNOW_BAND\t\tFALSE\n");
    }
    fprintf(stderr, "SKIPYEAR\t\t%d\n", global_param.skipyear);
    fprintf(stderr, "\n");
}
