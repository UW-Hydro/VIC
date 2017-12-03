/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine displays the current settings of options defined in the header
 * files and the global parameter file.
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

#include <vic_driver_python.h>

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

    print_version(VIC_DRIVER);

    if (mode == DISP_VERSION) {
        return;
    }

    fprintf(LOG_DEST, "\nCurrent Model Settings\n");

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "COMPILE-TIME OPTIONS (set in .h files)\n");
    fprintf(LOG_DEST, "----------------------------------------\n\n");

    fprintf(LOG_DEST, "VIC_DRIVER:\t\t%s\n", VIC_DRIVER);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Maximum Array Sizes:\n");
    fprintf(LOG_DEST, "MAX_FRONTS\t\t%2d\n", MAX_FRONTS);
    fprintf(LOG_DEST, "MAX_FROST_AREAS\t\t\t%2d\n", MAX_FROST_AREAS);
    fprintf(LOG_DEST, "MAX_LAKE_NODES\t\t%2d\n", MAX_LAKE_NODES);
    fprintf(LOG_DEST, "MAX_ZWTVMOIST\t\t%2d\n", MAX_ZWTVMOIST);
    fprintf(LOG_DEST, "MAX_LAYERS\t\t%2d\n", MAX_LAYERS);
    fprintf(LOG_DEST, "MAX_NODES\t\t%2d\n", MAX_NODES);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "MINSOILDEPTH\t\t%f\n", MINSOILDEPTH);
    fprintf(LOG_DEST, "MIN_FCANOPY\t\t%f\n", MIN_FCANOPY);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "MIN_SUBDAILY_STEPS_PER_DAY %d\n",
            MIN_SUBDAILY_STEPS_PER_DAY);
    fprintf(LOG_DEST, "MAX_SUBDAILY_STEPS_PER_DAY %d\n",
            MAX_SUBDAILY_STEPS_PER_DAY);
    fprintf(LOG_DEST, "\n");

    if (mode == DISP_COMPILE_TIME) {
        return;
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "RUN-TIME OPTIONS (set in global parameter file)\n");
    fprintf(LOG_DEST, "-----------------------------------------------\n");

    fprintf(LOG_DEST, "Simulation Dimensions:\n");
    fprintf(LOG_DEST, "NLAYER\t\t\t%zu\n", options.Nlayer);
    if (options.EQUAL_AREA) {
        fprintf(LOG_DEST, "EQUAL_AREA\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "EQUAL_AREA\t\tFALSE\n");
    }
    fprintf(LOG_DEST, "RESOLUTION\t\t%f\n", global_param.resolution);
    fprintf(LOG_DEST, "MODEL_DT\t\t%f\n", global_param.dt);
    fprintf(LOG_DEST, "SNOW_DT\t\t%f\n", global_param.snow_dt);
    fprintf(LOG_DEST, "RUNOFF_DT\t\t%f\n", global_param.runoff_dt);
    fprintf(LOG_DEST, "STARTYEAR\t\t%d\n", global_param.startyear);
    fprintf(LOG_DEST, "STARTMONTH\t\t%d\n", global_param.startmonth);
    fprintf(LOG_DEST, "STARTDAY\t\t%d\n", global_param.startday);
    fprintf(LOG_DEST, "STARTSEC\t\t%d\n", global_param.startsec);
    if (global_param.nrecs > 0) {
        fprintf(LOG_DEST, "NRECS\t\t%zu\n", global_param.nrecs);
    }
    else {
        fprintf(LOG_DEST, "ENDYEAR\t\t\t%d\n", global_param.endyear);
        fprintf(LOG_DEST, "ENDMONTH\t\t%d\n", global_param.endmonth);
        fprintf(LOG_DEST, "ENDDAY\t\t\t%d\n", global_param.endday);
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Simulation Parameters:\n");
    if (options.AERO_RESIST_CANSNOW == AR_406) {
        fprintf(LOG_DEST, "AERO_RESIST_CANSNOW\t\tAR_406\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_406_LS) {
        fprintf(LOG_DEST, "AERO_RESIST_CANSNOW\t\tAR_406_LS\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_406_FULL) {
        fprintf(LOG_DEST, "AERO_RESIST_CANSNOW\t\tAR_406_FULL\n");
    }
    else if (options.AERO_RESIST_CANSNOW == AR_410) {
        fprintf(LOG_DEST, "AERO_RESIST_CANSNOW\t\tAR_410\n");
    }
    if (options.BLOWING) {
        fprintf(LOG_DEST, "BLOWING\t\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "BLOWING\t\t\tFALSE\n");
    }
    if (options.CLOSE_ENERGY) {
        fprintf(LOG_DEST, "CLOSE_ENERGY\t\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "CLOSE_ENERGY\t\t\tFALSE\n");
    }
    if (options.COMPUTE_TREELINE) {
        fprintf(LOG_DEST, "COMPUTE_TREELINE\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "COMPUTE_TREELINE\t\tFALSE\n");
    }
    if (options.CONTINUEONERROR) {
        fprintf(LOG_DEST, "CONTINUEONERROR\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "CONTINUEONERROR\t\tFALSE\n");
    }
    if (options.CORRPREC) {
        fprintf(LOG_DEST, "CORRPREC\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "CORRPREC\t\tFALSE\n");
    }
    if (options.EXP_TRANS) {
        fprintf(LOG_DEST, "EXP_TRANS\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "EXP_TRANS\t\tFALSE\n");
    }
    if (options.FROZEN_SOIL) {
        fprintf(LOG_DEST, "FROZEN_SOIL\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "FROZEN_SOIL\t\tFALSE\n");
    }
    if (options.FULL_ENERGY) {
        fprintf(LOG_DEST, "FULL_ENERGY\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "FULL_ENERGY\t\tFALSE\n");
    }
    if (options.GRND_FLUX_TYPE == GF_406) {
        fprintf(LOG_DEST, "GRND_FLUX_TYPE\t\tGF_406\n");
    }
    else if (options.GRND_FLUX_TYPE == GF_410) {
        fprintf(LOG_DEST, "GRND_FLUX_TYPE\t\tGF_410\n");
    }
    if (options.IMPLICIT) {
        fprintf(LOG_DEST, "IMPLICIT\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "IMPLICIT\t\tFALSE\n");
    }
    if (options.NOFLUX) {
        fprintf(LOG_DEST, "NOFLUX\t\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "NOFLUX\t\t\tFALSE\n");
    }
    if (options.QUICK_FLUX) {
        fprintf(LOG_DEST, "QUICK_FLUX\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "QUICK_FLUX\t\tFALSE\n");
    }
    if (options.QUICK_SOLVE) {
        fprintf(LOG_DEST, "QUICK_SOLVE\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "QUICK_SOLVE\t\tFALSE\n");
    }
    if (options.SPATIAL_FROST) {
        fprintf(LOG_DEST, "SPATIAL_FROST\t\tTRUE\n");
        fprintf(LOG_DEST, "Nfrost\t\t%zu\n", options.Nfrost);
    }
    else {
        fprintf(LOG_DEST, "SPATIAL_FROST\t\tFALSE\n");
    }
    if (options.SPATIAL_SNOW) {
        fprintf(LOG_DEST, "SPATIAL_SNOW\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "SPATIAL_SNOW\t\tFALSE\n");
    }
    if (options.SNOW_DENSITY == DENS_BRAS) {
        fprintf(LOG_DEST, "SNOW_DENSITY\t\tDENS_BRAS\n");
    }
    else if (options.SNOW_DENSITY == DENS_SNTHRM) {
        fprintf(LOG_DEST, "SNOW_DENSITY\t\tDENS_SNTHRM\n");
    }
    if (options.TFALLBACK) {
        fprintf(LOG_DEST, "TFALLBACK\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "TFALLBACK\t\tFALSE\n");
    }
    fprintf(LOG_DEST, "WIND_H\t\t\t%f\n", global_param.wind_h);
    fprintf(LOG_DEST, "NODES\t\t\t%zu\n", options.Nnode);
    if (options.CARBON) {
        fprintf(LOG_DEST, "CARBON\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "CARBON\t\tFALSE\n");
    }
    if (options.SHARE_LAYER_MOIST) {
        fprintf(LOG_DEST, "SHARE_LAYER_MOIST\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "SHARE_LAYER_MOIST\t\tFALSE\n");
    }
    fprintf(LOG_DEST, "Ncanopy\t\t%zu\n", options.Ncanopy);

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "GRID_DECIMAL\t\t%hu\n", options.GRID_DECIMAL);

    fprintf(LOG_DEST, "\n");
    if (options.BASEFLOW == ARNO) {
        fprintf(LOG_DEST, "BASEFLOW\t\tARNO\n");
    }
    else if (options.BASEFLOW == NIJSSEN2001) {
        fprintf(LOG_DEST, "BASEFLOW\t\tNIJSSEN2001\n");
    }
    if (options.JULY_TAVG_SUPPLIED) {
        fprintf(LOG_DEST, "JULY_TAVG_SUPPLIED\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "JULY_TAVG_SUPPLIED\t\tFALSE\n");
    }
    if (options.ORGANIC_FRACT) {
        fprintf(LOG_DEST, "ORGANIC_FRACT\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "ORGANIC_FRACT\t\tFALSE\n");
    }

    fprintf(LOG_DEST, "\n");
    if (options.VEGLIB_PHOTO) {
        fprintf(LOG_DEST, "VEGLIB_PHOTO\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "VEGLIB_PHOTO\t\tFALSE\n");
    }
    fprintf(LOG_DEST, "ROOT_ZONES\t\t%zu\n", options.ROOT_ZONES);
    if (options.VEGPARAM_LAI) {
        fprintf(LOG_DEST, "VEGPARAM_LAI\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "VEGPARAM_LAI\t\tFALSE\n");
    }
    if (options.LAI_SRC == FROM_VEGHIST) {
        fprintf(LOG_DEST, "LAI_SRC\t\tFROM_VEGHIST\n");
    }
    else if (options.LAI_SRC == FROM_VEGPARAM) {
        fprintf(LOG_DEST, "LAI_SRC\t\tFROM_VEGPARAM\n");
    }
    else if (options.LAI_SRC == FROM_VEGLIB) {
        fprintf(LOG_DEST, "LAI_SRC\t\tFROM_VEGLIB\n");
    }
    if (options.VEGPARAM_FCAN) {
        fprintf(LOG_DEST, "VEGPARAM_FCAN\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "VEGPARAM_FCAN\t\tFALSE\n");
    }
    if (options.FCAN_SRC == FROM_VEGHIST) {
        fprintf(LOG_DEST, "FCAN_SRC\t\tFROM_VEGHIST\n");
    }
    else if (options.FCAN_SRC == FROM_VEGPARAM) {
        fprintf(LOG_DEST, "FCAN_SRC\t\tFROM_VEGPARAM\n");
    }
    else if (options.FCAN_SRC == FROM_VEGLIB) {
        fprintf(LOG_DEST, "FCAN_SRC\t\tFROM_VEGLIB\n");
    }
    else if (options.FCAN_SRC == FROM_DEFAULT) {
        fprintf(LOG_DEST, "FCAN_SRC\t\tFROM_DEFAULT\n");
    }
    if (options.VEGPARAM_ALB) {
        fprintf(LOG_DEST, "VEGPARAM_ALB\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "VEGPARAM_ALB\t\tFALSE\n");
    }
    if (options.ALB_SRC == FROM_VEGHIST) {
        fprintf(LOG_DEST, "ALB_SRC\t\tFROM_VEGHIST\n");
    }
    else if (options.ALB_SRC == FROM_VEGPARAM) {
        fprintf(LOG_DEST, "ALB_SRC\t\tFROM_VEGPARAM\n");
    }
    else if (options.ALB_SRC == FROM_VEGLIB) {
        fprintf(LOG_DEST, "ALB_SRC\t\tFROM_VEGLIB\n");
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Input Elevation Data:\n");
    if (options.SNOW_BAND > 1) {
        fprintf(LOG_DEST, "SNOW_BAND\t\t%zu\n", options.SNOW_BAND);
    }
    else if (options.SNOW_BAND == 1) {
        fprintf(LOG_DEST,
                "SNOW_BAND\t\t%zu\t(no input file needed for SNOW_BAND=1)\n",
                options.SNOW_BAND);
    }
    else {
        fprintf(LOG_DEST, "SNOW_BAND\t\t%zu\n", options.SNOW_BAND);
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Input Lake Data:\n");
    if (options.LAKES) {
        fprintf(LOG_DEST, "LAKES\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "LAKES\t\tFALSE\n");
    }
    if (options.LAKE_PROFILE) {
        fprintf(LOG_DEST, "LAKE_PROFILE\t\tTRUE\n");
    }
    else {
        fprintf(LOG_DEST, "LAKE_PROFILE\t\tFALSE\n");
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Input State File:\n");
    if (options.INIT_STATE) {
        fprintf(LOG_DEST, "INIT_STATE\t\tTRUE\n");
        if (options.STATE_FORMAT == BINARY) {
            fprintf(LOG_DEST, "STATE_FORMAT\tBINARY\n");
        }
        else if (options.STATE_FORMAT == ASCII) {
            fprintf(LOG_DEST, "STATE_FORMAT\tASCII\n");
        }
    }
    else {
        fprintf(LOG_DEST, "INIT_STATE\t\tFALSE\n");
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Output State File:\n");
    if (options.SAVE_STATE) {
        fprintf(LOG_DEST, "SAVE_STATE\t\tTRUE\n");
        fprintf(LOG_DEST, "STATEYEAR\t\t%d\n", global_param.stateyear);
        fprintf(LOG_DEST, "STATEMONTH\t\t%d\n", global_param.statemonth);
        fprintf(LOG_DEST, "STATEDAY\t\t%d\n", global_param.stateday);
        if (options.STATE_FORMAT == BINARY) {
            fprintf(LOG_DEST, "STATE_FORMAT\tBINARY\n");
        }
        else if (options.STATE_FORMAT == ASCII) {
            fprintf(LOG_DEST, "STATE_FORMAT\tASCII\n");
        }
    }
    else {
        fprintf(LOG_DEST, "SAVE_STATE\t\tFALSE\n");
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Output Data:\n");
    fprintf(LOG_DEST, "\n");
}
