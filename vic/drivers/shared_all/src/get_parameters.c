/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model parameters file
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Read the VIC model parameters file
 *****************************************************************************/
void
get_parameters(FILE *paramfile)
{
    char                     cmdstr[MAXSTRING];
    char                     optstr[MAXSTRING];

    extern parameters_struct param;

    /** Read through parameter file to find parameters **/

    rewind(paramfile);
    fgets(cmdstr, MAXSTRING, paramfile);

    while (!feof(paramfile)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#') {
                fgets(cmdstr, MAXSTRING, paramfile);
                continue;
            }

            /*************************************
               Get Model Parameters
            *************************************/
            // Lapse Rate
            if (strcasecmp("LAPSE_RATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAPSE_RATE);
            }
            // Precipitation Guage Height
            else if (strcasecmp("GAUGE_HEIGHT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.GAUGE_HEIGHT);
            }
            // Huge Resistance Term
            else if (strcasecmp("HUGE_RESIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.HUGE_RESIST);
            }
            // Surface Albedo Parameters
            else if (strcasecmp("ALBEDO_BARE_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ALBEDO_BARE_SOIL);
            }
            // Surface Emissivities
            else if (strcasecmp("EMISS_GRND", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_GRND);
            }
            else if (strcasecmp("EMISS_ICE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_ICE);
            }
            else if (strcasecmp("EMISS_VEG", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_VEG);
            }
            else if (strcasecmp("EMISS_SNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_SNOW);
            }
            else if (strcasecmp("EMISS_H2O", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.EMISS_H2O);
            }
            // Soil Constraints
            else if (strcasecmp("SOIL_RESID_MOIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SOIL_RESID_MOIST);
            }
            else if (strcasecmp("SOIL_SLAB_MOIST_FRACT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SOIL_SLAB_MOIST_FRACT);
            }
            // Vegetation Parameters
            else if (strcasecmp("VEG_LAI_SNOW_MULTIPLIER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.VEG_LAI_SNOW_MULTIPLIER);
            }
            else if (strcasecmp("VEG_MIN_INTERCEPTION_STORAGE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.VEG_MIN_INTERCEPTION_STORAGE);
            }
            else if (strcasecmp("VEG_LAI_WATER_FACTOR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.VEG_LAI_WATER_FACTOR);
            }
            // Canopy Parameters
            else if (strcasecmp("CANOPY_CLOSURE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_CLOSURE);
            }
            else if (strcasecmp("CANOPY_RSMAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_RSMAX);
            }
            else if (strcasecmp("CANOPY_VPDMINFACTOR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_VPDMINFACTOR);
            }
            // Lake Parameters
            else if (strcasecmp("LAKE_TMELT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_TMELT);
            }
            else if (strcasecmp("LAKE_MAX_SURFACE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_MAX_SURFACE);
            }
            else if (strcasecmp("LAKE_BETA", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_BETA);
            }
            else if (strcasecmp("LAKE_FRACMIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_FRACMIN);
            }
            else if (strcasecmp("LAKE_FRACLIM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_FRACLIM);
            }
            else if (strcasecmp("LAKE_DM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_DM);
            }
            else if (strcasecmp("LAKE_SNOWCRIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_SNOWCRIT);
            }
            else if (strcasecmp("LAKE_ZWATER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_ZWATER);
            }
            else if (strcasecmp("LAKE_ZSNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_ZSNOW);
            }
            else if (strcasecmp("LAKE_RHOSNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_RHOSNOW);
            }
            else if (strcasecmp("LAKE_CONDI", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_CONDI);
            }
            else if (strcasecmp("LAKE_CONDS", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_CONDS);
            }
            else if (strcasecmp("LAKE_LAMISW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMISW);
            }
            else if (strcasecmp("LAKE_LAMILW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMILW);
            }
            else if (strcasecmp("LAKE_LAMSSW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMSSW);
            }
            else if (strcasecmp("LAKE_LAMSLW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMSLW);
            }
            else if (strcasecmp("LAKE_LAMWSW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMWSW);
            }
            else if (strcasecmp("LAKE_LAMWLW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_LAMWLW);
            }
            else if (strcasecmp("LAKE_A1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_A1);
            }
            else if (strcasecmp("LAKE_A2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_A2);
            }
            else if (strcasecmp("LAKE_QWTAU", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.LAKE_QWTAU);
            }
            else if (strcasecmp("LAKE_MAX_ITER", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.LAKE_MAX_ITER);
            }
            // Saturation Vapor Pressure Parameters
            else if (strcasecmp("SVP_A", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_A);
            }
            else if (strcasecmp("SVP_B", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_B);
            }
            else if (strcasecmp("SVP_C", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SVP_C);
            }
            // Photosynthesis Parameters
            else if (strcasecmp("PHOTO_OMEGA", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_OMEGA);
            }
            else if (strcasecmp("PHOTO_LAIMAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_LAIMAX);
            }
            else if (strcasecmp("PHOTO_LAILIMIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_LAILIMIT);
            }
            else if (strcasecmp("PHOTO_LAIMIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_LAIMIN);
            }
            else if (strcasecmp("PHOTO_EPAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EPAR);
            }
            else if (strcasecmp("PHOTO_FCMAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FCMAX);
            }
            else if (strcasecmp("PHOTO_FCMIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FCMIN);
            }
            else if (strcasecmp("PHOTO_ZENITHMIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ZENITHMIN);
            }
            else if (strcasecmp("PHOTO_ZENITHMINPAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ZENITHMINPAR);
            }
            else if (strcasecmp("PHOTO_ALBSOIPARMIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ALBSOIPARMIN);
            }
            else if (strcasecmp("PHOTO_MINMAXETRANS", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_MINMAXETRANS);
            }
            else if (strcasecmp("PHOTO_MINSTOMCOND", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_MINSTOMCOND);
            }
            else if (strcasecmp("PHOTO_FCI1C3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FCI1C3);
            }
            else if (strcasecmp("PHOTO_FCI1C4", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FCI1C4);
            }
            else if (strcasecmp("PHOTO_OX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_OX);
            }
            else if (strcasecmp("PHOTO_KC", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_KC);
            }
            else if (strcasecmp("PHOTO_KO", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_KO);
            }
            else if (strcasecmp("PHOTO_EC", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EC);
            }
            else if (strcasecmp("PHOTO_EO", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EO);
            }
            else if (strcasecmp("PHOTO_EV", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EV);
            }
            else if (strcasecmp("PHOTO_ER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ER);
            }
            else if (strcasecmp("PHOTO_ALC3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ALC3);
            }
            else if (strcasecmp("PHOTO_FRDC3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FRDC3);
            }
            else if (strcasecmp("PHOTO_EK", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_EK);
            }
            else if (strcasecmp("PHOTO_ALC4", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_ALC4);
            }
            else if (strcasecmp("PHOTO_FRDC4", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FRDC4);
            }
            else if (strcasecmp("PHOTO_THETA", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_THETA);
            }
            else if (strcasecmp("PHOTO_FRLEAF", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FRLEAF);
            }
            else if (strcasecmp("PHOTO_FRGROWTH", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.PHOTO_FRGROWTH);
            }
            // Soil Respiration Parameters
            else if (strcasecmp("SRESP_E0_LT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_E0_LT);
            }
            else if (strcasecmp("SRESP_T0_LT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_T0_LT);
            }
            else if (strcasecmp("SRESP_WMINFM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_WMINFM);
            }
            else if (strcasecmp("SRESP_WMAXFM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_WMAXFM);
            }
            else if (strcasecmp("SRESP_WOPTFM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_WOPTFM);
            }
            else if (strcasecmp("SRESP_RHSAT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_RHSAT);
            }
            else if (strcasecmp("SRESP_RFACTOR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_RFACTOR);
            }
            else if (strcasecmp("SRESP_TAULITTER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_TAULITTER);
            }
            else if (strcasecmp("SRESP_TAUINTER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_TAUINTER);
            }
            else if (strcasecmp("SRESP_TAUSLOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_TAUSLOW);
            }
            else if (strcasecmp("SRESP_FAIR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_FAIR);
            }
            else if (strcasecmp("SRESP_FINTER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SRESP_FINTER);
            }
            // Snow Parameters
            else if (strcasecmp("SNOW_MAX_SURFACE_SWE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_MAX_SURFACE_SWE);
            }
            else if (strcasecmp("SNOW_LIQUID_WATER_CAPACITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_LIQUID_WATER_CAPACITY);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENSITY", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENSITY);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENS_MAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENS_MAX);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_DENS_MAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_DENS_MAX);
            }
            else if (strcasecmp("SNOW_DENS_DMLIMIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_DMLIMIT);
            }
            else if (strcasecmp("SNOW_DENS_DMLIMIT_FACTOR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_DMLIMIT_FACTOR);
            }
            else if (strcasecmp("SNOW_DENS_MAX_CHANGE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_MAX_CHANGE);
            }
            else if (strcasecmp("SNOW_DENS_ETA0", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_ETA0);
            }
            else if (strcasecmp("SNOW_DENS_C1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C1);
            }
            else if (strcasecmp("SNOW_DENS_C2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C2);
            }
            else if (strcasecmp("SNOW_DENS_C3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C3);
            }
            else if (strcasecmp("SNOW_DENS_C3_CONST", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C3_CONST);
            }
            else if (strcasecmp("SNOW_DENS_C4", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C4);
            }
            else if (strcasecmp("SNOW_DENS_C4WET", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C4WET);
            }
            else if (strcasecmp("SNOW_DENS_C5", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C5);
            }
            else if (strcasecmp("SNOW_DENS_C6", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C6);
            }
            else if (strcasecmp("SNOW_DENS_F", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_F);
            }
            else if (strcasecmp("SNOW_DENS_EXP", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_EXP);
            }
            else if (strcasecmp("SNOW_DENS_DENOM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_DENOM);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C1);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C2);
            }
            else if (strcasecmp("SNOW_NEW_SNT_C3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNT_C3);
            }
            else if (strcasecmp("SNOW_NEW_BRAS_DENOM", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_BRAS_DENOM);
            }
            else if (strcasecmp("SNOW_MIN_SWQ_EB_THRES", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_MIN_SWQ_EB_THRES);
            }
            else if (strcasecmp("SNOW_A1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_A1);
            }
            else if (strcasecmp("SNOW_A2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_A2);
            }
            else if (strcasecmp("SNOW_L1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_L1);
            }
            else if (strcasecmp("SNOW_L2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_L2);
            }
            else if (strcasecmp("SNOW_NEW_SNOW_ALB", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_NEW_SNOW_ALB);
            }
            else if (strcasecmp("SNOW_ALB_ACCUM_A", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_ALB_ACCUM_A);
            }
            else if (strcasecmp("SNOW_ALB_ACCUM_B", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_ALB_ACCUM_B);
            }
            else if (strcasecmp("SNOW_ALB_THAW_A", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_ALB_THAW_A);
            }
            else if (strcasecmp("SNOW_ALB_THAW_B", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_ALB_THAW_B);
            }
            else if (strcasecmp("SNOW_TRACESNOW", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_TRACESNOW);
            }
            else if (strcasecmp("SNOW_CONDUCT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_CONDUCT);
            }
            else if (strcasecmp("SNOW_MAX_SNOW_TEMP", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_MAX_SNOW_TEMP);
            }
            else if (strcasecmp("SNOW_MIN_RAIN_TEMP", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_MIN_RAIN_TEMP);
            }
            // Blowing Snow Parameters
            else if (strcasecmp("BLOWING_KA", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.BLOWING_KA);
            }
            else if (strcasecmp("BLOWING_CSALT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.BLOWING_CSALT);
            }
            else if (strcasecmp("BLOWING_UTHRESH", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.BLOWING_UTHRESH);
            }
            else if (strcasecmp("BLOWING_KIN_VIS", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.BLOWING_KIN_VIS);
            }
            else if (strcasecmp("BLOWING_MAX_ITER", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.BLOWING_MAX_ITER);
            }
            else if (strcasecmp("BLOWING_K", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.BLOWING_K);
            }
            else if (strcasecmp("BLOWING_SETTLING", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.BLOWING_SETTLING);
            }
            else if (strcasecmp("BLOWING_NUMINCS", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.BLOWING_NUMINCS);
            }
            // Treeline temperature
            else if (strcasecmp("TREELINE_TEMPERATURE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.TREELINE_TEMPERATURE);
            }
            // Iteration Bracket Widths
            else if (strcasecmp("SNOW_DT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DT);
            }
            else if (strcasecmp("SURF_DT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SURF_DT);
            }
            else if (strcasecmp("SOIL_DT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SOIL_DT);
            }
            else if (strcasecmp("CANOPY_DT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_DT);
            }
            else if (strcasecmp("CANOPY_VP", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CANOPY_VP);
            }
            // Convergence Tolerances
            else if (strcasecmp("TOL_GRND", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.TOL_GRND);
            }
            else if (strcasecmp("TOL_OVER", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.TOL_OVER);
            }
            // Frozen Soil Parameters
            else if (strcasecmp("FROZEN_MAXITER", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.FROZEN_MAXITER);
            }
            // Canopy Iterations
            else if (strcasecmp("MAX_ITER_GRND_CANOPY", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.MAX_ITER_GRND_CANOPY);
            }
            // Newton-Raphson Solver Parameters
            else if (strcasecmp("NEWT_RAPH_MAXTRIAL", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.NEWT_RAPH_MAXTRIAL);
            }
            else if (strcasecmp("NEWT_RAPH_TOLX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_TOLX);
            }
            else if (strcasecmp("NEWT_RAPH_TOLF", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_TOLF);
            }
            else if (strcasecmp("NEWT_RAPH_R_MAX", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_R_MAX);
            }
            else if (strcasecmp("NEWT_RAPH_R_MIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_R_MIN);
            }
            else if (strcasecmp("NEWT_RAPH_RELAX1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_RELAX1);
            }
            else if (strcasecmp("NEWT_RAPH_RELAX2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_RELAX2);
            }
            else if (strcasecmp("NEWT_RAPH_RELAX3", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_RELAX3);
            }
            else if (strcasecmp("NEWT_RAPH_EPS2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.NEWT_RAPH_EPS2);
            }
            // Root-Brent parameters
            else if (strcasecmp("ROOT_BRENT_MAXTRIES", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.ROOT_BRENT_MAXTRIES);
            }
            else if (strcasecmp("ROOT_BRENT_MAXITER", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &param.ROOT_BRENT_MAXITER);
            }
            else if (strcasecmp("ROOT_BRENT_TSTEP", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ROOT_BRENT_TSTEP);
            }
            else if (strcasecmp("ROOT_BRENT_T", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ROOT_BRENT_T);
            }
            else {
                log_warn("Unrecognized option in the parameter file:  %s "
                         "- check your spelling", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, paramfile);
    }
}

/******************************************************************************
 * @brief    Validate VIC model parameter values
 *****************************************************************************/
void
validate_parameters()
{
    extern parameters_struct param;

    // Validate Parameters
    // Lapse Rate
    if (!(param.LAPSE_RATE >= -1 && param.LAPSE_RATE <= 0)) {
        log_err("LAPSE_RATE must be defined on the interval [-1,0] (C/m)")
    }
    // Precipitation Guage Height
    if (!(param.GAUGE_HEIGHT >= 0 && param.GAUGE_HEIGHT <= 100)) {
        log_err("GAUGE_HEIGHT must be defined on the interval [0,100] (m)")
    }
    // Huge Resistance Term
    if (!(param.HUGE_RESIST >= 0.)) {
        log_err("HUGE_RESIST must be defined on the interval [0, inf) (s/m)");
    }
    // Surface Albedo Parameters
    if (!(param.ALBEDO_BARE_SOIL >= 0 && param.ALBEDO_BARE_SOIL <= 1)) {
        log_err("ALBEDO_BARE_SOIL must be defined on the interval [0,1] (-)")
    }
    // Surface Emissivities
    if (!(param.EMISS_GRND >= 0 && param.EMISS_GRND <= 1)) {
        log_err("EMISS_GRND must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_ICE >= 0 && param.EMISS_ICE <= 1)) {
        log_err("EMISS_ICE must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_VEG >= 0 && param.EMISS_VEG <= 1)) {
        log_err("EMISS_VEG must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_SNOW >= 0 && param.EMISS_SNOW <= 1)) {
        log_err("EMISS_SNOW must be defined on the interval [0,1] (-)")
    }
    if (!(param.EMISS_H2O >= 0 && param.EMISS_H2O <= 1)) {
        log_err("EMISS_H2O must be defined on the interval [0,1] (-)")
    }
    // Soil Constraints
    if (!(param.SOIL_RESID_MOIST >= 0.)) {
        log_err("SOIL_RESID_MOIST must be defined on the interval [0, inf)");
    }
    if (!(param.SOIL_SLAB_MOIST_FRACT >= 0 && param.SOIL_SLAB_MOIST_FRACT <=
          1)) {
        log_err(
            "SOIL_SLAB_MOIST_FRACT must be defined on the interval [0,1] (-)")
    }
    // Vegetation Parameters
    if (!(param.VEG_LAI_SNOW_MULTIPLIER >= 0.)) {
        log_err(
            "VEG_LAI_SNOW_MULTIPLIER must be defined on the interval [0, inf) (-)");
    }
    if (!(param.VEG_MIN_INTERCEPTION_STORAGE >= 0.)) {
        log_err(
            "VEG_MIN_INTERCEPTION_STORAGE must be defined on the interval [0, inf) (m)");
    }
    if (!(param.VEG_LAI_WATER_FACTOR >= 0.)) {
        log_err(
            "VEG_LAI_WATER_FACTOR must be defined on the interval [0, inf) (-)");
    }
    // Canopy Parameters
    if (!(param.CANOPY_CLOSURE >= 0.)) {
        log_err("CANOPY_CLOSURE must be defined on the interval [0, inf) (Pa)");
    }
    if (!(param.CANOPY_RSMAX >= 0.)) {
        log_err("CANOPY_RSMAX must be defined on the interval [0, inf) (s/m)");
    }
    if (!(param.CANOPY_VPDMINFACTOR >= 0.)) {
        log_err(
            "CANOPY_VPDMINFACTOR must be defined on the interval [0, inf) (-)");
    }

    // Lake Parameters
    // LAKE_TMELT - Currently, no constraints
    if (!(param.LAKE_MAX_SURFACE >= 0.)) {
        log_err("LAKE_MAX_SURFACE must be defined on the interval [0, inf) (m)");
    }
    // LAKE_BETA - Currently, no constraints
    // LAKE_FRACMIN - Currently, no constraints
    // LAKE_FRACLIM - Currently, no constraints
    // LAKE_DM - Currently, no constraints
    if (!(param.LAKE_SNOWCRIT >= 0.)) {
        log_err("LAKE_SNOWCRIT must be defined on the interval [0, inf) (m)");
    }
    // LAKE_ZWATER - Currently, no constraints
    // LAKE_ZSNOW - Currently, no constraints
    if (!(param.LAKE_RHOSNOW >= 0.)) {
        log_err("LAKE_RHOSNOW must be defined on the interval [0, inf) (kg m-3)");
    }
    // LAKE_CONDI - Currently, no constraints
    // LAKE_CONDS - Currently, no constraints
    // LAKE_LAMISW - Currently, no constraints
    // LAKE_LAMILW - Currently, no constraints
    // LAKE_LAMSSW - Currently, no constraints
    // LAKE_LAMSLW - Currently, no constraints
    // LAKE_LAMWSW - Currently, no constraints
    // LAKE_LAMWLW - Currently, no constraints
    // LAKE_A1 - Currently, no constraints
    // LAKE_A2 - Currently, no constraints
    // LAKE_QWTAU - Currently, no constraints
    if (!(param.LAKE_MAX_ITER >= 0.)) {
        log_err(
            "LAKE_MAX_ITER must be defined on the interval [0, inf) (iterations)");
    }
    // Saturation Vapor Pressure Parameters
    if (!(param.SVP_A >= 0.)) {
        log_err("SVP_A must be defined on the interval [0, inf) (kPa)");
    }
    if (!(param.SVP_B >= 0.)) {
        log_err("SVP_B must be defined on the interval [0, inf) (kPa)");
    }
    if (!(param.SVP_C >= 0.)) {
        log_err("SVP_C must be defined on the interval [0, inf) (kPa)");
    }
    // Photosynthesis Parameters
    // PHOTO_OMEGA - Currently, no constraints
    if (!(param.PHOTO_LAIMAX >= 0.)) {
        log_err("PHOTO_LAIMAX must be defined on the interval [0, inf) (-)");
    }
    // PHOTO_LAILIMIT - Currently, no constraints
    if (!(param.PHOTO_LAILIMIT >= 0.)) {
        log_err("PHOTO_LAILIMIT must be defined on the interval [0, inf) (-)");
    }
    if (!(param.PHOTO_LAIMIN >= 0.)) {
        log_err("PHOTO_LAIMIN must be defined on the interval [0, inf) (-)");
    }
    if (!(param.PHOTO_EPAR >= 0.)) {
        log_err(
            "PHOTO_EPAR must be defined on the interval [0, inf) (mol/MJ PAR)");
    }
    if (!(param.PHOTO_FCMAX >= 0 && param.PHOTO_FCMAX <= 1)) {
        log_err("PHOTO_FCMAX must be defined on the interval [0,1] (-)")
    }
    if (!(param.PHOTO_FCMIN >= 0 && param.PHOTO_FCMIN <= 1)) {
        log_err("PHOTO_FCMIN must be defined on the interval [0,1] (-)")
    }
    // PHOTO_ZENITHMIN - Currently, no constraints
    // PHOTO_ZENITHMINPAR - Currently, no constraints
    if (!(param.PHOTO_ALBSOIPARMIN >= 0 && param.PHOTO_ALBSOIPARMIN <= 1)) {
        log_err("PHOTO_ALBSOIPARMIN must be defined on the interval [0,1] (-)")
    }
    if (!(param.PHOTO_MINMAXETRANS >= 0.)) {
        log_err(
            "PHOTO_MINMAXETRANS must be defined on the interval [0, inf) (mol/(m^2 s))");
    }
    if (!(param.PHOTO_MINSTOMCOND >= 0.)) {
        log_err(
            "PHOTO_MINSTOMCOND must be defined on the interval [0, inf) (mol H2O/m2s)");
    }
    // PHOTO_FCI1C3 - Currently, no constraints
    // PHOTO_FCI1C4 - Currently, no constraints
    if (!(param.PHOTO_OX >= 0.)) {
        log_err(
            "PHOTO_OX must be defined on the interval [0, inf) (mol H2O/m2s)");
    }
    // PHOTO_KC - Currently, no constraints
    // PHOTO_KO - Currently, no constraints
    // PHOTO_EC - Currently, no constraints
    // PHOTO_EO - Currently, no constraints
    // PHOTO_EV - Currently, no constraints
    // PHOTO_ER - Currently, no constraints
    // PHOTO_ALC3 - Currently, no constraints
    // PHOTO_FRDC3 - Currently, no constraints
    // PHOTO_EK - Currently, no constraints
    // PHOTO_ALC4 - Currently, no constraints
    // PHOTO_FRDC4 - Currently, no constraints
    // PHOTO_THETA - Currently, no constraints
    // PHOTO_FRLEAF - Currently, no constraints
    // PHOTO_FRGROWTH - Currently, no constraints

    // Soil Respiration Parameters
    if (!(param.SRESP_E0_LT >= 0.)) {
        log_err(
            "SRESP_E0_LT must be defined on the interval [0, inf) (mol H2O/m2s)");
    }
    if (!(param.SRESP_T0_LT >= 0.)) {
        log_err(
            "SRESP_T0_LT must be defined on the interval [0, inf) (mol H2O/m2s)");
    }
    if (!(param.SRESP_WMINFM >= 0 && param.SRESP_WMINFM <= 1)) {
        log_err("SRESP_WMINFM must be defined on the interval [0,1] (-)")
    }
    if (!(param.SRESP_WMAXFM >= 0 && param.SRESP_WMAXFM <= 1)) {
        log_err("SRESP_WMAXFM must be defined on the interval [0,1] (-)")
    }
    if (!(param.SRESP_WOPTFM >= 0 && param.SRESP_WOPTFM <= 1)) {
        log_err("SRESP_WOPTFM must be defined on the interval [0,1] (-)")
    }
    // SRESP_RHSAT - Currently, no constraints
    // SRESP_RFACTOR - Currently, no constraints
    if (!(param.SRESP_TAULITTER >= 0.)) {
        log_err("SRESP_TAULITTER must be defined on the interval [0, inf) (y)");
    }
    if (!(param.SRESP_TAUINTER >= 0.)) {
        log_err("SRESP_TAUINTER must be defined on the interval [0, inf) (y)");
    }
    if (!(param.SRESP_TAUSLOW >= 0.)) {
        log_err("SRESP_TAUSLOW must be defined on the interval [0, inf) (y)");
    }
    if (!(param.SRESP_FAIR >= 0 && param.SRESP_FAIR <= 1)) {
        log_err("SRESP_FAIR must be defined on the interval [0,1] (-)")
    }
    if (!(param.SRESP_FINTER >= 0 && param.SRESP_FINTER <= 1)) {
        log_err("SRESP_FINTER must be defined on the interval [0,1] (-)")
    }

    // Snow Parameters
    if (!(param.SNOW_MAX_SURFACE_SWE >= 0.)) {
        log_err(
            "SNOW_MAX_SURFACE_SWE must be defined on the interval [0, inf) (m)");
    }
    if (!(param.SNOW_LIQUID_WATER_CAPACITY >= 0 &&
          param.SNOW_LIQUID_WATER_CAPACITY <= 1)) {
        log_err(
            "SNOW_LIQUID_WATER_CAPACITY must be defined on the interval [0,1] (-)")
    }
    if (!(param.SNOW_NEW_SNOW_DENSITY >= 0.)) {
        log_err(
            "SNOW_NEW_SNOW_DENSITY must be defined on the interval [0, inf) (kg/m^3)");
    }
    if (!(param.SNOW_DEPTH_THRES >= 0.)) {
        log_err(
            "SNOW_DEPTH_THRES must be defined on the interval [0, inf) (m)");
    }
    if (!(param.SNOW_DENS_DMLIMIT >= 0.)) {
        log_err(
            "SNOW_DENS_DMLIMIT must be defined on the interval [0, inf) (kg/m^3)");
    }
    if (!(param.SNOW_NEW_SNOW_DENS_MAX >= 0. &&
          param.SNOW_NEW_SNOW_DENS_MAX <= 700.)) {
        log_err(
            "SNOW_NEW_SNOW_DENS_MAX must be defined on the interval [0, 700) (kg/m^3)");
    }
    if (!(param.SNOW_DENS_MAX_CHANGE >= 0 && param.SNOW_DENS_MAX_CHANGE <= 1)) {
        log_err("SNOW_DENS_MAX_CHANGE must be defined on the interval [0,1] (-)")
    }
    // SNOW_DENS_ETA0 - Currently, no constraints
    // SNOW_DENS_C1 - Currently, no constraints
    // SNOW_DENS_C2 - Currently, no constraints
    // SNOW_DENS_C5 - Currently, no constraints
    // SNOW_DENS_C6 - Currently, no constraints
    // SNOW_DENS_F - Currently, no constraints
    if (!(param.SNOW_MIN_SWQ_EB_THRES >= 0.)) {
        log_err(
            "SNOW_MIN_SWQ_EB_THRES must be defined on the interval [0, inf) (m)");
    }
    if (!(param.SNOW_A1 >= 0.)) {
        log_err("SNOW_A1 must be defined on the interval [0, inf)");
    }
    if (!(param.SNOW_A2 >= 0.)) {
        log_err("SNOW_A2 must be defined on the interval [0, inf)");
    }
    if (!(param.SNOW_L1 >= 0.)) {
        log_err("SNOW_L1 must be defined on the interval [0, inf) (1/m)");
    }
    if (!(param.SNOW_L2 >= 0.)) {
        log_err("SNOW_L2 must be defined on the interval [0, inf) (1/m)");
    }
    if (!(param.SNOW_NEW_SNOW_ALB >= 0 && param.SNOW_NEW_SNOW_ALB <= 1)) {
        log_err("SNOW_NEW_SNOW_ALB must be defined on the interval [0,1] (-)");
    }
    if (!(param.SNOW_ALB_ACCUM_A >= 0.)) {
        log_err("SNOW_ALB_ACCUM_A must be defined on the interval [0, inf) (-)");
    }
    if (!(param.SNOW_ALB_ACCUM_B >= 0.)) {
        log_err("SNOW_ALB_ACCUM_B must be defined on the interval [0, inf) (-)");
    }
    if (!(param.SNOW_ALB_THAW_A >= 0.)) {
        log_err("SNOW_ALB_THAW_A must be defined on the interval [0, inf) (-)");
    }
    if (!(param.SNOW_ALB_THAW_B >= 0.)) {
        log_err("SNOW_ALB_THAW_B must be defined on the interval [0, inf) (-)");
    }
    if (!(param.SNOW_TRACESNOW >= 0.)) {
        log_err("SNOW_TRACESNOW must be defined on the interval [0, inf) (mm)");
    }
    if (!(param.SNOW_CONDUCT >= 0.)) {
        log_err("SNOW_CONDUCT must be defined on the interval [0, inf) (W/mK)");
    }
    if (!(param.SNOW_MAX_SNOW_TEMP >= -10 && param.SNOW_MAX_SNOW_TEMP <= 10)) {
        log_err(
            "SNOW_MAX_SNOW_TEMP must be defined on the interval [-10,10] (C)");
    }
    if (!(param.SNOW_MIN_RAIN_TEMP >= -10 && param.SNOW_MIN_RAIN_TEMP <= 10)) {
        log_err(
            "SNOW_MIN_RAIN_TEMP must be defined on the interval [-10,10] (C)");
    }
    if (!(param.SNOW_MIN_RAIN_TEMP < param.SNOW_MAX_SNOW_TEMP)) {
        log_err("SNOW_MIN_RAIN_TEMP > SNOW_MAX_SNOW_TEMP.");
    }
    // Blowing Snow Parameters
    if (!(param.BLOWING_KA >= 0.)) {
        log_err("BLOWING_KA must be defined on the interval [0, inf) (W/mK)");
    }
    if (!(param.BLOWING_CSALT >= 0.)) {
        log_err("BLOWING_CSALT must be defined on the interval [0, inf) (m/s)");
    }
    if (!(param.BLOWING_UTHRESH >= 0.)) {
        log_err(
            "BLOWING_UTHRESH must be defined on the interval [0, inf)  (m/s)");
    }
    if (!(param.BLOWING_KIN_VIS >= 0.)) {
        log_err(
            "BLOWING_KIN_VIS must be defined on the interval [0, inf)  (m2/s)");
    }
    if (!(param.BLOWING_MAX_ITER >= 1)) {
        log_err(
            "BLOWING_MAX_ITER must be defined on the interval [1, inf) (iterations");
    }
    if (!(param.BLOWING_K >= 0)) {
        log_err("BLOWING_K must be defined on the interval [0, inf)");
    }
    if (!(param.BLOWING_SETTLING >= 0.)) {
        log_err(
            "BLOWING_SETTLING must be defined on the interval [0, inf) (m/s)");
    }
    if (param.BLOWING_NUMINCS < 0) {
        log_err(
            "BLOWING_NUMINCS must be defined on the interval [0, inf) (intervals)");
    }
    // Treeline temperature
    if (!(param.TREELINE_TEMPERATURE >= -10 && param.TREELINE_TEMPERATURE <=
          20)) {
        log_warn(
            "TREELINE_TEMPERATURE must be defined on the interval [-10,20] (C)");
    }

    // Iteration Bracket Widths
    if (!(param.SNOW_DT >= 0.)) {
        log_err("SNOW_DT must be defined on the interval [0, inf) (C)");
    }
    if (!(param.SURF_DT >= 0.)) {
        log_err("SURF_DT must be defined on the interval [0, inf) (C)");
    }
    if (!(param.SOIL_DT >= 0.)) {
        log_err("SOIL_DT must be defined on the interval [0, inf) (C)");
    }
    if (!(param.CANOPY_DT >= 0.)) {
        log_err("CANOPY_DT must be defined on the interval [0, inf) (C)");
    }
    if (!(param.CANOPY_VP >= 0.)) {
        log_err("CANOPY_VP must be defined on the interval [0, inf) (Pa)");
    }
    // Convergence Tolerances
    if (!(param.TOL_GRND >= 0.)) {
        log_err("TOL_GRND must be defined on the interval [0, inf)");
    }
    if (!(param.TOL_OVER >= 0.)) {
        log_err("TOL_OVER must be defined on the interval [0, inf)");
    }
    // Frozen Soil Parameters
    if (!(param.FROZEN_MAXITER >= 0)) {
        log_err(
            "FROZEN_MAXITER must be defined on the interval [0, inf) (iterations");
    }
    // Canopy Iterations
    if (!(param.MAX_ITER_GRND_CANOPY >= 0)) {
        log_err(
            "MAX_ITER_GRND_CANOPY  must be defined on the interval [0, inf) (iterations");
    }
    // Newton-Raphson Solver Parameters
    if (!(param.NEWT_RAPH_MAXTRIAL >= 0)) {
        log_err(
            "NEWT_RAPH_MAXTRIAL must be defined on the interval [0, inf) (trials)");
    }
    if (!(param.NEWT_RAPH_TOLX >= 0.)) {
        log_err("NEWT_RAPH_TOLX must be defined on the interval [0, inf)");
    }
    if (!(param.NEWT_RAPH_TOLF >= 0.)) {
        log_err("NEWT_RAPH_TOLF must be defined on the interval [0, inf)");
    }
    // NEWT_RAPH_R_MAX - Currently, no constraints
    // NEWT_RAPH_R_MIN - Currently, no constraints
    if (!(param.NEWT_RAPH_RELAX1 >= 0.)) {
        log_err("NEWT_RAPH_RELAX1 must be defined on the interval [0, inf)");
    }
    if (!(param.NEWT_RAPH_RELAX2 >= 0.)) {
        log_err("NEWT_RAPH_RELAX2 must be defined on the interval [0, inf)");
    }
    if (!(param.NEWT_RAPH_RELAX3 >= 0.)) {
        log_err("NEWT_RAPH_RELAX3 must be defined on the interval [0, inf)");
    }
    if (!(param.NEWT_RAPH_EPS2 >= 0.)) {
        log_err("NEWT_RAPH_EPS2 must be defined on the interval [0, inf) (-)");
    }
    // Root-Brent parameters
    if (param.ROOT_BRENT_MAXTRIES < 0) {
        log_err("ROOT_BRENT_MAXTRIES must be defined on the interval [0, inf)");
    }
    if (param.ROOT_BRENT_MAXITER < 0) {
        log_err("ROOT_BRENT_MAXITER must be defined on the interval [0, inf)");
    }
    if (param.ROOT_BRENT_TSTEP < 0) {
        log_err("ROOT_BRENT_TSTEP must be defined on the interval [0, inf)");
    }
    if (!(param.ROOT_BRENT_T >= 0.)) {
        log_err("ROOT_BRENT_T must be defined on the interval [0, inf)");
    }
}
