/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model parameters file
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

#include <vic_driver_shared.h>

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
            // Default Wind Speed
            else if (strcasecmp("WIND_SPEED_DEFAULT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.WIND_SPEED_DEFAULT);
            }
            else if (strcasecmp("WIND_SPEED_MIN", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.WIND_SPEED_MIN);
            }
            // Huge Resistance Term
            else if (strcasecmp("HUGE_RESIST", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.HUGE_RESIST);
            }
            // Surface Albedo Parameters
            else if (strcasecmp("ALBEDO_BARE_SOIL", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ALBEDO_BARE_SOIL);
            }
            else if (strcasecmp("ALBEDO_H20_SURF", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.ALBEDO_H20_SURF);
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
            // MTCLIM Parameters
            else if (strcasecmp("MTCLIM_TDAYCOEF", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_TDAYCOEF);
            }
            else if (strcasecmp("MTCLIM_SOLAR_CONSTANT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SOLAR_CONSTANT);
            }
            else if (strcasecmp("MTCLIM_SNOW_TCRIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SNOW_TCRIT);
            }
            else if (strcasecmp("MTCLIM_SNOW_TRATE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SNOW_TRATE);
            }
            else if (strcasecmp("MTCLIM_TBASE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_TBASE);
            }
            else if (strcasecmp("MTCLIM_ABASE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_ABASE);
            }
            else if (strcasecmp("MTCLIM_C", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_C);
            }
            else if (strcasecmp("MTCLIM_B0", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_B0);
            }
            else if (strcasecmp("MTCLIM_B1", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_B1);
            }
            else if (strcasecmp("MTCLIM_B2", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_B2);
            }
            else if (strcasecmp("MTCLIM_RAIN_SCALAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_RAIN_SCALAR);
            }
            else if (strcasecmp("MTCLIM_DIF_ALB", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_DIF_ALB);
            }
            else if (strcasecmp("MTCLIM_SC_INT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SC_INT);
            }
            else if (strcasecmp("MTCLIM_SC_SLOPE", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SC_SLOPE);
            }
            else if (strcasecmp("MTCLIM_SRADDT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SRADDT);
            }
            else if (strcasecmp("MTCLIM_SW_PREC_THRESH", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.MTCLIM_SW_PREC_THRESH);
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
            // Carbon Parameters
            else if (strcasecmp("CARBON_CATMCURRENT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CARBON_CATMCURRENT);
            }
            else if (strcasecmp("CARBON_SW2PAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.CARBON_SW2PAR);
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
            else if (strcasecmp("SNOW_DENS_DMLIMIT", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_DMLIMIT);
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
            else if (strcasecmp("SNOW_DENS_C5", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C5);
            }
            else if (strcasecmp("SNOW_DENS_C6", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_C6);
            }
            else if (strcasecmp("SNOW_DENS_F", optstr) == 0) {
                sscanf(cmdstr, "%*s %lf", &param.SNOW_DENS_F);
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
                         "- check your spelling\n", optstr);
            }
        }
        fgets(cmdstr, MAXSTRING, paramfile);
    }
}
