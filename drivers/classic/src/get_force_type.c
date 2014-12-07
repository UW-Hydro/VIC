/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine determines the current forcing file data type and stores its
 * location in the description of the current forcing file.
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
 * @brief    This routine determines the current forcing file data type and
 *           stores its location in the description of the current forcing file.
 *****************************************************************************/
void
get_force_type(char *cmdstr,
               int   file_num,
               int  *field)
{
    extern param_set_struct param_set;

    char                    optstr[50];
    char                    flgstr[10];
    int                     type;

    type = SKIP;

    /** Initialize flgstr **/
    strcpy(flgstr, "NULL");

    if ((*field) >= (int)param_set.N_TYPES[file_num]) {
        log_err("Too many variables defined for forcing file %i.", file_num);
    }

    sscanf(cmdstr, "%*s %s", optstr);

    /***************************************
       Get meteorological data forcing info
    ***************************************/

    /* type 0: air temperature [C] (ALMA_INPUT: [K]) */
    if (strcasecmp("AIR_TEMP", optstr) == 0) {
        type = AIR_TEMP;
    }
    /* type 1: albedo [fraction] */
    else if (strcasecmp("ALBEDO", optstr) == 0) {
        type = ALBEDO;
    }
    /* type 2: atmospheric CO2 mixing ratio [ppm] */
    else if (strcasecmp("CATM", optstr) == 0) {
        type = CATM;
    }
    /* type 3: incoming channel flow [m3] (ALMA_INPUT: [m3/s]) */
    else if (strcasecmp("CHANNEL_IN", optstr) == 0) {
        type = CHANNEL_IN;
    }
    /* type 4: convective rainfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("CRAINF", optstr) == 0) {
        type = CRAINF;
    }
    /* type 5: convective snowfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("CSNOWF", optstr) == 0) {
        type = CSNOWF;
    }
    /* type 6: air density [kg/m3] */
    else if (strcasecmp("DENSITY", optstr) == 0) {
        type = DENSITY;
    }
    /* type 7: direct fraction of shortwave [fraction] */
    else if (strcasecmp("FDIR", optstr) == 0) {
        type = FDIR;
    }
    /* type 8: LAI [m2/m2] */
    else if (strcasecmp("LAI_IN", optstr) == 0) {
        type = LAI_IN;
    }
    /* type 9: incoming longwave radiation [W/m2] */
    else if (strcasecmp("LONGWAVE",
                        optstr) == 0 || strcasecmp("LWDOWN", optstr) == 0) {
        type = LONGWAVE;
    }
    /* type 10: large-scale rainfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("LSRAINF", optstr) == 0) {
        type = LSRAINF;
    }
    /* type 11: large-scale snowfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("LSSNOWF", optstr) == 0) {
        type = LSSNOWF;
    }
    /* type 12: photosynthetically active radiation [uE/m2s] */
    else if (strcasecmp("PAR", optstr) == 0) {
        type = PAR;
    }
    /* type 13: precipitation [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("PREC", optstr) == 0) {
        type = PREC;
    }
    /* type 14: air pressure [kPa] (ALMA_INPUT: [Pa]) */
    else if (strcasecmp("PRESSURE", optstr) == 0) {
        type = PRESSURE;
    }
    /* type 15: specific humidity [kg/kg] */
    else if (strcasecmp("QAIR", optstr) == 0) {
        type = QAIR;
    }
    /* type 16: rainfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("RAINF", optstr) == 0) {
        type = RAINF;
    }
    /* type 17: relative humidity [fraction] */
    else if (strcasecmp("REL_HUMID", optstr) == 0) {
        type = REL_HUMID;
    }
    /* type 18: shortwave radiation [W/m2] */
    else if (strcasecmp("SHORTWAVE",
                        optstr) == 0 || strcasecmp("SWDOWN", optstr) == 0) {
        type = SHORTWAVE;
    }
    /* type 19: snowfall [mm] (ALMA_INPUT: [mm/s]) */
    else if (strcasecmp("SNOWF", optstr) == 0) {
        type = SNOWF;
    }
    /* type 20: maximum daily temperature [C] (ALMA_INPUT: [K]) */
    else if (strcasecmp("TMAX", optstr) == 0) {
        type = TMAX;
    }
    /* type 21: minimum daily temperature [C] (ALMA_INPUT: [K]) */
    else if (strcasecmp("TMIN", optstr) == 0) {
        type = TMIN;
    }
    /* type 22: cloud cover fraction */
    else if (strcasecmp("TSKC", optstr) == 0) {
        type = TSKC;
    }
    /* type 23: vegetation cover fraction */
    else if (strcasecmp("VEGCOVER", optstr) == 0) {
        type = VEGCOVER;
    }
    /* type 24: vapor pressure [kPa] (ALMA_INPUT: [Pa]) */
    else if (strcasecmp("VP", optstr) == 0) {
        type = VP;
    }
    /* type 25: wind speed [m/s] */
    else if (strcasecmp("WIND", optstr) == 0) {
        type = WIND;
    }
    /* type 26: zonal component of wind speed [m/s] */
    else if (strcasecmp("WIND_E", optstr) == 0) {
        type = WIND_E;
    }
    /* type 27: meridional component of wind speed [m/s] */
    else if (strcasecmp("WIND_N", optstr) == 0) {
        type = WIND_N;
    }
    /* type 28: unused (blank) data */
    else if (strcasecmp("SKIP", optstr) == 0) {
        type = SKIP;
    }
    /** Undefined variable type **/
    else {
        log_err("Undefined forcing variable type %s in file %i.",
                optstr, file_num);
    }

    param_set.TYPE[type].SUPPLIED = file_num + 1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if (type == SKIP) {
        param_set.TYPE[type].multiplier = 1;
        param_set.TYPE[type].SIGNED = false;
    }
    else {
        sscanf(cmdstr, "%*s %*s %s %lf", flgstr,
               &param_set.TYPE[type].multiplier);
        if (strcasecmp("SIGNED", flgstr) == 0) {
            param_set.TYPE[type].SIGNED = true;
        }
        else {
            param_set.TYPE[type].SIGNED = false;
        }
    }
    param_set.TYPE[type].N_ELEM = 1;

    (*field)++;
}
