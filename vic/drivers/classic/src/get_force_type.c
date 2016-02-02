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

    char                    optstr[MAXSTRING];
    char                    flgstr[MAXSTRING];
    int                     type;

    type = SKIP;

    /** Initialize flgstr **/
    strcpy(flgstr, "NULL");

    if ((*field) >= (int) param_set.N_TYPES[file_num]) {
        log_err("Too many variables defined for forcing file %i.", file_num);
    }

    sscanf(cmdstr, "%*s %s", optstr);

    /***************************************
       Get meteorological data forcing info
    ***************************************/

    /* type 0: air temperature [C] */
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
    /* type 3: incoming channel flow [m3] */
    else if (strcasecmp("CHANNEL_IN", optstr) == 0) {
        type = CHANNEL_IN;
    }
    /* type 4: direct fraction of shortwave [fraction] */
    else if (strcasecmp("FDIR", optstr) == 0) {
        type = FDIR;
    }
    /* type 5: LAI [m2/m2] */
    else if (strcasecmp("LAI_IN", optstr) == 0) {
        type = LAI_IN;
    }
    /* type 6: incoming longwave radiation [W/m2] */
    else if (strcasecmp("LONGWAVE",
                        optstr) == 0 || strcasecmp("LWDOWN", optstr) == 0) {
        type = LONGWAVE;
    }
    /* type 7: photosynthetically active radiation [uE/m2s] */
    else if (strcasecmp("PAR", optstr) == 0) {
        type = PAR;
    }
    /* type 8: precipitation [mm] */
    else if (strcasecmp("PREC", optstr) == 0) {
        type = PREC;
    }
    /* type 9: air pressure [kPa] */
    else if (strcasecmp("PRESSURE", optstr) == 0) {
        type = PRESSURE;
    }
    /* type 10: vapor pressure [kPa] */
    else if (strcasecmp("VP", optstr) == 0) {
        type = VP;
    }
    /* type 11: rainfall [mm] */
    else if (strcasecmp("SHORTWAVE",
                        optstr) == 0 || strcasecmp("SWDOWN", optstr) == 0) {
        type = SHORTWAVE;
    }
    /* type 12: vegetation cover fraction */
    else if (strcasecmp("VEGCOVER", optstr) == 0) {
        type = VEGCOVER;
    }
    /* type 13: wind speed [m/s] */
    else if (strcasecmp("WIND", optstr) == 0) {
        type = WIND;
    }
    /* type 14: unused (blank) data */
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


/******************************************************************************
 * @brief    This routine determines the counts the number of forcing variables
             in each forcing file specified in the global parameter file.
 *****************************************************************************/
size_t
count_force_vars(FILE *gp)
{
    size_t        nvars;
    unsigned long start_position;
    char          cmdstr[MAXSTRING];
    char          optstr[MAXSTRING];

    // Figure out where we are in the input file
    fflush(gp);
    start_position = ftell(gp);

    // read the first line
    fgets(cmdstr, MAXSTRING, gp);

    // initalize nvars
    nvars = 0;

    // Loop through the lines
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            // line is not blank or a comment
            sscanf(cmdstr, "%s", optstr);

            // if the line starts with FORCE_TYPE
            if (strcasecmp("FORCE_TYPE", optstr) == 0) {
                nvars++;
            }
            // else if we arive at another forcing file break out of loop
            else if (strcasecmp("FORCING1", optstr) == 0 ||
                     strcasecmp("FORCING2", optstr) == 0) {
                break;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    // put the position in the file back to where we started
    fseek(gp, start_position, SEEK_SET);

    return nvars;
}
