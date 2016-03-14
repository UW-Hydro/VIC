/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting information
 * for output variables list (if any).
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Get output info from global parameter file.
 *****************************************************************************/
int
parse_output_info(FILE             *gp,
                  out_data_struct **out_data)
{
    char   cmdstr[MAXSTRING];
    char   optstr[MAXSTRING];
    char   varname[MAXSTRING];
    bool   found;
    int    outvarnum;
    size_t i;

    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);
    outvarnum = 0;
    while (!feof(gp)) {
        found = false;
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);
            if (strcasecmp("OUTVAR", optstr) == 0) {
                sscanf(cmdstr, "%*s %s", varname);
                for (i = 0; i < N_OUTVAR_TYPES; i++) {
                    if (strcmp(out_data[0][i].varname, varname) == 0) {
                        found = true;
                        out_data[0][i].write = true;
                    }
                }
                if (!found) {
                    log_err("\"%s\" was not found in the list of supported "
                            "output variable names.  Please use "
                            "the exact name listed in vic_def.h.",
                            varname);
                }
                outvarnum++;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    return outvarnum;
}
