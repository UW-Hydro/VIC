/******************************************************************************
 * @section DESCRIPTION
 *
 * This program builds the files names for input and output of grided data
 * files.
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
 * @brief    Build files names for input and output of grided data files.
 *****************************************************************************/
void
make_in_and_outfiles(filep_struct     *filep,
                     filenames_struct *filenames,
                     soil_con_struct  *soil,
                     stream_struct   **streams,
                     dmy_struct       *dmy)
{
    extern option_struct    options;
    extern param_set_struct param_set;
    extern FILE *open_file(char string[], char type[]);

    char                    latchar[20], lngchar[20], junk[6];
    size_t                  filenum;

    sprintf(junk, "%%.%if", options.GRID_DECIMAL);
    sprintf(latchar, junk, soil->lat);
    sprintf(lngchar, junk, soil->lng);

    /********************************
       Input Forcing Files
    ********************************/

    strcpy(filenames->forcing[0], filenames->f_path_pfx[0]);
    strcat(filenames->forcing[0], latchar);
    strcat(filenames->forcing[0], "_");
    strcat(filenames->forcing[0], lngchar);
    if (param_set.FORCE_FORMAT[0] == BINARY) {
        filep->forcing[0] = open_file(filenames->forcing[0], "rb");
    }
    else {
        filep->forcing[0] = open_file(filenames->forcing[0], "r");
    }

    filep->forcing[1] = NULL;
    if (strcasecmp(filenames->f_path_pfx[1], "MISSING") != 0) {
        strcpy(filenames->forcing[1], filenames->f_path_pfx[1]);
        strcat(filenames->forcing[1], latchar);
        strcat(filenames->forcing[1], "_");
        strcat(filenames->forcing[1], lngchar);
        if (param_set.FORCE_FORMAT[0] == BINARY) {
            filep->forcing[1] = open_file(filenames->forcing[1], "rb");
        }
        else {
            filep->forcing[1] = open_file(filenames->forcing[1], "r");
        }
    }

    /********************************
       Output Files
    ********************************/

    for (filenum = 0; filenum < options.Noutstreams; filenum++) {
        strcpy((*streams)[filenum].filename, filenames->result_dir);
        strcat((*streams)[filenum].filename, "/");
        strcat((*streams)[filenum].filename,
               (*streams)[filenum].prefix);
        strcat((*streams)[filenum].filename, "_");
        strcat((*streams)[filenum].filename, latchar);
        strcat((*streams)[filenum].filename, "_");
        strcat((*streams)[filenum].filename, lngchar);
        if ((*streams)[filenum].file_format == BINARY) {
            strcat((*streams)[filenum].filename, ".bin");
            (*streams)[filenum].fh = open_file(
                (*streams)[filenum].filename, "wb");
        }
        else if ((*streams)[filenum].file_format == ASCII) {
            strcat((*streams)[filenum].filename, ".txt");
            (*streams)[filenum].fh = open_file(
                (*streams)[filenum].filename, "w");
        }
        else {
            log_err("Unrecognized OUT_FORMAT option");
        }
    }
    /** Write output file headers **/
    write_header(streams, dmy);
}
