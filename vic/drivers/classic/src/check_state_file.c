/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine opens a model state file and verifys that the starting date,
 * number of layers and number of thermal nodes in the file agrees with what
 * was defined in the model global control file.
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
 * @brief    This subroutine opens a model state file and verifys that the
             starting date, number of layers and number of thermal nodes in the
             file agrees with what was defined in the model global control file.
 *****************************************************************************/
FILE *
check_state_file(char  *init_state_name,
                 size_t Nlayer,
                 size_t Nnodes,
                 int   *startrec)
{
    extern option_struct options;

    FILE                *init_state;
    size_t               tmp_Nlayer;
    size_t               tmp_Nnodes;
    unsigned short int   startday, startmonth, startyear;

    /* open state file */
    if (options.STATE_FORMAT == BINARY) {
        init_state = open_file(init_state_name, "rb");
    }
    else {
        init_state = open_file(init_state_name, "r");
    }

    /* Initialize startrec */
    *startrec = 0;

    /* Check state date information */
    if (options.STATE_FORMAT == BINARY) {
        fread(&startyear, sizeof(int), 1, init_state);
        fread(&startmonth, sizeof(int), 1, init_state);
        fread(&startday, sizeof(int), 1, init_state);
    }
    else {
        fscanf(init_state, "%hu %hu %hu\n", &startyear, &startmonth, &startday);
    }

    /* Check simulation options */
    if (options.STATE_FORMAT == BINARY) {
        fread(&tmp_Nlayer, sizeof(size_t), 1, init_state);
        fread(&tmp_Nnodes, sizeof(size_t), 1, init_state);
    }
    else {
        fscanf(init_state, "%zu %zu\n", &tmp_Nlayer, &tmp_Nnodes);
    }
    if (tmp_Nlayer != Nlayer) {
        log_err("The number of soil moisture layers in the model state file "
                "(%zu) does not equal that defined in the global control file "
                "(%zu).  Check your input files.", tmp_Nlayer, Nlayer);
    }
    if (tmp_Nnodes != Nnodes) {
        log_err("The number of soil thermal nodes in the model state file "
                "(%zu) does not equal that defined in the global control file "
                "(%zu).  Check your input files.", tmp_Nnodes, Nnodes);
    }

    return(init_state);
}
