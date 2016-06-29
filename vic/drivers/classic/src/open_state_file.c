/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine opens the model state file for output.
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
 * @brief    Open state file to write to.
 *****************************************************************************/
FILE *
open_state_file(global_param_struct *global,
                filenames_struct     filenames,
                size_t               Nlayer,
                size_t               Nnodes)
{
    extern option_struct options;

    FILE                *statefile;
    char                 filename[MAXSTRING];

    /* open state file */
    sprintf(filename, "%s", filenames.statefile);
    if (options.STATE_FORMAT == BINARY) {
        statefile = open_file(filename, "wb");
    }
    else {
        statefile = open_file(filename, "w");
    }

    /* Write save state date information */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&global->stateyear, sizeof(int), 1, statefile);
        fwrite(&global->statemonth, sizeof(int), 1, statefile);
        fwrite(&global->stateday, sizeof(int), 1, statefile);
    }
    else {
        fprintf(statefile, "%i %i %i\n", global->stateyear,
                global->statemonth, global->stateday);
    }

    /* Write simulation flags */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&Nlayer, sizeof(size_t), 1, statefile);
        fwrite(&Nnodes, sizeof(size_t), 1, statefile);
    }
    else {
        fprintf(statefile, "%zu %zu\n", Nlayer, Nnodes);
    }

    return(statefile);
}
