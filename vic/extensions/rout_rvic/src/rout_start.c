/******************************************************************************
 * @section DESCRIPTION
 *
 * Read global parameters for routing.
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

#include <rout.h>

/******************************************************************************
 * @brief    Wrapper function for RVIC startup.
 *****************************************************************************/
void
rout_start(void)
{
    extern filenames_struct filenames;
    extern filep_struct     filep;
    extern int              mpi_rank;

        log_info("CHECK-1...\n");

    if (mpi_rank == VIC_MPI_ROOT) {
        // read global settings
        filep.globalparam = open_file(filenames.global, "r");
            log_info("CHECK0...\n");

        get_global_param_rout(filep.globalparam);
    }
}

/******************************************************************************
 * @brief    Read the VIC model global control file, getting values for
 *           global parameters specifically for RVIC.
 *****************************************************************************/
void
get_global_param_rout(FILE *gp)
{
    extern rout_struct rout;
    char               cmdstr[MAXSTRING];
    char               optstr[MAXSTRING];
    
    log_info("CHECK1...\n");

    /** Read through global control file to find parameters **/
    rewind(gp);
    fgets(cmdstr, MAXSTRING, gp);

    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);
    log_info("CHECK2...\n");

            /* Handle case of comment line in which '#' is indented */
            if (optstr[0] == '#') {
                fgets(cmdstr, MAXSTRING, gp);
                continue;
            }

            /*************************************
               Get Model Global Parameters
            *************************************/
            if (strcasecmp("ROUT_PARAM", optstr) == 0) {
                    log_info("CHECK3...\n");

                sscanf(cmdstr, "%*s %s", rout.param_filename);
                log_info("ROUT_PARAM: %s\n", rout.param_filename);
                break;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }
}
