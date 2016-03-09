/******************************************************************************
 * @section DESCRIPTION
 *
 * Restore model state.
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
#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    This function handles tasks related to restoring model state.
 *****************************************************************************/
void
vic_restore(char *runtype_str)
{
    extern all_vars_struct *all_vars;
    extern domain_struct    local_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern veg_con_struct **veg_con;
    extern filenames_struct filenames;


    double                  surf_temp;
    size_t                  i;
    size_t                  nveg;
    unsigned short int      runtype;

    debug("In vic_restore");

    runtype = start_type_from_char(runtype_str);

    // read first forcing timestep (used in restoring model state)
    // reset the forcing offset to what it was before
    vic_force();

    // read the model state from the netcdf file if there is one
    if (runtype == CESM_RUNTYPE_RESTART || runtype == CESM_RUNTYPE_BRANCH) {
        // Get restart file from rpointer file
        read_rpointer_file(filenames.init_state);

        // TODO: read initial state file
    }
    else if (runtype == CESM_RUNTYPE_CLEANSTART) {
        // run type is clean start
        for (i = 0; i < local_domain.ncells_active; i++) {
            // TBD: do something sensible for surf_temp
            surf_temp = 0.;
            nveg = veg_con[i][0].vegetat_type_num;
            initialize_model_state(&(all_vars[i]), nveg, options.Nnode,
                                   surf_temp, &(soil_con[i]), veg_con[i]);
        }
    }
}

/******************************************************************************
 * @brief Convert runtype string to enum integer
 * @note  See CESM's seq_infodata_mod.F90 for more information
 *****************************************************************************/
unsigned short int
start_type_from_char(char *start_str)
{
    if (strcasecmp("startup", start_str) == 0) {
        return CESM_RUNTYPE_CLEANSTART;
    }
    else if (strcasecmp("continue", start_str) == 0) {
        return CESM_RUNTYPE_RESTART;
    }
    else if (strcasecmp("branch", start_str) == 0) {
        return CESM_RUNTYPE_BRANCH;
    }
    else {
        log_err("Unknown calendar specified: %s", start_str);
    }
}

/******************************************************************************
 * @brief Read rpointer file
 *****************************************************************************/
void
read_rpointer_file(char *fname)
{
    FILE *fp = NULL;

    fname = NULL;

    fp = fopen(RPOINTER, "r");

    if (fp == NULL) {
        log_err("Error reading rpointer file %s", RPOINTER);
    }

    fgets(fname, MAXSTRING, fp);

    fclose(fp);

    if (fname == NULL) {
        log_err("Error reading rpointer file %s", RPOINTER);
    }
}
