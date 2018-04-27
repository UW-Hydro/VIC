/******************************************************************************
 * @section DESCRIPTION
 *
 * Populate model state.
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

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    This function handles tasks related to populating model state.
 *****************************************************************************/
void
vic_populate_model_state(char       *runtype_str,
                         dmy_struct *dmy_current)
{
    extern all_vars_struct *all_vars;
    extern lake_con_struct *lake_con;
    extern domain_struct    local_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern veg_con_struct **veg_con;
    extern filenames_struct filenames;


    size_t                  i;
    unsigned short int      runtype;

    debug("In vic_populate_model_state");

    runtype = start_type_from_char(runtype_str);

    // read the model state from the netcdf file
    if (runtype == CESM_RUNTYPE_RESTART || runtype == CESM_RUNTYPE_BRANCH) {
        // Get restart file from rpointer file
        read_rpointer_file(filenames.init_state.nc_filename);

        // set options.INIT_STATE to true since we have found a state file in
        // the rpointer file.
        options.INIT_STATE = true;

        // read initial state file -- specified in rpointer file
        vic_restore();
    }
    else if (runtype == CESM_RUNTYPE_CLEANSTART) {
        if (options.INIT_STATE) {
            // read initial state file -- specified in global param file
            vic_restore();
        }
        else {
            // no initial state file specified - generate default state
            for (i = 0; i < local_domain.ncells_active; i++) {
                generate_default_state(&(all_vars[i]), &(soil_con[i]),
                                       veg_con[i], dmy_current);
                if (options.LAKES) {
                    generate_default_lake_state(&(all_vars[i]), &(soil_con[i]),
                                                lake_con[i]);
                }
            }
        }
    }

    // compute those state variables that are derived from the others
    for (i = 0; i < local_domain.ncells_active; i++) {
        compute_derived_state_vars(&(all_vars[i]), &(soil_con[i]), veg_con[i]);
        if (options.LAKES) {
            compute_derived_lake_dimensions(&(all_vars[i].lake_var),
                                            lake_con[i]);
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
    char  linestr[MAXSTRING];

    fp = open_file(RPOINTER, "r");

    fgets(linestr, MAXSTRING, fp);

    // Read through rpointer file file to find state file name
    while (!feof(fp)) {
        if (linestr[0] != '#' && linestr[0] != '\n' && linestr[0] != '\0') {
            sscanf(linestr, "%s", fname);
            break;
        }
        fgets(linestr, MAXSTRING, fp);
    }
    fclose(fp);

    // remove trailing newline if present
    fname[strcspn(fname, "\n")] = 0;
}

/******************************************************************************
 * @brief Write rpointer file
 *****************************************************************************/
void
write_rpointer_file(char *fname)
{
    FILE *fp = NULL;
    char *header = "# VIC CESM Driver restart pointer file\n";

    fp = open_file(RPOINTER, "w");

    fprintf(fp, header, fname);
    fprintf(fp, "%s\n", fname);

    fclose(fp);
}
