/******************************************************************************
 * @section DESCRIPTION
 *
 * These routine handles the startup tasks for the CESM driver.
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

#include <ctype.h>  // Do we need this line?
#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Wrapper function for VIC startup tasks.
 *****************************************************************************/
void
vic_cesm_start(vic_clock     *vclock,
               case_metadata *cmeta)
{
    extern filep_struct        filep;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern int                 mpi_rank;

    // Initialize structures
    initialize_global_structures();

    // Driver specific settings
    if (mpi_rank == VIC_MPI_ROOT) {
        strcpy(filenames.global, GLOBALPARAM);

        // assign case name to state file name
        strncpy(filenames.statefile, trimstr(cmeta->caseid),
                sizeof(filenames.statefile));

        // read global settings
        filep.globalparam = open_file(filenames.global, "r");
        get_global_param(filep.globalparam);

        // Unpack the vic_clock structure
        // Model timestep
        global_param.dt = (double) vclock->timestep;
        global_param.snow_dt = (double) vclock->timestep;
        global_param.runoff_dt = (double) vclock->timestep;
        global_param.atmos_dt = (double) vclock->timestep;

        global_param.model_steps_per_day =
            (int) ((double) SEC_PER_DAY / global_param.dt);
        global_param.snow_steps_per_day = global_param.model_steps_per_day;
        global_param.runoff_steps_per_day = global_param.model_steps_per_day;
        global_param.atmos_steps_per_day = global_param.model_steps_per_day;

        // Start date/time
        global_param.startyear = vclock->current_year;
        global_param.startmonth = vclock->current_month;
        global_param.startday = vclock->current_day;
        global_param.startsec = vclock->current_dayseconds;
        global_param.nrecs = 1;

        // Calendar
        global_param.calendar = str_to_calendar(trimstr(vclock->calendar));
        // set NR and NF
        NF = global_param.snow_steps_per_day / global_param.model_steps_per_day;
        if (NF == 1) {
            NR = 0;
        }
        else {
            NR = NF;
        }
    }

    // initialize image mode structures and settings
    vic_start();

    // Check that model parameters are valid
    validate_parameters();
    validate_filenames(&filenames);
    validate_global_param(&global_param);
    validate_options(&options);
}

/******************************************************************************
 * @brief    C equivalent of the Fortran TRIM function
 * @note     This function returns a pointer to a substring of the original
 *           string. If the given string was allocated dynamically, the caller
 *           must not overwrite that pointer with the returned value, since the
 *           original pointer must be deallocated using the same allocator with
 *           which it was allocated.  The return value must NOT be deallocated
 *           using free() etc.
 *****************************************************************************/
char *
trimstr(char *str)
{
    char *end;

    // Trim leading space
    while (isspace(*str)) {
        str++;
    }

    if (*str == 0) { // All spaces?
        return str;
    }

    // Trim trailing space
    end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) {
        end--;
    }

    // Write new null terminator
    *(end + 1) = '\0';

    return str;
}
