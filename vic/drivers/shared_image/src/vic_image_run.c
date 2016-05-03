/******************************************************************************
 * @section DESCRIPTION
 *
 * Run function for image mode driver.
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
 * @brief    Run VIC for one timestep and store output data
 *****************************************************************************/
void
vic_image_run(dmy_struct *dmy_current)
{
    extern size_t              current;
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern domain_struct       local_domain;
    extern global_param_struct global_param;
    extern lake_con_struct     lake_con;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;

    char                       dmy_str[MAXSTRING];
    size_t                     i;

    // Print the current timestep info before running vic_run
    sprint_dmy(dmy_str, dmy_current);
    debug("Running timestep %zu: %s", current, dmy_str);

    for (i = 0; i < local_domain.ncells_active; i++) {
        // Set global reference string (for debugging inside vic_run)
        sprintf(vic_run_ref_str, "Gridcell io_idx: %zu, timestep info: %s",
                local_domain.locations[i].io_idx, dmy_str);

        update_step_vars(&(all_vars[i]), veg_con[i], veg_hist[i]);
        vic_run(&(atmos[i]), &(all_vars[i]), dmy_current, &global_param,
                &lake_con, &(soil_con[i]), veg_con[i], veg_lib[i]);
        put_data(&(all_vars[i]), &(atmos[i]), &(soil_con[i]), veg_con[i],
                 veg_lib[i], &lake_con, out_data[i], &(save_data[i]),
                 current);
    }
}
