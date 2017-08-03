/******************************************************************************
 * @section DESCRIPTION
 *
 * This function populates the model state.
 *
 * If a state file is provided to the model then its contents are checked
 * to see if it agrees with the current simulation set-up, if so it is used
 * to initialize the model state.  If no state file is provided then the
 * model initializes all variables with defaults and the user should expect
 * to throw out the beginning of the simulation period as model spin-up.
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

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    This function handles tasks related to populating model state.
 *****************************************************************************/
void
vic_populate_model_state(dmy_struct *dmy_current)
{
    extern all_vars_struct *all_vars;
    extern lake_con_struct *lake_con;
    extern domain_struct    local_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern veg_con_struct **veg_con;

    size_t                  i;

    // read the model state from the netcdf file if there is one
    if (options.INIT_STATE) {
        vic_restore();
    }
    else {
        // else generate a default state
        for (i = 0; i < local_domain.ncells_active; i++) {
            generate_default_state(&(all_vars[i]), &(soil_con[i]), veg_con[i],
                                   dmy_current);
            if (options.LAKES) {
                generate_default_lake_state(&(all_vars[i]), &(soil_con[i]),
                                            lake_con[i]);
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
