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

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    This function handles tasks related to restoring model state.
 *****************************************************************************/
void
vic_restore(void)
{
    extern all_vars_struct    *all_vars;
    extern domain_struct       local_domain;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;

    double                     surf_temp;
    int                        store_offset;
    size_t                     i;
    size_t                     nveg;

    // read first forcing timestep (used in restoring model state)
    // reset the forcing offset to what it was before
    store_offset = global_param.forceoffset[0];
    vic_force();
    global_param.forceoffset[0] = store_offset;

    // read the model state from the netcdf file if there is one
    if (options.INIT_STATE) {
    }

    // run through the remaining VIC initialization routines
    for (i = 0; i < local_domain.ncells_active; i++) {
        // TBD: do something sensible for surf_temp
        surf_temp = 0.;
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_model_state(&(all_vars[i]), nveg, options.Nnode,
                               surf_temp, &(soil_con[i]), veg_con[i]);
    }
}
