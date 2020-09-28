/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the lake model state (energy balance, water balance,
 * snow components) to default values.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize the lake model state (energy balance, water balance,
 *           and snow components) to default values.
 *****************************************************************************/
void
generate_default_lake_state(lake_var_struct *lake,
                            soil_con_struct *soil_con,
                            lake_con_struct  lake_con)
{
    extern option_struct options;

    size_t               k;

    /************************************************************************
       Initialize lake state variables
       TBD: currently setting depth to depth_in from parameter file, but
            in future we should initialize to mindepth as default and
            eliminate depth_in (require user to use a state file if they
            want control over initial depth)
    ************************************************************************/
    if (options.LAKES) {
        lake->ldepth = lake_con.depth_in;
        for (k = 0; k < lake->activenod; k++) {
            // lake model requires FULL_ENERGY set to true
            lake->temp[k] = soil_con->avg_temp;
        }
    }
}
