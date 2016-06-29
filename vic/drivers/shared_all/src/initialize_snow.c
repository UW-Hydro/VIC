/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the snow variable arrays for each new grid cell.
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
 * @brief    Initialize the snow variable arrays for each new grid cell.
 *****************************************************************************/
void
initialize_snow(snow_data_struct **snow,
                size_t             veg_num)
{
    extern option_struct options;
    size_t               i, j;

    for (i = 0; i <= veg_num; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Prognostic states
            snow[i][j].albedo = 0.0;
            snow[i][j].canopy_albedo = 0.0;
            snow[i][j].coldcontent = 0.0;
            snow[i][j].coverage = 0.0;
            snow[i][j].density = 0.0;
            snow[i][j].depth = 0.0;
            snow[i][j].last_snow = 0;
            snow[i][j].max_snow_depth = 0.0;
            snow[i][j].MELTING = false;
            snow[i][j].pack_temp = 0.0;
            snow[i][j].pack_water = 0.0;
            snow[i][j].snow = false;
            snow[i][j].snow_canopy = 0.0;
            snow[i][j].store_coverage = 0.0;
            snow[i][j].store_snow = false;
            snow[i][j].store_swq = 0.0;
            snow[i][j].surf_temp = 0.0;
            snow[i][j].surf_temp_fbcount = 0;
            snow[i][j].surf_temp_fbflag = false;
            snow[i][j].surf_water = 0.0;
            snow[i][j].swq = 0.0;
            snow[i][j].snow_distrib_slope = 0.0;
            snow[i][j].tmp_int_storage = 0.0;
            // Fluxes
            snow[i][j].blowing_flux = 0.0;
            snow[i][j].canopy_vapor_flux = 0.0;
            snow[i][j].mass_error = 0.0;
            snow[i][j].melt = 0.0;
            snow[i][j].Qnet = 0.0;
            snow[i][j].surface_flux = 0.0;
            snow[i][j].transport = 0.0;
            snow[i][j].vapor_flux = 0.0;
        }
    }
}
