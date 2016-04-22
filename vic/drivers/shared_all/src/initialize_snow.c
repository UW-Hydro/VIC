/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the snow variable arrays for each new grid cell.
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
    double               dummy_double;
    unsigned int         dummy_unsigned_int;
    bool                 dummy_bool;

    for (i = 0; i <= veg_num; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Prognostic states
            snow[i][j].albedo = dummy_double;
            snow[i][j].canopy_albedo = dummy_double;
            snow[i][j].coldcontent = dummy_double;
            snow[i][j].coverage = dummy_double;
            snow[i][j].density = dummy_double;
            snow[i][j].depth = dummy_double;
            snow[i][j].last_snow = dummy_unsigned_int;
            snow[i][j].max_snow_depth = dummy_double;
            snow[i][j].MELTING = dummy_bool;
            snow[i][j].pack_temp = dummy_double;
            snow[i][j].pack_water = dummy_double;
            snow[i][j].snow = dummy_bool;
            snow[i][j].snow_canopy = dummy_double;
            snow[i][j].store_coverage = dummy_double;
            snow[i][j].store_snow = dummy_bool;
            snow[i][j].store_swq = dummy_double;
            snow[i][j].surf_temp = dummy_double;
            snow[i][j].surf_temp_fbcount = dummy_unsigned_int;
            snow[i][j].surf_temp_fbflag = dummy_bool;
            snow[i][j].surf_water = dummy_double;
            snow[i][j].swq = dummy_double;
            snow[i][j].snow_distrib_slope = dummy_double;
            snow[i][j].tmp_int_storage = dummy_double;
            // Fluxes
            snow[i][j].blowing_flux = dummy_double;
            snow[i][j].canopy_vapor_flux = dummy_double;
            snow[i][j].mass_error = dummy_double;
            snow[i][j].melt = dummy_double;
            snow[i][j].Qnet = dummy_double;
            snow[i][j].surface_flux = dummy_double;
            snow[i][j].transport = dummy_double;
            snow[i][j].vapor_flux = dummy_double;
        }
    }
}
