/******************************************************************************
* @section DESCRIPTION
*
* This subroutine computes grid cell area.
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
******************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
* @brief        This subroutine computes grid cell area.
******************************************************************************/
void
compute_cell_area(soil_con_struct *soil_con)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    int                        i;
    double                     lat;
    double                     lon;
    double                     start_lat;
    double                     right_lon;
    double                     left_lon;
    double                     delta;
    double                     area;

    if (options.EQUAL_AREA) {
        soil_con->cell_area = global_param.resolution * M_PER_KM * M_PER_KM; /* Grid cell area in m^2. */
    }
    else {
        lat = fabs(soil_con->lat);
        lon = fabs(soil_con->lng);

        start_lat = lat - global_param.resolution / 2;
        right_lon = lon + global_param.resolution / 2;
        left_lon = lon - global_param.resolution / 2;

        delta = get_dist(lat, lon, lat + global_param.resolution / 10., lon);

        area = 0.;

        for (i = 0; i < 10; i++) {
            area += get_dist(start_lat, left_lon, start_lat, right_lon) * delta;
            start_lat += global_param.resolution / 10;
        }

        soil_con->cell_area = area; /* Grid cell area in m^2. */
    }
}
