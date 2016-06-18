/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate distance between two locations.
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
 * @brief    Get distance between two locations.
 *****************************************************************************/
double
get_dist(double lat1,
         double long1,
         double lat2,
         double long2)
{
    double theta1;
    double phi1;
    double theta2;
    double phi2;
    double dtor;
    double term1;
    double term2;
    double term3;
    double temp;
    double distance;

    dtor = 2.0 * CONST_PI / 360.0;
    theta1 = dtor * long1;
    phi1 = dtor * lat1;
    theta2 = dtor * long2;
    phi2 = dtor * lat2;
    term1 = cos(phi1) * cos(theta1) * cos(phi2) * cos(theta2);
    term2 = cos(phi1) * sin(theta1) * cos(phi2) * sin(theta2);
    term3 = sin(phi1) * sin(phi2);
    temp = term1 + term2 + term3;
    temp = (double) (1.0 < temp) ? 1.0 : temp;
    distance = CONST_REARTH * acos(temp);

    return distance;
}
