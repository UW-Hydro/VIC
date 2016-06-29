/******************************************************************************
 * @section DESCRIPTION
 *
 * Determines from the air temperature what fraction of incoming precipitation
 * is frozen and unfrozen (snow and rain).
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

#include <vic_run.h>

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
double
calc_rainonly(double air_temp,
              double prec,
              double MAX_SNOW_TEMP,
              double MIN_RAIN_TEMP)
{
    double rainonly;

    rainonly = 0.;
    if (MAX_SNOW_TEMP <= MIN_RAIN_TEMP) {
        log_err("MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
    }
    if (air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
        rainonly = (air_temp - MIN_RAIN_TEMP) /
                   (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
    }
    else if (air_temp >= MAX_SNOW_TEMP) {
        rainonly = prec;
    }

    return(rainonly);
}
