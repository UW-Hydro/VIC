/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initalizes all global parameters before they are called by
 * the model.
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
 * @brief    Initialize all global parameters before they are called by the
 *           model.
 *****************************************************************************/
void
initialize_global()
{
    extern global_param_struct global_param;

    size_t                     i;

    global_param.dt = MISSING;
    global_param.snow_dt = MISSING;
    global_param.runoff_dt = MISSING;
    global_param.atmos_dt = MISSING;
    global_param.model_steps_per_day = 0;
    global_param.snow_steps_per_day = 0;
    global_param.runoff_steps_per_day = 0;
    global_param.atmos_steps_per_day = 0;
    global_param.nrecs = 0;
    global_param.startyear = 0;
    global_param.startmonth = 0;
    global_param.startday = 0;
    global_param.startsec = 0;
    global_param.endyear = 0;
    global_param.endmonth = 0;
    global_param.endday = 0;
    global_param.resolution = MISSING;
    global_param.wind_h = 10.0;
    for (i = 0; i < 2; i++) {
        global_param.forceyear[i] = 0;
        global_param.forcemonth[i] = 1;
        global_param.forceday[i] = 1;
        global_param.forcesec[i] = 0;
        global_param.forceskip[i] = 0;
        global_param.forceoffset[i] = 0;
    }
    global_param.stateyear = 0;
    global_param.statemonth = 0;
    global_param.stateday = 0;
    global_param.statesec = 0;
    global_param.calendar = CALENDAR_STANDARD;
    global_param.time_units = TIME_UNITS_DAYS;
    global_param.time_origin_num = MISSING;
}
