/******************************************************************************
 * @section DESCRIPTION
 *
 * Function to check whether model state should be saved for the current
 * time step
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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief   Function to check whether model state should be saved for the
 *          current time step
 *****************************************************************************/
bool
check_save_state_flag(dmy_struct *dmy,
                      size_t      current)
{
    extern global_param_struct global_param;

    double                     offset;
    double                     time_num;
    dmy_struct                 dmy_offset;

    // Advance dmy by one timestep because dmy is the "timestep-beginning"
    // timestamp, but we want to check whether the end of the current
    // time step is the user-specified output state time
    offset = global_param.dt / (double) SEC_PER_DAY;
    time_num = date2num(global_param.time_origin_num, &dmy[current], 0,
                        global_param.calendar, TIME_UNITS_DAYS);
    time_num += offset;
    num2date(global_param.time_origin_num, time_num, 0,
             global_param.calendar, TIME_UNITS_DAYS,
             &dmy_offset);

    // Check if the end of the current time step is equal to the state output
    // timestep specified by user
    if (dmy_offset.year == global_param.stateyear &&
        dmy_offset.month == global_param.statemonth &&
        dmy_offset.day == global_param.stateday &&
        dmy_offset.dayseconds == global_param.statesec) {
        return true;
    }
    else {
        return false;
    }
}
