/******************************************************************************
 * @section DESCRIPTION
 *
 * This file includes routines that calculate and raise alarms for writing
 * history and state files.
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
 * @brief   This routine resets an alarm
 *****************************************************************************/
void
reset_alarm(alarm_struct *alarm,
            dmy_struct   *dmy_current)
{
    extern global_param_struct global_param;

    alarm->count = 0;
    alarm->next = MISSING;

    if (alarm->freq == FREQ_NEVER) {
        ;  // Do nothing, already set
    }
    else if (alarm->freq == FREQ_NSTEPS) {
        alarm->next = alarm->n;
    }
    else if (alarm->freq == FREQ_DATE) {
        ;  // Do nothing, already set
    }
    else if (alarm->freq == FREQ_END) {
        ;  // Do nothing, already set
    }
    else {
        alarm->next = global_param.model_steps_per_day * time_delta(dmy_current,
                                                                    alarm->freq,
                                                                    alarm->n);
    }
}

/******************************************************************************
 * @brief   This routine raises an alarm
 *****************************************************************************/
bool
raise_alarm(alarm_struct *alarm,
            dmy_struct   *dmy_current)
{
    if ((int) alarm->count == alarm->next) {
        return true;
    }
    else if ((alarm->freq == FREQ_DATE) &&
             (dmy_equal(dmy_current, &(alarm->date)))) {
        return true;
    }
    else {
        return false;
    }
}
