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

    if ((alarm->freq == FREQ_NEVER) || (alarm->freq == FREQ_NSTEPS) ||
        (alarm->freq == FREQ_DATE) || (alarm->freq == FREQ_END)) {
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

/******************************************************************************
 * @brief   This routine sets an alarm
 *****************************************************************************/
void
set_alarm(dmy_struct   *dmy_current,
          unsigned int  freq,
          void         *value,
          alarm_struct *alarm)
{
    extern global_param_struct global_param;

    alarm->count = 0;
    alarm->next = MISSING;
    alarm->freq = freq;
    alarm->n = MISSING;

    if ((freq == FREQ_NSTEPS) || (freq == FREQ_NSECONDS) ||
        (freq == FREQ_NMINUTES) || (freq == FREQ_NHOURS) ||
        (freq == FREQ_NDAYS) || (freq == FREQ_NMONTHS) ||
        (freq == FREQ_NYEARS)) {
        alarm->n = *((int*) value);
    }
    else if (freq == FREQ_DATE) {
        alarm->date = *((dmy_struct*) value);
    }
    else if ((freq == FREQ_NEVER) || (freq == FREQ_END)) {
        ;  // Do nothing
    }
    else {
        log_err("Did not recognize the frequency value %u", freq);
    }

    // Set alarm->next via reset_alarm
    reset_alarm(alarm, dmy_current);

    // Set subdaily attribute
    if (alarm->next < (int) global_param.model_steps_per_day) {
        alarm->is_subdaily = true;
    }
    else {
        alarm->is_subdaily = false;
    }
}
