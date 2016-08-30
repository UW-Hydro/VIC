/******************************************************************************
 * @section DESCRIPTION
 *
 * This file includes routines that calculate and raise alarms for writing
 * history and state files.
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
 * @brief   This routine resets an alarm
 *****************************************************************************/
void
reset_alarm(alarm_struct *alarm,
            dmy_struct   *dmy_current,
            double        offset)
{
    extern global_param_struct global_param;

    double                     delta;
    double                     current;
    double                     next;

    alarm->count = 0;

    if ((alarm->freq == FREQ_NEVER) || (alarm->freq == FREQ_NSTEPS) ||
        (alarm->freq == FREQ_DATE) || (alarm->freq == FREQ_END)) {
        ;  // Do nothing, already set
    }
    else {
        delta = time_delta(dmy_current, alarm->freq, alarm->n);
        current = date2num(global_param.time_origin_num, dmy_current, 0,
                           global_param.calendar, TIME_UNITS_DAYS);
        next = delta + current + offset;
        num2date(global_param.time_origin_num, next, 0,
                 global_param.calendar, TIME_UNITS_DAYS,
                 &(alarm->next_dmy));
    }
}

/******************************************************************************
 * @brief   This routine raises an alarm
 *****************************************************************************/
bool
raise_alarm(alarm_struct *alarm,
            dmy_struct   *dmy_current)
{
    if (((alarm->freq == FREQ_NSTEPS) &&
         (alarm->next_count == (int) alarm->count)) ||
        (dmy_equal(dmy_current, &(alarm->next_dmy)))) {
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
          alarm_struct *alarm,
          double        offset)
{
    extern global_param_struct global_param;

    alarm->count = 0;
    alarm->freq = freq;
    alarm->n = MISSING;
    alarm->next_count = MISSING;

    if (freq == FREQ_NSTEPS) {
        alarm->n = *((int*) value);
        alarm->next_count = alarm->n;
        if (alarm->n <= 0) {
            log_err("invalid n (%d) provided to set_alarm", alarm->n);
        }
    }
    else if ((freq == FREQ_NSECONDS) || (freq == FREQ_NMINUTES) ||
             (freq == FREQ_NHOURS) || (freq == FREQ_NDAYS) ||
             (freq == FREQ_NMONTHS) || (freq == FREQ_NYEARS)) {
        alarm->n = *((int*) value);
        if (alarm->n <= 0) {
            log_err("invalid n (%d) provided to set_alarm", alarm->n);
        }
    }
    else if (freq == FREQ_DATE) {
        alarm->next_dmy = *((dmy_struct*) value);
    }
    else if ((freq == FREQ_NEVER) || (freq == FREQ_END)) {
        ;  // Do nothing
    }
    else {
        log_err("Did not recognize the frequency value %u", freq);
    }

    // Set alarm->next via reset_alarm
    reset_alarm(alarm, dmy_current, offset);

    // Set subdaily attribute
    if (((freq == FREQ_NSTEPS) &&
         (alarm->next_count < (int) global_param.model_steps_per_day)) ||
        ((freq == FREQ_NSECONDS) && (alarm->n < SEC_PER_DAY)) ||
        ((freq == FREQ_NMINUTES) && (alarm->n < MIN_PER_DAY)) ||
        ((freq == FREQ_NHOURS) && (alarm->n < HOURS_PER_DAY))) {
        alarm->is_subdaily = true;
    }
    else {
        alarm->is_subdaily = false;
    }
}
