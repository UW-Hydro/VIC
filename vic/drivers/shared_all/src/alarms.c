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
            dmy_struct   *dmy_current)
{
    extern global_param_struct global_param;

    double                     delta;
    double                     current;
    double                     next;
    double                     offset;
    dmy_struct                 dmy_current_offset;

    alarm->count = 0;

    if ((alarm->freq == FREQ_NEVER) || (alarm->freq == FREQ_NSTEPS) ||
        (alarm->freq == FREQ_DATE) || (alarm->freq == FREQ_END)) {
        ;  // Do nothing, already set
    }
    else if (alarm->freq == FREQ_NMONTHS) {
        // If aggregation frequency is NMONTHS, then first shift dmy_current
        // forward for one timestep, advance month(s), then shift backward for
        // one timestep. This is to avoid the following issue: the most common
        // usage of AGGFREQ = NMONTHS is to start a simulation from 00:00:00
        // of the first day of a certain month; since dmy_current is
        // "timestep-beginning", the aggregation window timestamp will be the
        // beginning of the last timestep of a month. This will cause problem
        // because of different number of days in each month. This shift here
        // will avoid this problem.
        // NOTE: if a simulation does not start from the beginning of a month,
        // there might be a problem if startday > 28 !!!

        // Shift forward by one time step
        offset = global_param.dt / (double) SEC_PER_DAY;
        current = date2num(global_param.time_origin_num, dmy_current, 0,
                           global_param.calendar, TIME_UNITS_DAYS);
        current += offset;
        num2date(global_param.time_origin_num, current, 0,
                 global_param.calendar, TIME_UNITS_DAYS,
                 &dmy_current_offset);
        // Advance
        delta = time_delta(&dmy_current_offset, alarm->freq, alarm->n);
        current = date2num(global_param.time_origin_num, &dmy_current_offset,
                           0, global_param.calendar, TIME_UNITS_DAYS);
        next = delta + current;
        // Shift backward by one time step
        next -= offset;
        num2date(global_param.time_origin_num, next, 0,
                 global_param.calendar, TIME_UNITS_DAYS,
                 &(alarm->next_dmy));
    }
    else {
        // If other frequency types, directly advance without shifting
        delta = time_delta(dmy_current, alarm->freq, alarm->n);
        current = date2num(global_param.time_origin_num, dmy_current, 0,
                           global_param.calendar, TIME_UNITS_DAYS);
        next = delta + current;
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
          alarm_struct *alarm)
{
    extern global_param_struct global_param;
    dmy_struct                 dmy_current_offset;
    double                     delta;
    double                     current;

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
    // set_alarm only gets called in the initialization step. Therefore,
    // shift dmy_current to one timestep earlier before advancing for one
    // aggregation window, because dmy_current does not change after running
    // the first timestep (so it would still be the beginning of the first
    // timestep)
    if ((alarm->freq == FREQ_NEVER) || (alarm->freq == FREQ_NSTEPS) ||
        (alarm->freq == FREQ_DATE) || (alarm->freq == FREQ_END)) {
        ;  // Do nothing, already set
    }
    else {
        delta = time_delta(dmy_current, FREQ_NSECONDS,
                           (int) global_param.dt);
        current = date2num(global_param.time_origin_num, dmy_current, 0,
                           global_param.calendar, TIME_UNITS_DAYS);
        current -= delta;
        num2date(global_param.time_origin_num, current, 0,
                 global_param.calendar, TIME_UNITS_DAYS,
                 &dmy_current_offset);
    }
    // set alarm->next
    reset_alarm(alarm, &dmy_current_offset);

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
