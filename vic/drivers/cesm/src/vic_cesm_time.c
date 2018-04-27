/******************************************************************************
 * @section DESCRIPTION
 *
 * CESM driver time handling
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

#include <vic_driver_cesm.h>

double dt_time_units = MISSING;
double numdate = MISSING;


/******************************************************************************
 * @brief    Initialize cesm time
 *****************************************************************************/
void
initialize_cesm_time(void)
{
    extern size_t              current;
    extern dmy_struct          dmy_current;
    extern global_param_struct global_param;

    debug("In initialize_cesm_time");

    current = 0;

    // initialize time
    initialize_time();

    // Set dmy_current using global param start date/time
    dmy_current.year = global_param.startyear;
    dmy_current.month = global_param.startmonth;
    dmy_current.day = global_param.startday;
    dmy_current.dayseconds = global_param.startsec;

    // initialize module level globals
    dt_seconds_to_time_units(global_param.time_units, global_param.dt,
                             &dt_time_units);

    // initialize numdate
    numdate = date2num(global_param.time_origin_num, &dmy_current, 0.,
                       global_param.calendar, global_param.time_units);

    num2date(global_param.time_origin_num, numdate, 0., global_param.calendar,
             global_param.time_units, &dmy_current);
}

/******************************************************************************
 * @brief    Finalize cesm time
 *****************************************************************************/
void
finalize_cesm_time(vic_clock *vclock)
{
    extern size_t              current;
    extern global_param_struct global_param;

    // populate fields in global_param needed for timing tables
    global_param.nrecs = current;
    global_param.endyear = vclock->current_year;
    global_param.endmonth = vclock->current_month;
    global_param.endday = vclock->current_day;
}

/******************************************************************************
 * @brief    Advance one timestep
 *****************************************************************************/
void
advance_vic_time(void)
{
    extern size_t              current;
    extern dmy_struct          dmy_current;
    extern global_param_struct global_param;

    char                       dmy_string[MAXSTRING];

    current++;

    numdate += dt_time_units;

    num2date(global_param.time_origin_num, numdate, 0., global_param.calendar,
             global_param.time_units, &dmy_current);

    if (invalid_date(global_param.calendar, &dmy_current)) {
        sprint_dmy(dmy_string, &dmy_current);
        log_err("Invalid date encountered while advancing VIC clock in "
                "timestep %zu.\n%s", current, dmy_string);
    }
}

/******************************************************************************
 * @brief    Raise an error if vclock and dmy don't match
 *****************************************************************************/
void
assert_time_insync(vic_clock  *vclock,
                   dmy_struct *dmy)
{
    if ((vclock->current_year != dmy->year) ||
        (vclock->current_month != dmy->month) ||
        (vclock->current_day != dmy->day) ||
        (vclock->current_dayseconds != dmy->dayseconds)) {
        log_err("Clocks out of sync");
    }
}
