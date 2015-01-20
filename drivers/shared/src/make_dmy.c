/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine creates an array of structures that contain information
 * about the day, month and year of each time step.
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_shared.h>

#ifndef _LEAPYR
#define LEAPYR(y) (!((y) % 400) || (!((y) % 4) && ((y) % 100)))
#endif

/******************************************************************************
 * @brief    This subroutine creates an array of structures that contain
 *           information about the day, month and year of each time step.
 *****************************************************************************/
dmy_struct *
make_dmy(global_param_struct *global)
{
    extern param_set_struct param_set;

    dmy_struct             *temp;
    unsigned short          hr, year, day, month, jday, ii;
    unsigned short          days[MONTHS_PER_YEAR] = {
        31, 28, 31, 30, 31, 30,
        31, 31, 30, 31, 30, 31
    };
    unsigned short          endmonth, endday, endyear, skiprec, i, offset;
    unsigned short          tmpmonth, tmpday, tmpyear, tmphr, tmpjday;
    char                    DONE;

    hr = global->starthour;
    year = global->startyear;
    day = global->startday;
    month = global->startmonth;

    /** Check if user defined end date instead of number of records **/
    if (global->nrecs == 0) {
        if ((global->endyear == 0) || (global->endmonth == 0) ||
            (global->endday == 0)) {
            log_err("The model global file MUST define EITHER the number of "
                    "records to simulate (NRECS), or the year (ENDYEAR), "
                    "month (ENDMONTH), and day (ENDDAY) of the last full "
                    "simulation day");
        }
        endday = global->endday;
        endmonth = global->endmonth;
        endyear = global->endyear;
        if (LEAPYR(endyear)) {
            days[1] = 29;
        }
        else {
            days[1] = 28;
        }
        if (endday < days[global->endmonth - 1]) {
            endday++;
        }
        else {
            endday = 1;
            endmonth++;
            if (endmonth > MONTHS_PER_YEAR) {
                endmonth = 1;
                endyear++;
            }
        }

        DONE = false;
        ii = 0;

        tmpyear = year;
        tmpmonth = month;
        tmpday = day;
        tmphr = hr;
        while (!DONE) {
            get_next_time_step(&tmpyear, &tmpmonth, &tmpday, &tmphr,
                               &tmpjday, global->dt);
            ii++;
            if (tmpyear == endyear) {
                if (tmpmonth == endmonth) {
                    if (tmpday == endday) {
                        DONE = true;
                    }
                }
            }
        }
        global->nrecs = ii;
    }
    else {
        offset = 0;
        tmphr = hr;
        while (tmphr != 0) {
            tmphr += global->dt;
            offset++;
            if (tmphr >= HOURS_PER_DAY) {
                tmphr = 0;
            }
        }
        if (((global->dt * (global->nrecs - offset)) % HOURS_PER_DAY) != 0) {
            log_err("Nrecs must be defined such that the model ends after "
                    "completing a full day.  Currently Nrecs is set to %i, "
                    "while %i and %i are allowable values.", global->nrecs,
                    ((global->dt * (global->nrecs - offset)) / HOURS_PER_DAY) * HOURS_PER_DAY,
                    ((global->dt *
                      (global->nrecs -
                       offset)) /
                     HOURS_PER_DAY) * HOURS_PER_DAY + HOURS_PER_DAY);
        }
    }

    // allocate dmy struct
    temp = (dmy_struct*) calloc(global->nrecs, sizeof(dmy_struct));

    /** Create Date Structure for each Modeled Time Step **/
    jday = day;
    if (LEAPYR(year)) {
        days[1] = 29;
    }
    else {
        days[1] = 28;
    }
    for (ii = 0; ii < month - 1; ii++) {
        jday += days[ii];
    }

    DONE = false;
    ii = 0;

    while (!DONE) {
        temp[ii].hour = hr;
        temp[ii].day = day;
        temp[ii].month = month;
        temp[ii].year = year;
        temp[ii].day_in_year = jday;

        get_next_time_step(&year, &month, &day, &hr, &jday, global->dt);

        ii++;
        if (ii == global->nrecs) {
            DONE = true;
        }
    }

    /** Determine number of forcing records to skip before model start time **/
    for (i = 0; i < 2; i++) {
        if (param_set.FORCE_DT[i] != MISSING) {
            if (global->forceyear[i] > 0) {
                tmpyear = global->forceyear[i];
                tmpmonth = global->forcemonth[i];
                tmpday = global->forceday[i];
                tmphr = global->forcehour[i];
                tmpjday = tmpday;
                if (LEAPYR(tmpyear)) {
                    days[1] = 29;
                }
                else {
                    days[1] = 28;
                }
                for (ii = 0; ii < tmpmonth - 1; ii++) {
                    tmpjday += days[ii];
                }

                while (tmpyear < temp[0].year ||
                       (tmpyear == temp[0].year && tmpjday <
                        temp[0].day_in_year) ||
                       (tmpyear == temp[0].year && tmpjday ==
                        temp[0].day_in_year && tmphr < temp[0].hour)) {
                    get_next_time_step(&tmpyear, &tmpmonth, &tmpday, &tmphr,
                                       &tmpjday, global->dt);

                    global->forceskip[i]++;
                }
                global->forceoffset[i] = global->forceskip[i];
            }
        }
    }

    /** Determine the number of records to skip before starting output files **/
    skiprec = 0;
    for (i = 0; i < global->skipyear; i++) {
        if (LEAPYR(temp[skiprec].year)) {
            skiprec += DAYS_PER_LYEAR * HOURS_PER_DAY / global->dt;
        }
        else {
            skiprec += DAYS_PER_YEAR * HOURS_PER_DAY / global->dt;
        }
    }
    global->skipyear = skiprec;

    return temp;
}

/******************************************************************************
 * @brief    Get the next timestep.
 *****************************************************************************/
void
get_next_time_step(unsigned short *year,
                   unsigned short *month,
                   unsigned short *day,
                   unsigned short *hr,
                   unsigned short *jday,
                   unsigned short  dt)
{
    unsigned short days[MONTHS_PER_YEAR] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };

    *hr += dt;
    if (*hr >= HOURS_PER_DAY) {
        *hr = 0;
        *day += 1;
        *jday += 1;

        if (LEAPYR(*year)) {
            days[1] = 29;
        }
        else {
            days[1] = 28;
        }

        if (*day > days[*month - 1]) {
            *day = 1;
            *month += 1;
            if (*month == 13) {
                *month = 1;
                *jday = 1;
                *year += 1;
            }
        }
    }
}

/******************************************************************************
 * @brief    Free the dmy array.
 *****************************************************************************/
void
free_dmy(dmy_struct **dmy)
{
    if (*dmy == NULL) {
        return;
    }

    free(*dmy);
}
