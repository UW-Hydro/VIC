/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine creates an array of structures that contain information
 * about the day, month and year of each time step.
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
 * @brief    This subroutine creates an array of structures that contain
 *           information about the day, month and year of each time step.
 *****************************************************************************/
dmy_struct *
make_dmy(global_param_struct *global)
{
    extern param_set_struct param_set;

    dmy_struct             *temp;
    dmy_struct              start_dmy, end_dmy, force_dmy;
    size_t                  i;
    unsigned int            offset;
    double                  dt_time_units, start_num, end_num, force_num,
                            numdate;

    start_dmy.dayseconds = global->startsec;
    start_dmy.year = global->startyear;
    start_dmy.day = global->startday;
    start_dmy.month = global->startmonth;

    start_num = date2num(global->time_origin_num, &start_dmy, 0.,
                         global->calendar, global->time_units);

    /** Check if user defined end date instead of number of records **/
    if (global->nrecs == 0) {
        if ((global->endyear == 0) || (global->endmonth == 0) ||
            (global->endday == 0)) {
            log_err("The model global file MUST define EITHER the number of "
                    "records to simulate (NRECS), or the year (ENDYEAR), "
                    "month (ENDMONTH), and day (ENDDAY) of the last full "
                    "simulation day");
        }
        end_dmy.day = global->endday;
        end_dmy.month = global->endmonth;
        end_dmy.year = global->endyear;
        end_dmy.dayseconds = SEC_PER_DAY - global->dt;

        end_num = date2num(global->time_origin_num, &end_dmy, 0.,
                           global->calendar, global->time_units);
        global->nrecs =
            (unsigned int) ((end_num -
                             start_num) * global->model_steps_per_day) + 1;
    }
    else {
        offset =
            (unsigned int) ((double) (SEC_PER_DAY -
                                      start_dmy.dayseconds) / global->dt);
        if ((((unsigned int) global->dt *
              (global->nrecs - offset)) % SEC_PER_DAY) !=
            0) {
            log_err("Nrecs must be defined such that the model ends after "
                    "completing a full day.  Currently Nrecs is set to %zu.",
                    global->nrecs);
        }
    }

    /** Determine number of forcing records to skip before model start time **/
    for (i = 0; i < 2; i++) {
        if (param_set.force_steps_per_day[i] != 0) {
            force_dmy.dayseconds = global->forcesec[i];
            force_dmy.year = global->forceyear[i];
            force_dmy.day = global->forceday[i];
            force_dmy.month = global->forcemonth[i];

            force_num = date2num(global->time_origin_num, &force_dmy, 0.,
                                 global->calendar, global->time_units);

            global->forceskip[i] =
                (unsigned int) round((start_num - force_num) *
                                     (double) param_set.force_steps_per_day[i]);
        }
    }

    // allocate dmy struct
    temp = calloc(global->nrecs, sizeof(*temp));

    /** Create Date Structure for each Model Time Step **/
    for (i = 0; i < global->nrecs; i++) {
        dt_seconds_to_time_units(global->time_units, i * global->dt,
                                 &dt_time_units);
        numdate = start_num + dt_time_units;
        num2date(global->time_origin_num, numdate, 0., global->calendar,
                 global->time_units, &temp[i]);
    }

    return temp;
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
