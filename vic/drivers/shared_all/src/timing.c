/******************************************************************************
 * @section DESCRIPTION
 *
 * Routines to calculate and store model runtime timing.
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
 * @brief    Get wall time
 *****************************************************************************/
double
get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        log_err("Unable to get time of day")
    }
    return (double) time.tv_sec + (double) time.tv_usec * 0.000001;
}

/******************************************************************************
 * @brief    Get CPU time
 *****************************************************************************/
double
get_cpu_time()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

/******************************************************************************
 * @brief    Initialize timer values
 *****************************************************************************/
void
timer_init(timer_struct *t)
{
    t->start_wall = 0;
    t->start_cpu = 0;

    t->delta_wall = 0;
    t->delta_cpu = 0;
}

/******************************************************************************
 * @brief    Start timer
 *****************************************************************************/
void
timer_start(timer_struct *t)
{
    timer_init(t);

    t->start_wall = get_wall_time();
    t->start_cpu = get_cpu_time();
}

/******************************************************************************
 * @brief    Stop timer
 *****************************************************************************/
void
timer_stop(timer_struct *t)
{
    t->stop_wall = get_wall_time();
    t->stop_cpu = get_cpu_time();

    t->delta_wall += t->stop_wall - t->start_wall;
    t->delta_cpu += t->stop_cpu - t->start_cpu;
}

/******************************************************************************
 * @brief    Continue timer without resetting counters
 *****************************************************************************/
void
timer_continue(timer_struct *t)
{
    t->start_wall = get_wall_time();
    t->start_cpu = get_cpu_time();
}
