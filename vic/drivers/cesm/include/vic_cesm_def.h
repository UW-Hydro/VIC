/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for defining CESM and ESMF types in C
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

#ifndef VIC_CESM_DEF_H
#define VIC_CESM_DEF_H

/******************************************************************************
 * @brief   This structure stores clock information. See also ESMF_Clock.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_interface_f.F90
 *****************************************************************************/
typedef struct {
    int       timestep; // timestep in seconds
    short int start_year;  // start year
    short int start_month;  // start month
    short int start_day;  // start day
    int       start_dayseconds;  // start dayseconds
    short int stop_year;  // stop year
    short int stop_month;  // stop month
    short int stop_day;  // stop day
    int       stop_dayseconds;  // stop dayseconds
    short int current_year;  // current year
    short int current_month;  // current month
    short int current_day;  // current day
    int       current_dayseconds;  // current dayseconds
    short int previous_year;  // previous year
    short int previous_month;  // previous month
    short int previous_day;  // previous day
    int       previous_dayseconds;  // previous dayseconds
    bool      state_flag; // state flag
    bool      stop_flag; // stop flag
} vic_clock_struct;

/******************************************************************************
 * @brief   This structure stores mpi and domain information.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_interface_f.F90
 *****************************************************************************/
// typedef struct {
//     int id; // component id
//     int mpicom; //mpi communicator

//     sometype dom; // domain info
//     sometype gsmap; // decomp info
//     sometype infodata; // input init object

// } vic_seq_cdata_struct;

#endif

int vic_cesm_init(vic_clock_struct vic_clock,
                  char *vic_global_param_file);
int vic_cesm_run(vic_clock_struct vic_clock);
int vic_cesm_final(void);
