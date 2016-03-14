/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for rout_stub routines
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
#ifndef ROUT_STUB_H
#define ROUT_STUB_H

#define ROUT_EXT "rout_stub"

/******************************************************************************
 * @brief   main routing Struct
 *****************************************************************************/
typedef struct {
} rout_struct;

/******************************************************************************
 * @brief   prototypes for dummy functions of the rout_stub extension
 *****************************************************************************/
void rout_start(void);      // read global parameters for routing
void rout_alloc(void);      // allocate memory
void rout_init(void);       // initialize model parameters from parameter files
void rout_run(void);        // run routing over the domain
void rout_write(void);      // write routine for routing
void rout_finalize(void);   // clean up routine for routing
#endif
