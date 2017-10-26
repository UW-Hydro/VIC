/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for rvic routing routines
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
#ifndef ROUT_RVIC_H
#define ROUT_RVIC_H

#define ROUT_EXT "rout_rvic"

#include <vic_def.h>
#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief   Routing Structs
 *****************************************************************************/
typedef struct {
    size_t full_time_length;                  /*scalar - total number of timesteps*/
    size_t n_timesteps;                        /*scalar - number of timesteps*/
    size_t n_sources;                          /*scalar - number of sources*/
    size_t n_outlets;                          /*scalar - length of subset*/
    size_t *source2outlet_ind;                /*1d array - source to outlet mapping*/
    int *source_y_ind;                        /*1d array - source y location*/
    int *source_x_ind;                        /*1d array - source x location*/
    double *source_lat;                       /*1d array - Latitude coordinate of source grid cell*/
    double *source_lon;                       /*1d array - Longitude coordinate of source grid cell*/
    double *outlet_lat;                       /*1d array - Latitude coordinate of outlet grid cell*/
    double *outlet_lon;                       /*1d array - Longitude coordinate of outlet grid cell*/
    int *source_VIC_index;                    /*1d array - mapping of routing-source index to VIC index*/
    int *outlet_VIC_index;                    /*1d array - mapping of routing-outlet index to VIC index*/
    int *source_time_offset;                  /*1d array - source time offset*/
    double *unit_hydrograph;                  /*2d array[times][sources] - unit hydrographs*/
    double *aggrunin;                         /*2d array[ysize][xsize] - vic runoff flux*/
} rout_param_struct;

/******************************************************************************
 * @brief   main routing Struct
 *****************************************************************************/
typedef struct {
    rout_param_struct rout_param;
    double *ring;
    double *discharge;
} rout_struct;

/******************************************************************************
 * @brief   Function prototypes for the rout_rvic extension
 *****************************************************************************/
void rout_alloc(void);                 // allocate memory
void rout_init(void);                  // initialize model parameters from parameter files
void rout_run(void);                   // run routing over the domain
void rout_finalize(void);              // clean up routine for routing
void convolution(double *, double *);  // convolution over the domain

/******************************************************************************
 * @brief   MPI Function prototypes for the rout_rvic extension
 *****************************************************************************/
void scatter_var_double(double *, double *);
void gather_var_double(double *, double *);

/******************************************************************************
 * @brief   Convolution function adapted from the RVIC scheme
 *****************************************************************************/
void get_global_param_rout(FILE *gp);
void cshift(double *, int, int, int, int);
void vic_store_rout_extension(nc_file_struct *);
void vic_restore_rout_extension(nameid_struct *, metadata_struct *);
void state_metadata_rout_extension();
void set_nc_state_file_info_rout_extension(nc_file_struct *);
void set_nc_state_var_info_rout_extension(nc_file_struct *);
void initialize_state_file_rout_extension(char *, nc_file_struct *);

/******************************************************************************
 * @brief   Output state variable.
 *****************************************************************************/
enum
{
    STATE_ROUT_RING,                   /**<  routing ring: rout_ring[routing_timestep, outlet] */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_STATE_VARS_EXT                       /**< used as a loop counter*/
};

#endif
