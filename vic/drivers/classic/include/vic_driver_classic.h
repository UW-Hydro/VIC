/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_classic routines
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

#ifndef VIC_DRIVER_CLASSIC_H
#define VIC_DRIVER_CLASSIC_H

#include <vic_driver_shared_all.h>

#define VIC_DRIVER "Classic"

#define BINHEADERSIZE 256
#define MAX_VEGPARAM_LINE_LENGTH 500
#define ASCII_STATE_FLOAT_FMT "%.16g"

void alloc_atmos(int, atmos_data_struct **);
void alloc_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
void check_files(filep_struct *, filenames_struct *);
bool check_save_state_flag(dmy_struct *, size_t);
FILE  *check_state_file(char *, size_t, size_t, int *);
void close_files(filep_struct *filep, stream_struct **streams);
void compute_cell_area(soil_con_struct *);
void free_atmos(int nrecs, atmos_data_struct **atmos);
void free_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void free_veglib(veg_lib_struct **);
double get_dist(double lat1, double long1, double lat2, double long2);
void get_force_type(char *, int, int *);
void get_global_param(FILE *);
void initialize_forcing_files(void);
void make_in_and_outfiles(filep_struct *filep, filenames_struct *filenames,
                          soil_con_struct *soil, stream_struct **streams,
                          dmy_struct *dmy);
FILE *open_state_file(global_param_struct *, filenames_struct, size_t, size_t);
void print_atmos_data(atmos_data_struct *atmos, size_t nr);
void parse_output_info(FILE *gp, stream_struct **output_streams,
                       dmy_struct *dmy_current);
void read_atmos_data(FILE *, global_param_struct, int, int, double **,
                     double ***);
double **read_forcing_data(FILE **, global_param_struct, double ****);
void read_initial_model_state(FILE *, all_vars_struct *, int, int, int,
                              soil_con_struct *, lake_con_struct);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
void read_snowband(FILE *, soil_con_struct *);
soil_con_struct read_soilparam(FILE *, char *, char *);
veg_lib_struct *read_veglib(FILE *, size_t *);
veg_con_struct *read_vegparam(FILE *, int, size_t);
void vic_force(atmos_data_struct *, dmy_struct *, FILE **, veg_con_struct *,
               veg_hist_struct **, soil_con_struct *);
void vic_populate_model_state(all_vars_struct *, filep_struct, size_t,
                              soil_con_struct *, veg_con_struct *,
                              lake_con_struct);
void write_data(stream_struct *streams, dmy_struct *dmy);
void write_header(stream_struct **streams, dmy_struct *dmy);
void write_model_state(all_vars_struct *, int, int, filep_struct *,
                       soil_con_struct *);
void write_output(stream_struct **streams, dmy_struct *dmy);
#endif
