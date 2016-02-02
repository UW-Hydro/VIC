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

#include <vic_driver_shared.h>

#define VIC_DRIVER "Classic"

#define BINHEADERSIZE 256
#define MAX_VEGPARAM_LINE_LENGTH 500

void alloc_atmos(int, atmos_data_struct **);
void alloc_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
void check_files(filep_struct *, filenames_struct *);
FILE  *check_state_file(char *, size_t, size_t, int *);
void close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
size_t count_force_vars(FILE *gp);
size_t count_n_outfiles(FILE *gp);
size_t count_outfile_nvars(FILE *gp);
out_data_struct *create_output_list();
void free_atmos(int nrecs, atmos_data_struct **atmos);
void free_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void free_veglib(veg_lib_struct **);
void get_force_type(char *, int, int *);
void get_global_param(FILE *);
void init_output_list(out_data_struct *, int, char *, int, double);
void initialize_forcing_files(void);
int initialize_model_state(all_vars_struct *, global_param_struct *,
                           filep_struct, size_t, size_t, size_t, double,
                           soil_con_struct *, veg_con_struct *,
                           lake_con_struct);
void make_in_and_outfiles(filep_struct *, filenames_struct *, soil_con_struct *,
                          out_data_file_struct *);
FILE *open_state_file(global_param_struct *, filenames_struct, int, int);
void print_atmos_data(atmos_data_struct *atmos, size_t nr);
void parse_output_info(FILE *, out_data_file_struct **, out_data_struct *);
void read_atmos_data(FILE *, global_param_struct, int, int, double **,
                     double ***);
double **read_forcing_data(FILE **, global_param_struct, double ****);
void read_initial_model_state(FILE *, all_vars_struct *, int, int, int,
                              soil_con_struct *);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
void read_snowband(FILE *, soil_con_struct *);
soil_con_struct read_soilparam(FILE *, char *, char *);
veg_lib_struct *read_veglib(FILE *, size_t *);
veg_con_struct *read_vegparam(FILE *, int, size_t);
out_data_file_struct *set_output_defaults(out_data_struct *);
void vic_force(atmos_data_struct *, dmy_struct *, FILE **, veg_con_struct *,
               veg_hist_struct **, soil_con_struct *);
void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *,
                double);
void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *,
                        out_data_struct *, dmy_struct *);
void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *,
                  global_param_struct);
void write_model_state(all_vars_struct *, int, int, filep_struct *,
                       soil_con_struct *);
void write_output(out_data_struct *out_data,
                  out_data_file_struct *out_data_files, dmy_struct *dmy,
                  int rec);
#endif
