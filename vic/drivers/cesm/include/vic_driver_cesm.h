/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for CESM-VIC routines
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

#ifndef VIC_DRIVER_CESM_H
#define VIC_DRIVER_CESM_H

#include <stdbool.h>
#include <vic_mpi.h>
#include <netcdf.h>
#include <vic_driver_shared.h>
#include <vic_cesm_def.h>

void add_nveg_to_global_domain(char *nc_name, domain_struct *global_domain);
void advance_time(void);
void alloc_atmos(atmos_data_struct *atmos);
void alloc_veg_hist(veg_hist_struct *veg_hist);
double air_density(double t, double p);
void assert_time_insync(vic_clock *vclock, dmy_struct *dmy);
double average(double *ar, size_t n);
out_data_struct *create_output_list(void);
void free_atmos(atmos_data_struct *atmos);
void free_out_data(out_data_struct **out_data);
void free_veg_hist(veg_hist_struct *veg_hist);
size_t get_global_domain(char *fname, domain_struct *global_domain);
void get_global_param(FILE *);
size_t get_nc_dimension(char *nc_name, char *dim_name);
int get_nc_field_double(char *nc_name, char *var_name, size_t *start,
                        size_t *count, double *var);
int get_nc_field_float(char *nc_name, char *var_name, size_t *start,
                       size_t *count, float *var);
int get_nc_field_int(char *nc_name, char *var_name, size_t *start,
                     size_t *count, int *var);
void init_output_list(out_data_struct *out_data, int write, char *format,
                      int type, double mult);
void initialize_cesm_time(void);
void initialize_domain(domain_struct *domain);
void initialize_energy(energy_bal_struct **energy, size_t nveg);
void initialize_history_file(nc_file_struct *nc);
void initialize_l2x_data(void);
void initialize_location(location_struct *location);
int initialize_model_state(all_vars_struct *all_vars, size_t Nveg,
                           size_t Nnodes, double surf_temp,
                           soil_con_struct *soil_con, veg_con_struct *veg_con);
void initialize_soil_con(soil_con_struct *soil_con);
void initialize_state_file(nc_file_struct *nc);
void initialize_veg_con(veg_con_struct *veg_con);
void initialize_x2l_data(void);
FILE *open_file(char *string, char *type);
int parse_output_info(FILE *gp, out_data_struct **out_data);
void print_atmos_data(atmos_data_struct *atmos);
void print_domain(domain_struct *domain, bool print_loc);
void print_l2x_data(l2x_data_struct *l2x);
void print_location(location_struct *location);
void print_nc_file(nc_file_struct *nc);
void print_nc_var(nc_var_struct *nc_var, size_t ndims);
void print_vic_clock(vic_clock *vclock);
void print_veg_con_map(veg_con_map_struct *veg_con_map);
int put_nc_field_double(char *nc_name, bool *open, int *nc_id, double fillval,
                        int *dimids, int ndims, char *var_name, size_t *start,
                        size_t *count, double *var);
int put_nc_field_int(char *nc_name, bool *open, int *nc_id, int fillval,
                     int *dimids, int ndims, char *var_name, size_t *start,
                     size_t *count, int *var);
void print_x2l_data(x2l_data_struct *x2l);
double q_to_vp(double q, double p);
void read_rpointer_file(char *fname);
void sprint_location(char *str, location_struct *loc);
unsigned short int start_type_from_char(char *start_str);
void vic_alloc(void);
int vic_cesm_init_mpi(int MPI_COMM_VIC_F);
int vic_cesm_init(char *vic_global_param_file, char *caseid, char *runtype,
                  vic_clock vclock);
int vic_cesm_final(void);
int vic_cesm_run(vic_clock vclock);
void vic_nc_info(nc_file_struct *nc_hist_file, out_data_struct **out_data,
                 nc_var_struct *nc_vars);
void vic_finalize(void);
void vic_force(void);
void vic_cesm_put_data(void);
void vic_cesm_run_model(void);
void vic_init(void);
void vic_init_output(void);
void vic_restore(char *runtype_str);
void vic_start(void);
void vic_store(void);
void vic_write(void);
char will_it_snow(double *t, double t_offset, double max_snow_temp,
                  double *prcp, size_t n);
void write_rpointer_file(char *fname);

#endif
