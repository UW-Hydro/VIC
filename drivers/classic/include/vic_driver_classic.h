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

#include <stdio.h>
#include <vic_driver_shared.h>

void alloc_atmos(int, atmos_data_struct **);
void alloc_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void calc_longwave(double *, double, double, double);
void calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
void check_files(filep_struct *, filenames_struct *);
FILE  *check_state_file(char *, size_t, size_t, int *);
void close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
out_data_struct *create_output_list();
void free_atmos(int nrecs, atmos_data_struct **atmos);
void free_out_data_files(out_data_file_struct **);
void free_out_data(out_data_struct **);
void free_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
void free_veglib(veg_lib_struct **);
void get_force_type(char *, int, int *);
global_param_struct get_global_param(filenames_struct *, FILE *);
double hermint(double, int, double *, double *, double *, double *, double *);
void hermite(int, double *, double *, double *, double *, double *);
void HourlyT(int, int, int *, double *, int *, double *, double *);
void init_output_list(out_data_struct *, int, char *, int, double);
void initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **,
                      veg_lib_struct *, veg_con_struct *, veg_hist_struct **,
                      soil_con_struct *, out_data_file_struct *,
                      out_data_struct *);
int initialize_model_state(all_vars_struct *, global_param_struct *,
                           filep_struct, size_t, size_t, size_t, double,
                           soil_con_struct *, veg_con_struct *,
                           lake_con_struct);
void make_in_and_outfiles(filep_struct *, filenames_struct *, soil_con_struct *,
                          out_data_file_struct *);
void mtclim_wrapper(int, int, double, double, double, double, double, double,
                    double, double, size_t, dmy_struct *, double *, double *,
                    double *, double *, double *, double *, double *);
FILE *open_state_file(global_param_struct *, filenames_struct, int,
                      int);
void parse_output_info(FILE *, out_data_file_struct **,
                       out_data_struct *);
void print_atmos_data(atmos_data_struct *atmos, size_t nr);
void print_cell_data(cell_data_struct *cell, size_t nlayers, size_t nfrost,
                     size_t npet);
void print_dmy(dmy_struct *dmy);
void print_energy_bal(energy_bal_struct *eb, size_t nnodes, size_t nfronts);
void print_filenames(filenames_struct *fnames);
void print_filep(filep_struct *fp);
void print_force_type(force_type_struct *force_type);
void print_global_param(global_param_struct *gp);
void print_lake_con(lake_con_struct *lcon, size_t nlnodes);
void print_lake_var(lake_var_struct *lvar, size_t nlnodes, size_t nfronts,
                    size_t nlayers, size_t nnodes, size_t nfrost, size_t npet);
void print_layer_data(layer_data_struct *ldata, size_t nfrost);
void print_option(option_struct *option);
void print_out_data(out_data_struct *out, size_t nelem);
void print_out_data_file(out_data_file_struct *outf);
void print_param_set(param_set_struct *param_set);
void print_save_data(save_data_struct *save);
void print_snow_data(snow_data_struct *snow);
void print_soil_con(soil_con_struct *scon, size_t nlayers, size_t nnodes,
                    size_t nfrost, size_t nbands, size_t nzwt);
void print_veg_con(veg_con_struct *vcon, size_t nroots, char blowing, char lake,
                   char carbon, size_t ncanopy);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
void print_veg_var(veg_var_struct *vvar, size_t ncanopy);
void read_atmos_data(FILE *, global_param_struct, int, int, double **, double ***);
double **read_forcing_data(FILE **, global_param_struct, double ****);
void read_initial_model_state(FILE *, all_vars_struct *, int, int, int,
                              soil_con_struct *);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
void read_snowband(FILE *, soil_con_struct *);
soil_con_struct read_soilparam(FILE *, char *, char *);
veg_lib_struct *read_veglib(FILE *, size_t *);
veg_con_struct *read_vegparam(FILE *, int, size_t);
void set_max_min_hour(double *, int, int *, int *);
out_data_file_struct *set_output_defaults(out_data_struct *);
int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *,
                   int, char *, int, double);
void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *, int);
void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *,
                        out_data_struct *);
void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *,
                  global_param_struct);
void write_model_state(all_vars_struct *, int, int, filep_struct *,
                       soil_con_struct *);
void write_output(out_data_struct *out_data,
                  out_data_file_struct *out_data_files, dmy_struct *dmy,
                  int rec);
#endif
