/* header file for vic_driver_classic routines */

#ifndef VIC_DRIVER_CLASSIC_H
#define VIC_DRIVER_CLASSIC_H

#include <stdio.h>

void alloc_atmos(int, atmos_data_struct **);
void   calc_longwave(double *, double, double, double);
void   calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
void   calc_root_fractions(veg_con_struct *, soil_con_struct *);
void   check_files(filep_struct *, filenames_struct *);
FILE  *check_state_file(char *, dmy_struct *, global_param_struct *, int, int, 
                        int *);
void   close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
void   cmd_proc(int argc, char **argv, char *globalfilename);
void   compress_files(char string[]);
void   compute_treeline(atmos_data_struct *, dmy_struct *, double, double *, char *);
out_data_struct *create_output_list();
void   display_current_settings(int, filenames_struct *, global_param_struct *);
void   free_atmos(int nrecs, atmos_data_struct **atmos);
void   free_all_vars(all_vars_struct *, int);
void   free_dmy(dmy_struct **dmy);
void   free_out_data_files(out_data_file_struct **);
void   free_out_data(out_data_struct **);
void   free_vegcon(veg_con_struct **);
void   free_veglib(veg_lib_struct **);
double get_dist(double, double, double, double);
void   get_force_type(char *, int, int *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void   get_next_time_step(int *, int *, int *, int *, int *, int);
double hermint(double, int, double *, double *, double *, double *, double *);
void   hermite(int, double *, double *, double *, double *, double *);
void   HourlyT(int, int, int *, double *, int *, double *, double *);
void   init_output_list(out_data_struct *, int, char *, int, float);
void   initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **,
                        soil_con_struct *, out_data_file_struct *,
                        out_data_struct *);
void   initialize_global();
int    initialize_model_state(all_vars_struct *, dmy_struct,
                              global_param_struct *, filep_struct, 
                              int, int, int,  
                              double, soil_con_struct *,
                              veg_con_struct *, lake_con_struct);
void   initialize_snow(snow_data_struct **, int, int);
void   initialize_soil(cell_data_struct **, soil_con_struct *, veg_con_struct *, int);
void   initialize_veg(veg_var_struct **, veg_con_struct *,
                      global_param_struct *, int);
cell_data_struct **make_cell_data(int, int);
all_vars_struct make_all_vars(int);
dmy_struct *make_dmy(global_param_struct *);
energy_bal_struct **make_energy_bal(int);
void make_in_and_outfiles(filep_struct *, filenames_struct *, 
                          soil_con_struct *, out_data_file_struct *);
snow_data_struct **make_snow_data(int);
veg_var_struct **make_veg_var(int);
void mtclim_wrapper(int, int, double, double, double, double,
                    double, double, double, double,
                    int, dmy_struct *, double *,
                    double *, double *, double *, double *, double *, double *);
FILE  *open_file(char string[], char type[]);
FILE  *open_state_file(global_param_struct *, filenames_struct, int, int);
void parse_output_info(filenames_struct *, FILE *, out_data_file_struct **, out_data_struct *);
void   read_atmos_data(FILE *, global_param_struct, int, int, double **);
double **read_forcing_data(FILE **, global_param_struct);
void   read_initial_model_state(FILE *, all_vars_struct *, 
                                global_param_struct *, int, int, int, 
                                soil_con_struct *, lake_con_struct);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
void   read_snowband(FILE *, soil_con_struct *);
soil_con_struct read_soilparam(FILE *, char *, char *);
veg_lib_struct *read_veglib(FILE *, int *);
veg_con_struct *read_vegparam(FILE *, int, int);
void set_max_min_hour(double *, int, int *, int *);
out_data_file_struct *set_output_defaults(out_data_struct *);
int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *, int, char *, int, float);
int update_thermal_nodes(all_vars_struct *, 
                         int, int, soil_con_struct *, veg_con_struct *);
void usage(char *);
void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *, int);
void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *, out_data_struct *);
void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *, global_param_struct);
void write_model_state(all_vars_struct *, global_param_struct *, int, 
                       int, filep_struct *, soil_con_struct *, lake_con_struct);

#endif
