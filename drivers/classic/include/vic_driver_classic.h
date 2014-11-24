/* header file for vic_driver_classic routines */

#ifndef VIC_DRIVER_CLASSIC_H
#define VIC_DRIVER_CLASSIC_H

#include <stdio.h>

/*******************************************************
   Stores forcing file input information.
*******************************************************/
typedef struct {
    size_t N_ELEM; // number of elements per record; for LAI and ALBEDO,
                   // 1 element per veg tile; for others N_ELEM = 1;
    char SIGNED;
    int SUPPLIED;
    double multiplier;
} force_type_struct;

/******************************************************************
   This structure records the parameters set by the forcing file
   input routines.  Those filled, are used to estimate the paramters
   needed for the model run in initialize_atmos.c.
******************************************************************/
typedef struct {
    force_type_struct TYPE[N_FORCING_TYPES];
    unsigned short FORCE_DT[2];    /* forcing file time step */
    int FORCE_ENDIAN[2]; /* endian-ness of input file, used for
                            DAILY_BINARY format */
    int FORCE_FORMAT[2]; /* ASCII or BINARY */
    int FORCE_INDEX[2][N_FORCING_TYPES];
    unsigned N_TYPES[2];
} param_set_struct;

void alloc_atmos(int, atmos_data_struct **);
void alloc_veg_hist(int, int, veg_hist_struct ***);
void calc_longwave(double *, double, double, double);
void calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
void calc_root_fractions(veg_con_struct *, soil_con_struct *);
void check_files(filep_struct *, filenames_struct *);
FILE  *check_state_file(char *, size_t, size_t, int *);
void close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
void cmd_proc(int argc, char **argv, char *globalfilename);
void compress_files(char string[]);
void compute_treeline(atmos_data_struct *, dmy_struct *, double, double *,
                      char *);
out_data_struct *create_output_list();
void display_current_settings(int, filenames_struct *, global_param_struct *);
void free_atmos(int nrecs, atmos_data_struct **atmos);
void free_all_vars(all_vars_struct *, int);
void free_dmy(dmy_struct **dmy);
void free_out_data_files(out_data_file_struct **);
void free_out_data(out_data_struct **);
void free_vegcon(veg_con_struct **);
void free_veglib(veg_lib_struct **);
void free_veg_hist(int nrecs, int nveg, veg_hist_struct ***veg_hist);
double get_dist(double, double, double, double);
void get_force_type(char *, int, int *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void get_next_time_step(unsigned short *, unsigned short *, unsigned short *,
                        unsigned short *, unsigned short *, unsigned short);
double hermint(double, int, double *, double *, double *, double *, double *);
void hermite(int, double *, double *, double *, double *, double *);
void HourlyT(int, int, int *, double *, int *, double *, double *);
void init_output_list(out_data_struct *, int, char *, int, double);
void initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **,
                      veg_lib_struct *, veg_con_struct *, veg_hist_struct **,
                      soil_con_struct *, out_data_file_struct *,
                      out_data_struct *);
void initialize_global();
int initialize_model_state(all_vars_struct *, global_param_struct *,
                           filep_struct, size_t, size_t, size_t, double,
                           soil_con_struct *, veg_con_struct *,
                           lake_con_struct);
void initialize_snow(snow_data_struct **, size_t);
void initialize_soil(cell_data_struct **, soil_con_struct *, size_t);
void initialize_veg(veg_var_struct **, size_t);
cell_data_struct **make_cell_data(size_t);
all_vars_struct make_all_vars(size_t);
dmy_struct *make_dmy(global_param_struct *);
energy_bal_struct **make_energy_bal(size_t);
void make_in_and_outfiles(filep_struct *, filenames_struct *, soil_con_struct *,
                          out_data_file_struct *);
snow_data_struct **make_snow_data(size_t);
veg_var_struct **make_veg_var(size_t);
void mtclim_wrapper(int, int, double, double, double, double, double, double,
                    double, double, size_t, dmy_struct *, double *, double *,
                    double *, double *, double *, double *, double *);
FILE  *open_file(char string[], char type[]);
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
void usage(char *);
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
