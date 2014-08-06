/* header file for vic_driver_image routines */

#ifndef VIC_DRIVER_IMAGE_H
#define VIC_DRIVER_IMAGE_H

#include <stdbool.h>
#include <netcdf.h>
#include <stdio.h>

#define MAXDIMS 10

/*******************************************************
   Stores forcing file input information.
*******************************************************/
typedef struct {
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
    int FORCE_DT[2];    /* forcing file time step */
    int FORCE_ENDIAN[2]; /* endian-ness of input file, used for
                            DAILY_BINARY format */
    int FORCE_FORMAT[2]; /* ASCII or BINARY */
    int FORCE_INDEX[2][N_FORCING_TYPES];
    int N_TYPES[2];
} param_set_struct;

/*
   Structure to store location information for individual grid cells.

   The global and local indices show the position of the grid cell within the
   global and local (processor) domains. When the model is run on a single
   processor, the glonal and local domains are identical. The model is run over a
   list of cells.

   cellidx - varies from 0 to total number of cells (-1)

   latidx - index of latitude of grid cell in 2-D image

   lonidx - index of longitude of grid cell in 2-D image.

 */
typedef struct {
    double latitude; // latitude of grid cell center
    double longitude; // longitude of grid cell center
    double area; // area of grid cell
    double frac; // fraction of grid cell that is active
    size_t global_cell_idx; // index of grid cell in global list of grid cells
    size_t global_x_idx; // index of x-dimension in global domain
    size_t global_y_idx; // index of y-dimension in global domain
    size_t local_cell_idx; // index of grid cell in local list of grid cells
    size_t local_x_idx; // index of x-dimension in local domain
    size_t local_y_idx; // index of y-dimension in local domain
} location_struct;


/*
   Structure to store local and global domain information. If the model is run on
   a single processor, then the two are identical. Note that this
 */
typedef struct {
    size_t ncells_global; // number of active grid cell on global domain
    size_t n_nx; // size of x-index;
    size_t n_ny; // size of y-index
    size_t ncells_local; // number of active grid cell on local domain
    location_struct *locations; // locations structs for local domain
} domain_struct;

/* Structure for netcdf file information. Initially to store information for the
   output files (state and history) */
typedef struct {
    char fname[MAXSTRING + 1];
    char c_fillvalue;
    int i_fillvalue;
    double d_fillvalue;
    float f_fillvalue;
    int nc_id;
    int band_dimid;
    int front_dimid;
    int frost_dimid;
    int layer_dimid;
    int ni_dimid;
    int nj_dimid;
    int node_dimid;
    int root_zone_dimid;
    int time_dimid;
    int veg_dimid;
    size_t band_size;
    size_t front_size;
    size_t frost_size;
    size_t layer_size;
    size_t ni_size;
    size_t nj_size;
    size_t node_size;
    size_t root_zone_size;
    size_t time_size;
    size_t veg_size;
    bool open;
} nc_file_struct;

/* Structure for netcdf variable information */
typedef struct {
    char nc_var_name[MAXSTRING]; // variable name
    char nc_units[MAXSTRING]; // variable name
    int nc_dimids[MAXDIMS]; // ids of dimensions
    int nc_counts[MAXDIMS]; // size of dimid
    int nc_type; // netcdf type
    int nc_aggtype; // aggregation type as defined in vic_def.h
    int nc_dims; // number of dimensions
    int nc_write; // TRUE: write to file; FALSE: don't
} nc_var_struct;

/* Structure for mapping the vegetation types for each grid cell as stored
   in VIC's veg_con_struct to a regular array. This is convoluted, but will
   allow us to keep using VIC's existing memory layout */
typedef struct {
    size_t nv_types; // total number of vegetation types
                     // size of vidx and Cv arrays
    size_t nv_active; // number of active vegetation types. Because of the
                      // way that VIC defines nveg, this is nveg+1
                      // (for bare soil) or nveg+2 (if the treeline option
                      // is active as well)
    int *vidx;      // array of indices for active vegetation types
    double *Cv;     // array of fractional coverage for nc_types
} veg_con_map_struct;

void alloc_atmos(atmos_data_struct *atmos);
double air_density(double t, double p, double vp);
double average(double *ar, size_t n);
void calc_root_fractions(veg_con_struct *veg_con, soil_con_struct *soil_con);
void cmd_proc(int argc, char **argv, char *globalfilename);
out_data_struct *create_output_list(void);
void display_current_settings(int, filenames_struct *, global_param_struct *);
void free_all_vars(all_vars_struct *all_vars, int Nveg);
void free_atmos(atmos_data_struct *atmos);
void free_dmy(dmy_struct **dmy);
void free_out_data(out_data_struct **out_data);
size_t get_global_domain(char *fname, domain_struct *global_domain);
size_t get_global_idx(domain_struct *domain, size_t i);
global_param_struct get_global_param(filenames_struct *, FILE *);
size_t get_nc_dimension(char *nc_name, char *dim_name);
int get_nc_field_double(char *nc_name, char *var_name, size_t *start,
                        size_t *count, double *var);
int get_nc_field_float(char *nc_name, char *var_name, size_t *start,
                       size_t *count, float *var);
int get_nc_field_int(char *nc_name, char *var_name, size_t *start,
                     size_t *count, int *var);
void get_next_time_step(int *year, int *month, int *day, int *hr, int *jday,
                        int dt);
void init_output_list(out_data_struct *out_data, int write, char *format,
                      int type, float mult);
void initialize_domain(domain_struct *domain);
void initialize_energy(energy_bal_struct **energy, soil_con_struct *soil_con,
                       int nveg);
void initialize_history_file(nc_file_struct *nc);
void initialize_location(location_struct *location);
void initialize_global(void);
int initialize_model_state(all_vars_struct *all_vars, int Nveg, int Nnodes,
                           double surf_temp, soil_con_struct *soil_con,
                           veg_con_struct *veg_con, veg_lib_struct *veg_lib);
void initialize_snow(snow_data_struct **snow, int veg_num);
void initialize_soil(cell_data_struct **cell, soil_con_struct *soil_con,
                     veg_con_struct *veg_con, int veg_num);
void initialize_soil_con(soil_con_struct *soil_con);
void initialize_state_file(nc_file_struct *nc);
void initialize_veg(veg_var_struct **veg_var, veg_con_struct *veg_con,
                    int nveg);
void initialize_veg_con(veg_con_struct *veg_con);
all_vars_struct make_all_vars(int nveg);
cell_data_struct **make_cell_data(int veg_type_num, int Nlayer);
dmy_struct *make_dmy(global_param_struct *);
energy_bal_struct **make_energy_bal(int nveg);
snow_data_struct **make_snow_data(int nveg);
veg_var_struct **make_veg_var(int veg_type_num);
FILE *open_file(char *string, char *type);
int parse_output_info(FILE *gp, out_data_struct **out_data);
void print_domain(domain_struct *domain, bool print_loc);
void print_force_type(force_type_struct *force_type);
void print_location(location_struct *location);
void print_param_set(param_set_struct *param_set);
void print_veg_con(veg_con_struct *vcon, size_t nroots, char blowing, char lake,
                   char carbon, size_t ncanopy);
void print_veg_con_map(veg_con_map_struct *veg_con_map);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
int put_nc_field_double(char *nc_name, bool *open, int *nc_id,
                        double fillval, int *dimids, int ndims,
                        char *var_name, size_t *start, size_t *count,
                        double *var);
int put_nc_field_int(char *nc_name, bool *open, int *nc_id, int fillval,
                     int *dimids, int ndims, char *var_name, size_t *start,
                     size_t *count, int *var);
double q_to_vp(double q, double p);
void soil_moisture_from_water_table(soil_con_struct *soil_con, int nlayers);
int update_thermal_nodes(all_vars_struct *all_vars, int Nveg, int Nnodes,
                         soil_con_struct *soil_con, veg_con_struct  *veg_con,
                         veg_lib_struct *veg_lib);
void vic_alloc(void);
void vic_nc_info(nc_file_struct *nc_hist_file, out_data_struct **out_data,
                 nc_var_struct *nc_vars); 
void vic_finalize(void);
void vic_force(void);
void vic_image_run(void);
void vic_init(void);
void vic_init_output(void);
void vic_restore(void);
void vic_start(void);
void vic_store(void);
void vic_write(void);
char will_it_snow(double *t, double t_offset, double max_snow_temp,
                  double *prcp, int n);
void usage(char *);

#endif
