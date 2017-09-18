/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_shared_image routines
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

#ifndef VIC_DRIVER_SHARED_IMAGE_H
#define VIC_DRIVER_SHARED_IMAGE_H

#include <vic_driver_shared_all.h>
#include <vic_image_log.h>
#include <vic_mpi.h>

#include <netcdf.h>

#define MAXDIMS 10
#define AREA_SUM_ERROR_THRESH 1e-20

/******************************************************************************
 * @brief   NetCDF file types
 *****************************************************************************/
enum
{
    NC_HISTORY_FILE,
    NC_STATE_FILE,
};

/******************************************************************************
 * @brief    Structure to store location information for individual grid cells.
 * @details  The global and local indices show the position of the grid cell
 *           within the global and local (processor) domains. When the model is
 *           run on a single processor, the glonal and local domains are
 *           identical. The model is run over a list of cells.
 *****************************************************************************/
typedef struct {
    bool run; /**< TRUE: run grid cell. FALSE: do not run grid cell */
    double latitude; /**< latitude of grid cell center */
    double longitude; /**< longitude of grid cell center */
    double area; /**< area of grid cell */
    double frac; /**< fraction of grid cell that is active */
    size_t nveg; /**< number of vegetation type according to parameter file */
    size_t global_idx; /**< index of grid cell in global list of grid cells */
    size_t io_idx; /**< index of cell in 1-D I/O arrays */
    size_t local_idx; /**< index of grid cell in local list of grid cells */
} location_struct;

/******************************************************************************
 * @brief    Structure to store information about the domain file.
 *****************************************************************************/
typedef struct {
    char lat_var[MAXSTRING]; /**< latitude variable name in the domain file */
    char lon_var[MAXSTRING];  /**< longitude variable name in the domain file */
    char mask_var[MAXSTRING]; /**< mask variable name in the domain file */
    char area_var[MAXSTRING]; /**< area variable name in the domain file */
    char frac_var[MAXSTRING]; /**< fraction variable name in the domain file */
    char y_dim[MAXSTRING]; /**< y dimension name in the domain file */
    char x_dim[MAXSTRING]; /**< x dimension name in the domain file */
    size_t n_coord_dims; /**< number of x/y coordinates */
} domain_info_struct;

/******************************************************************************
 * @brief    Structure to store local and global domain information. If the
 *           model is run on a single processor, then the two are identical.
 *****************************************************************************/
typedef struct {
    size_t ncells_total; /**< total number of grid cells on domain */
    size_t ncells_active; /**< number of active grid cells on domain */
    size_t n_nx; /**< size of x-index; */
    size_t n_ny; /**< size of y-index */
    location_struct *locations; /**< locations structs for local domain */
    domain_info_struct info; /**< structure storing domain file info */
} domain_struct;

/******************************************************************************
 * @brief    Structure for netcdf variable information
 *****************************************************************************/
typedef struct {
    int nc_varid;                   /**< variable netcdf id */
    int nc_type;                    /**< variable netcdf type */
    int nc_dimids[MAXDIMS];         /**< ids of dimensions */
    size_t nc_counts[MAXDIMS];      /**< size of dimid */
    size_t nc_dims;                 /**< number of dimensions */
} nc_var_struct;

/******************************************************************************
 * @brief    Structure for netcdf file information. Initially to store
 *           information for the output files (state and history)
 *****************************************************************************/
typedef struct {
    char c_fillvalue;
    int i_fillvalue;
    double d_fillvalue;
    float f_fillvalue;
    short int s_fillvalue;
    int nc_id;
    int band_dimid;
    int front_dimid;
    int frost_dimid;
    int lake_node_dimid;
    int layer_dimid;
    int ni_dimid;
    int nj_dimid;
    int node_dimid;
    int outlet_dimid;
    int routing_timestep_dimid;
    int root_zone_dimid;
    int time_dimid;
    int time_bounds_dimid;
    int veg_dimid;
    int time_varid;
    int time_bounds_varid;
    size_t band_size;
    size_t front_size;
    size_t frost_size;
    size_t lake_node_size;
    size_t layer_size;
    size_t ni_size;
    size_t nj_size;
    size_t node_size;
    size_t outlet_size;
    size_t routing_timestep_size;
    size_t root_zone_size;
    size_t time_size;
    size_t veg_size;
    bool open;
    nc_var_struct *nc_vars;
} nc_file_struct;

/******************************************************************************
 * @brief    Structure for mapping the vegetation types for each grid cell as
 *           stored in VIC's veg_con_struct to a regular array.
 *****************************************************************************/
typedef struct {
    size_t nv_types; /**< total number of vegetation types */
                     /**< size of vidx and Cv arrays */
    size_t nv_active; /**< number of active vegetation types. Because of the */
                      /**< way that VIC defines nveg, this is nveg+1 */
                      /**< (for bare soil) or nveg+2 (if the treeline option */
                      /**< is active as well) */
    int *vidx;     /**< array of indices for active vegetation types */
    double *Cv;    /**< array of fractional coverage for nc_types */
} veg_con_map_struct;

/******************************************************************************
 * @brief   file structures
 *****************************************************************************/
typedef struct {
    FILE *globalparam;  /**< global parameters file */
    FILE *constants;    /**< model constants parameter file */
    FILE *logfile;      /**< log file */
} filep_struct;

/******************************************************************************
 * @brief   This structure stores input and output filenames.
 *****************************************************************************/
typedef struct {
    nameid_struct forcing[MAX_FORCE_FILES];  /**< atmospheric forcing files */
    char f_path_pfx[MAX_FORCE_FILES][MAXSTRING]; /**< path and prefix for
                                                    atmospheric forcing files */
    char global[MAXSTRING];     /**< global control file name */
    nameid_struct domain;       /**< domain file name and nc_id*/
    char constants[MAXSTRING];  /**< model constants file name */
    nameid_struct params;       /**< model parameters file name and nc_id */
    nameid_struct rout_params;  /**< routing parameters file name and nc_id */
    nameid_struct init_state;   /**< initial model state file name and nc_id */
    char result_dir[MAXSTRING]; /**< result directory */
    char statefile[MAXSTRING];  /**< name of model state file */
    char log_path[MAXSTRING];   /**< Location to write log file to */
} filenames_struct;

void add_nveg_to_global_domain(nameid_struct *nc_nameid,
                               domain_struct *global_domain);
void alloc_force(force_data_struct *force);
void alloc_veg_hist(veg_hist_struct *veg_hist);
double air_density(double t, double p);
double average(double *ar, size_t n);
void check_init_state_file(void);
void compare_ncdomain_with_global_domain(nameid_struct *nc_nameid);
void free_force(force_data_struct *force);
void free_veg_hist(veg_hist_struct *veg_hist);
void get_domain_type(char *cmdstr);
size_t get_global_domain(nameid_struct *domain_nc_nameid,
                         nameid_struct *param_nc_nameid,
                         domain_struct *global_domain);
void copy_domain_info(domain_struct *domain_from, domain_struct *domain_to);
void get_nc_latlon(nameid_struct *nc_nameid, domain_struct *nc_domain);
size_t get_nc_dimension(nameid_struct *nc_nameid, char *dim_name);
void get_nc_var_attr(nameid_struct *nc_nameid, char *var_name, char *attr_name,
                     char **attr);
int get_nc_var_type(nameid_struct *nc_nameid, char *var_name);
int get_nc_varndimensions(nameid_struct *nc_nameid, char *var_name);
int get_nc_field_double(nameid_struct *nc_nameid, char *var_name, size_t *start,
                        size_t *count, double *var);
int get_nc_field_float(nameid_struct *nc_nameid, char *var_name, size_t *start,
                       size_t *count, float *var);
int get_nc_field_int(nameid_struct *nc_nameid, char *var_name, size_t *start,
                     size_t *count, int *var);
int get_nc_dtype(unsigned short int dtype);
int get_nc_mode(unsigned short int format);
void initialize_domain(domain_struct *domain);
void initialize_domain_info(domain_info_struct *info);
void initialize_filenames(void);
void initialize_fileps(void);
void initialize_global_structures(void);
void initialize_history_file(nc_file_struct *nc, stream_struct *stream);
void initialize_state_file(char *filename, nc_file_struct *nc_state_file,
                           dmy_struct *dmy_state);
void initialize_location(location_struct *location);
int initialize_model_state(all_vars_struct *all_vars, size_t Nveg,
                           size_t Nnodes, double surf_temp,
                           soil_con_struct *soil_con, veg_con_struct *veg_con);
void initialize_nc_file(nc_file_struct *nc_file, size_t nvars,
                        unsigned int *varids, unsigned short int *dtypes);
void initialize_soil_con(soil_con_struct *soil_con);
void initialize_veg_con(veg_con_struct *veg_con);
void parse_output_info(FILE *gp, stream_struct **output_streams,
                       dmy_struct *dmy_current);
void print_force_data(force_data_struct *force);
void print_domain(domain_struct *domain, bool print_loc);
void print_location(location_struct *location);
void print_nc_file(nc_file_struct *nc);
void print_nc_var(nc_var_struct *nc_var);
void print_veg_con_map(veg_con_map_struct *veg_con_map);
void put_nc_attr(int nc_id, int var_id, const char *name, const char *value);
void set_force_type(char *cmdstr, int file_num, int *field);
void set_global_nc_attributes(int ncid, unsigned short int file_type);
void set_state_meta_data_info();
void set_nc_var_dimids(unsigned int varid, nc_file_struct *nc_hist_file,
                       nc_var_struct *nc_var);
void set_nc_var_info(unsigned int varid, unsigned short int dtype,
                     nc_file_struct *nc_hist_file, nc_var_struct *nc_var);
void set_nc_state_file_info(nc_file_struct *nc_state_file);
void set_nc_state_var_info(nc_file_struct *nc_state_file);
void sprint_location(char *str, location_struct *loc);
void vic_alloc(void);
void vic_finalize(void);
void vic_image_run(dmy_struct *dmy_current);
void vic_init(void);
void vic_init_output(dmy_struct *dmy_current);
void vic_restore(void);
void vic_start(void);
void vic_store(dmy_struct *dmy_state, char *state_filename);
void vic_write(stream_struct *stream, nc_file_struct *nc_hist_file,
               dmy_struct *dmy_current);
void vic_write_output(dmy_struct *dmy);
void write_vic_timing_table(timer_struct *timers, char *driver);
#endif
