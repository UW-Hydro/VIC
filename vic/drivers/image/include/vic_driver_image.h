/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_image routines
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

#ifndef VIC_DRIVER_IMAGE_H
#define VIC_DRIVER_IMAGE_H

#include <vic_mpi.h>
#include <netcdf.h>
#include <vic_driver_shared.h>

#define VIC_DRIVER "Image"

#define MAXDIMS 10

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
 * @brief    Structure to store local and global domain information. If the
 *           model is run on a single processor, then the two are identical.
 *****************************************************************************/
typedef struct {
    size_t ncells_total; /**< total number of grid cells on domain */
    size_t ncells_active; /**< number of active grid cells on domain */
    size_t n_nx; /**< size of x-index; */
    size_t n_ny; /**< size of y-index */
    location_struct *locations; /**< locations structs for local domain */
} domain_struct;

/******************************************************************************
 * @brief    Structure for netcdf file information. Initially to store
 *           information for the output files (state and history)
 *****************************************************************************/
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

/******************************************************************************
 * @brief    Structure for netcdf variable information
 *****************************************************************************/
typedef struct {
    char nc_var_name[MAXSTRING]; /**< variable name */
    char nc_units[MAXSTRING]; /**< variable name */
    int nc_dimids[MAXDIMS]; /**< ids of dimensions */
    int nc_counts[MAXDIMS]; /**< size of dimid */
    int nc_type; /**< netcdf type */
    int nc_aggtype; /**< aggregation type as defined in vic_def.h */
    int nc_dims; /**< number of dimensions */
    int nc_write; /**< TRUE: write to file; FALSE: don't */
} nc_var_struct;

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
    int *vidx;      /**< array of indices for active vegetation types */
    double *Cv;     /**< array of fractional coverage for nc_types */
} veg_con_map_struct;

void add_nveg_to_global_domain(char *nc_name, domain_struct *global_domain);
void alloc_atmos(atmos_data_struct *atmos);
void alloc_veg_hist(veg_hist_struct *veg_hist);
double air_density(double t, double p);
double average(double *ar, size_t n);
out_data_struct *create_output_list(void);
void free_atmos(atmos_data_struct *atmos);
void free_veg_hist(veg_hist_struct *veg_hist);
void get_forcing_file_info(param_set_struct *param_set, size_t file_num);
size_t get_global_domain(char *fname, domain_struct *global_domain);
void get_global_param(FILE *);
size_t get_nc_dimension(char *nc_name, char *dim_name);
void get_nc_var_attr(char *nc_name, char *var_name, char *attr_name,
                     char **attr);
int get_nc_varndimensions(char *nc_name, char *var_name);
int get_nc_field_double(char *nc_name, char *var_name, size_t *start,
                        size_t *count, double *var);
int get_nc_field_float(char *nc_name, char *var_name, size_t *start,
                       size_t *count, float *var);
int get_nc_field_int(char *nc_name, char *var_name, size_t *start,
                     size_t *count, int *var);
void initialize_domain(domain_struct *domain);
void initialize_energy(energy_bal_struct **energy, size_t nveg);
void initialize_history_file(nc_file_struct *nc);
void initialize_location(location_struct *location);
int initialize_model_state(all_vars_struct *all_vars, size_t Nveg,
                           size_t Nnodes, double surf_temp,
                           soil_con_struct *soil_con, veg_con_struct *veg_con);
void initialize_soil_con(soil_con_struct *soil_con);
void initialize_state_file(nc_file_struct *nc);
void initialize_veg_con(veg_con_struct *veg_con);
FILE *open_file(char *string, char *type);
int parse_output_info(FILE *gp, out_data_struct **out_data);
void print_atmos_data(atmos_data_struct *atmos);
void print_domain(domain_struct *domain, bool print_loc);
void print_location(location_struct *location);
void print_nc_file(nc_file_struct *nc);
void print_nc_var(nc_var_struct *nc_var, size_t ndims);
void print_veg_con_map(veg_con_map_struct *veg_con_map);
int put_nc_field_double(char *nc_name, bool *open, int *nc_id, double fillval,
                        int *dimids, int ndims, char *var_name, size_t *start,
                        size_t *count, double *var);
int put_nc_field_int(char *nc_name, bool *open, int *nc_id, int fillval,
                     int *dimids, int ndims, char *var_name, size_t *start,
                     size_t *count, int *var);
void sprint_location(char *str, location_struct *loc);
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

#endif
