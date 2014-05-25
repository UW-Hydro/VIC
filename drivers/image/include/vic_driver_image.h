/* header file for vic_driver_image routines */

#ifndef VIC_DRIVER_IMAGE_H
#define VIC_DRIVER_IMAGE_H

#include <stdbool.h>
#include <netcdf.h>
#include <stdio.h>


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

void calc_root_fractions(veg_con_struct *veg_con, soil_con_struct *soil_con);
void cmd_proc(int argc, char **argv, char *globalfilename);
void display_current_settings(int, filenames_struct *, global_param_struct *);
void free_dmy(dmy_struct **dmy);
size_t get_global_domain(char *fname, domain_struct *global_domain);
global_param_struct get_global_param(filenames_struct *, FILE *);
size_t get_nc_dimension(char *nc_name, char *dim_name);
int get_nc_field_double(char *nc_name, char *var_name, size_t *start,
                        size_t *count, double *var);
int get_nc_field_int(char *nc_name, char *var_name, size_t *start,
                     size_t *count, int *var);
void initialize_domain(domain_struct *domain);
void initialize_location(location_struct *location);
void initialize_global(void);
dmy_struct *make_dmy(global_param_struct *);
FILE *open_file(char *string, char *type);
void print_domain(domain_struct *domain, bool print_loc);
void print_location(location_struct *location);
void vic_alloc(void);
void vic_finalize(void);
void vic_init(void);
void vic_start(void);
void usage(char *);

#endif
