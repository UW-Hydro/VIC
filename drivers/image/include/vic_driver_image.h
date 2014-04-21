/* header file for vic_driver_image routines */

#ifndef VIC_DRIVER_IMAGE_H
#define VIC_DRIVER_IMAGE_H

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
    float latitude; // latitude of grid cell center
    float longitude; // longitude of grid cell center
    long int global_cellidx; // index of grid cell in global list of grid cells
    long int global_latidx; // index of latitude in global domain
    long int global_lonidx; // index of longitude in global domain
    long int local_cellidx; // index of grid cell in local list of grid cells
    long int local_latidx; // index of latitude in local domain
    long int local_lonidx; // index of longitude in local domain
} location_struct;


/*
   Structure to store local and global domain information. If the model is run on
   a single processor, then the two are identical.
 */
typedef struct {
    long int ncells_global; // number of active grid cell on global domain
    float global_ul_lat; // lat of upper left grid cell center of global domain
    float global_ul_lon; // lon of upper left grid cell center of global domain
    float global_lr_lat; // lat of lower right grid cell center of global domain
    float global_lr_lon; // lon of lower right grid cell center of global domain
    long int ncells_local; // number of active grid cell on local domain
    float local_ul_lat; // lat of upper left grid cell center of local domain
    float local_ul_lon; // lon of upper left grid cell center of local domain
    float local_lr_lat; // lat of lower right grid cell center of local domain
    float local_lr_lon; // lon of lower right grid cell center of local domain
} domain_struct;


void cmd_proc(int argc, char **argv, char *globalfilename);
void display_current_settings(int, filenames_struct *, global_param_struct *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void initialize_global();
FILE *open_file(char *string, char *type);
void vic_alloc(void);
void vic_start(void);
void usage(char *);

#endif
