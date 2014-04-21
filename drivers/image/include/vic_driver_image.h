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

void cmd_proc(int argc, char **argv, char *globalfilename);
void display_current_settings(int, filenames_struct *, global_param_struct *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void initialize_global();
FILE *open_file(char *string, char *type);
void vic_alloc(void);
void vic_start(void);
void usage(char *);

#endif
