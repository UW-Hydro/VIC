/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for rout_stub routines
 *****************************************************************************/

#ifndef ROUT_STUB_H
#define ROUT_STUB_H

#define ROUT_EXT "rout_stub"

#include <vic_def.h>
#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief   Routing Structs
 *****************************************************************************/
typedef struct {
    size_t full_time_length;                        /*scalar - number of timesteps*/
    size_t n_outlets;                          /*scalar - length of subset*/
} rout_param_struct;

/******************************************************************************
 * @brief   main routing Struct
 *****************************************************************************/
typedef struct {
    rout_param_struct rout_param;
    double ring[1];
} rout_struct;

/******************************************************************************
 * @brief   prototypes for dummy functions of the rout_stub extension
 *****************************************************************************/
void rout_alloc(void);      // allocate memory
void rout_init(void);       // initialize model parameters from parameter files
void rout_run(void);        // run routing over the domain
void rout_finalize(void);   // clean up routine for routing
void vic_store_rout_extension(nc_file_struct *);
void vic_restore_rout_extension(nameid_struct *, metadata_struct *);
void state_metadata_rout_extension();
void set_nc_state_file_info_rout_extension(nc_file_struct *);
void set_nc_state_var_info_rout_extension(nc_file_struct *);
void initialize_state_file_rout_extension(char *, nc_file_struct *);

/******************************************************************************
 * @brief   Output state variable.
 *****************************************************************************/
enum
{
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_STATE_VARS_EXT                       /**< used as a loop counter*/
};

#endif
