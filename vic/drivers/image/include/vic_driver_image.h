/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_image routines
 *****************************************************************************/

#ifndef VIC_DRIVER_IMAGE_H
#define VIC_DRIVER_IMAGE_H

#include <vic_driver_shared_image.h>

#define VIC_DRIVER "Image"

bool check_save_state_flag(size_t, dmy_struct *dmy_offset);
void display_current_settings(int);
void get_forcing_file_info(param_set_struct *param_set, size_t file_num);
void get_global_param(FILE *);
void vic_force(void);
void vic_image_init(void);
void vic_image_finalize();
void vic_image_start(void);
void vic_populate_model_state(dmy_struct *dmy_current);

#endif
