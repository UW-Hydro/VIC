/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_image routines
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
