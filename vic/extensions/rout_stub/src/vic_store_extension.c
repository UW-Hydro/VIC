/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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

#include <vic_driver_shared_image.h>
#include <rout.h>

/******************************************************************************
 * @brief    Save model state.
 *****************************************************************************/
void
vic_store_rout_extension(nc_file_struct *nc_state_file)
{
}

/******************************************************************************
 * @brief   Setup state file netcdf structure
 *****************************************************************************/
void
set_nc_state_file_info_rout_extension(nc_file_struct *nc_state_file)
{
}

/******************************************************************************
 * @brief   Setup state variable dimensions, types, etc.
 *****************************************************************************/
void
set_nc_state_var_info_rout_extension(nc_file_struct *nc)
{
}

/******************************************************************************
 * @brief   Initialize state file by creating dimensions, variables,
            and adding metadata.
 *****************************************************************************/
void
initialize_state_file_rout_extension(char           *filename,
                                     nc_file_struct *nc_state_file)
{
}
