/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes atmospheric variables for both the model time step,
 * and the time step used by the snow algorithm (if different). Air temperature
 * is estimated using MTCLIM (see routine for reference), atmospheric moisture
 * is estimated using Kimball's algorithm (see routine for reference), and
 * radiation is estimated using Bras's algorithms (see routines for reference).
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize atmospheric variables for both the model time step.
 *****************************************************************************/
void
vic_force(atmos_data_struct    *atmos,
          dmy_struct           *dmy,
          FILE                **infile,
          veg_lib_struct       *veg_lib,
          veg_con_struct       *veg_con,
          veg_hist_struct     **veg_hist,
          soil_con_struct      *soil_con,
          out_data_file_struct *out_data_files,
          out_data_struct      *out_data)
{
    ;
}
