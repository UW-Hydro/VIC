/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF dimension.
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
#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Get netCDF dimension.
 *****************************************************************************/
int
get_nc_varndimensions(char *nc_name,
                      char *var_name)
{
    int nc_id;
    int var_id;
    int ndims;
    int status;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        log_err("Error opening %s", nc_name);
    }

    // get variable id
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        log_err("Error getting variable id %s in %s", var_name, nc_name);
    }

    // get number of dimensions
    status = nc_inq_varndims(nc_id, var_id, &ndims);
    if (status != NC_NOERR) {
        log_err("Error getting number of dimensions for var %s in %s", var_name,
                nc_name);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        log_err("Error closing %s", nc_name);
    }

    return ndims;
}
