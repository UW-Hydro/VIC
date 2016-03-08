/******************************************************************************
 * @section DESCRIPTION
 *
 * Functions to support reading from netCDF field.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Read double precision netCDF field from file.
 *****************************************************************************/
int
get_nc_field_double(char   *nc_name,
                    char   *var_name,
                    size_t *start,
                    size_t *count,
                    double *var)
{
    int nc_id;
    int status;
    int var_id;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        log_err("Error opening %s", nc_name);
    }

    /* get NetCDF variable */
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        log_err("Error getting variable id for %s in %s", var_name, nc_name);
    }

    status = nc_get_vara_double(nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        log_err("Error getting values for %s in %s", var_name, nc_name);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        log_err("Error closing %s", nc_name);
    }

    return status;
}

/******************************************************************************
 * @brief    Read single precision netCDF field from file.
 *****************************************************************************/
int
get_nc_field_float(char   *nc_name,
                   char   *var_name,
                   size_t *start,
                   size_t *count,
                   float  *var)
{
    int nc_id;
    int status;
    int var_id;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        log_err("Error opening %s", nc_name);
    }

    /* get NetCDF variable */
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        log_err("Error getting variable id for %s in %s", var_name, nc_name);
    }

    status = nc_get_vara_float(nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        log_err("Error getting values for %s in %s", var_name, nc_name);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        log_err("Error closing %s", nc_name);
    }

    return status;
}

/******************************************************************************
 * @brief    Read integer netCDF field from file.
 *****************************************************************************/
int
get_nc_field_int(char   *nc_name,
                 char   *var_name,
                 size_t *start,
                 size_t *count,
                 int    *var)
{
    int nc_id;
    int status;
    int var_id;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        log_err("Error opening %s", nc_name);
    }

    /* get NetCDF variable */
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        log_err("Error getting variable id for %s in %s", var_name, nc_name);
    }

    status = nc_get_vara_int(nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        log_err("Error getting values for %s in %s", var_name, nc_name);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        log_err("Error closing %s", nc_name);
    }

    return status;
}
