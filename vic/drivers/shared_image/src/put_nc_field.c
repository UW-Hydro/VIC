/******************************************************************************
 * @section DESCRIPTION
 *
 * Put data to netCDF.
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
 * @brief    Put double precision data field to netCDF.
 *****************************************************************************/
int
put_nc_field_double(char   *nc_name,
                    bool   *open,
                    int    *nc_id,
                    double  fillval,
                    int    *dimids,
                    int     ndims,
                    char   *var_name,
                    size_t *start,
                    size_t *count,
                    double *var)
{
    int old_fill_mode;
    int status;
    int var_id;

    if (!(*open)) {
        // open the netcdf file
        status = nc_open(nc_name,
                         NC_WRITE | NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL,
                         nc_id);
        if (status != NC_NOERR) {
            log_err("Error opening %s", nc_name);
        }
        *open = true;

        // set the NC_FILL attribute
        status = nc_set_fill(*nc_id, NC_FILL, &old_fill_mode);
        if (status != NC_NOERR) {
            log_err("Error setting fill value in %s", nc_name);
        }
    }

    /* get NetCDF variable */
    status = nc_inq_varid(*nc_id, var_name, &var_id);
    if (status == NC_ENOTVAR) {
        // enter define mode
        status = nc_redef(*nc_id);
        if (status != NC_NOERR) {
            log_err("Error entering define mode in %s", nc_name);
        }
        // define the variable
        status = nc_def_var(*nc_id, var_name, NC_DOUBLE, ndims, dimids,
                            &var_id);
        if (status != NC_NOERR) {
            log_err("Error defining variable %s in %s", var_name, nc_name);
        }
        // set the fill value attribute
        status = nc_put_att_double(*nc_id, var_id, "_FillValue", NC_DOUBLE, 1,
                                   &fillval);
        if (status != NC_NOERR) {
            log_err("Error putting _FillValue attribute to %s in %s", var_name,
                    nc_name)
        }
        // leave define mode
        status = nc_enddef(*nc_id);
        if (status != NC_NOERR) {
            log_err("Error ending define mode for %s in %s", var_name, nc_name);
        }
    }
    else if (status != NC_NOERR) {
        log_err("Error getting variable id for %s in %s", var_name, nc_name);
    }

    status = nc_put_vara_double(*nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        log_err("Error writing values to %s in %s", var_name, nc_name);
    }

    // keep the file open

    return status;
}

/******************************************************************************
 * @brief    Put interger data field to netCDF.
 *****************************************************************************/
int
put_nc_field_int(char   *nc_name,
                 bool   *open,
                 int    *nc_id,
                 int     fillval,
                 int    *dimids,
                 int     ndims,
                 char   *var_name,
                 size_t *start,
                 size_t *count,
                 int    *var)
{
    int old_fill_mode;
    int status;
    int var_id;

    if (!(*open)) {
        // open the netcdf file
        status = nc_open(nc_name,
                         NC_WRITE | NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL,
                         nc_id);
        if (status != NC_NOERR) {
            log_err("Error opening %s", nc_name);
        }
        *open = true;

        // set the NC_FILL attribute
        status = nc_set_fill(*nc_id, NC_FILL, &old_fill_mode);
        if (status != NC_NOERR) {
            log_err("Error setting fill value in %s", nc_name);
        }
    }

    /* get NetCDF variable */
    status = nc_inq_varid(*nc_id, var_name, &var_id);
    if (status == NC_ENOTVAR) {
        // enter define mode
        status = nc_redef(*nc_id);
        if (status != NC_NOERR) {
            log_err("Error entering define mode in %s", nc_name);
        }
        // define the variable
        status = nc_def_var(*nc_id, var_name, NC_INT, ndims, dimids,
                            &var_id);
        if (status != NC_NOERR) {
            log_err("Error defining variable %s in %s", var_name, nc_name);
        }
        // set the fill value attribute
        status = nc_put_att_int(*nc_id, var_id, "_FillValue", NC_INT, 1,
                                &fillval);
        if (status != NC_NOERR) {
            log_err("Error putting _FillValue attribute to %s in %s", var_name,
                    nc_name);
        }
        // leave define mode
        status = nc_enddef(*nc_id);
        if (status != NC_NOERR) {
            log_err("Error ending define mode for %s in %s", var_name, nc_name);
        }
    }
    else if (status != NC_NOERR) {
        log_err("Error getting variable id for %s in %s", var_name, nc_name);
    }

    status = nc_put_vara_int(*nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        log_err("Error writing values to %s in %s", var_name, nc_name);
    }

    // keep the file open

    return status;
}
