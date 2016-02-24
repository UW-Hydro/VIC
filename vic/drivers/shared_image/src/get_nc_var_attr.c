/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF variable attribute.
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
 * @brief    Get netCDF dimension.
 *****************************************************************************/
void
get_nc_var_attr(char  *nc_name,
                char  *var_name,
                char  *attr_name,
                char **attr)
{
    int    nc_id;
    int    var_id;
    int    status;
    size_t attr_len;

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

    // get size of the attribute
    status = nc_inq_attlen(nc_id, var_id, attr_name, &attr_len);
    if (status != NC_NOERR) {
        log_err("Error getting attribute length for %s:%s in %s", var_name,
                attr_name, nc_name);
    }

    // allocate memory for attribute
    *attr = malloc((attr_len + 1) * sizeof(**attr));

    // read attribute text
    status = nc_get_att_text(nc_id, var_id, attr_name, *attr);
    if (status != NC_NOERR) {
        log_err("Error getting netCDF attribute %s for var %s in %s", attr_name,
                var_name, nc_name);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        log_err("Error closing %s", nc_name);
    }
}
