/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF variable attribute.
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

/******************************************************************************
 * @brief    Get netCDF variable attributes.
 *****************************************************************************/
void
get_nc_var_attr(nameid_struct *nc_nameid,
                char          *var_name,
                char          *attr_name,
                char         **attr)
{
    int    var_id;
    int    status;
    size_t attr_len;

    // get variable id
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id %s in %s", var_name,
                    nc_nameid->nc_filename);

    // get size of the attribute
    status = nc_inq_attlen(nc_nameid->nc_id, var_id, attr_name, &attr_len);
    check_nc_status(status, "Error getting attribute length for %s:%s in %s",
                    var_name,
                    attr_name, nc_nameid->nc_filename);

    // allocate memory for attribute
    *attr = malloc((attr_len + 1) * sizeof(**attr));
    check_alloc_status(*attr, "Memory allocation error.");

    // read attribute text
    status = nc_get_att_text(nc_nameid->nc_id, var_id, attr_name, *attr);
    check_nc_status(status,
                    "Error getting netCDF attribute %s for var %s in %s",
                    attr_name,
                    var_name, nc_nameid->nc_filename);

    // we need to null terminate the string ourselves according to NetCDF docs
    (*attr)[attr_len] = '\0';
}
