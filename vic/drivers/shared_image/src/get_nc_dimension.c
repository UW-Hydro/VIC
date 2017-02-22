/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF dimension.
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
 * @brief    Get netCDF dimension.
 *****************************************************************************/
size_t
get_nc_dimension(nameid_struct *nc_nameid,
                 char          *dim_name)
{
    int    dim_id;
    size_t dim_size;
    int    status;

    // get dimension id
    status = nc_inq_dimid(nc_nameid->nc_id, dim_name, &dim_id);
    check_nc_status(status, "Error getting dimension id %s in %s", dim_name,
                    nc_nameid->nc_filename);

    // get dimension size
    status = nc_inq_dimlen(nc_nameid->nc_id, dim_id, &dim_size);
    check_nc_status(status, "Error getting dimension size for dim %s in %s",
                    dim_name,
                    nc_nameid->nc_filename);

    return dim_size;
}
