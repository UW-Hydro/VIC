/******************************************************************************
 * @section DESCRIPTION
 *
 * Functions related to opening and closing netCDF files.
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
 * @brief    Open a netCDF file
 *****************************************************************************/
int
open_nc(char *nc_name)
{
    int    nc_id;
    int    status;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    check_nc_status(status, "Error opening %s", nc_name);

    return(nc_id);
}

/******************************************************************************
 * @brief    Close a netCDF file
 *****************************************************************************/
int
close_nc(nameid_struct nc_nameid)
{
    int status;

    // close the netcdf file
    status = nc_close(nc_nameid.nc_id);
    check_nc_status(status, "Error closing %s", nc_nameid.nc_file);

    return(status);
}
