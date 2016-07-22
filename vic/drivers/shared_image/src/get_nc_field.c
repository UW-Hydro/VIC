/******************************************************************************
 * @section DESCRIPTION
 *
 * Functions to support reading from netCDF field.
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
 * @brief    Read double precision netCDF field from file.
 *****************************************************************************/
int
get_nc_field_double(char   *nc_name,
                    char   *var_name,
                    size_t *start,
                    size_t *count,
                    double *var)
{
    int              status;
    int              var_id;

    extern nc_struct netcdf;

    if (!netcdf.name) {
        netcdf.name = nc_name;
        status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
        check_nc_status(status, "Error opening %s", netcdf.name);
    }
    else {
        // check if another netcdf has to be opened
        if (strcmp(netcdf.name, nc_name) != false) {
            // close the old netcdf file
            status = nc_close(netcdf.id);
            check_nc_status(status, "Error closing %s", netcdf.name);
            netcdf.id = 0;
            netcdf.name = '\0';

            // open the new netcdf file
            status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
            check_nc_status(status, "Error opening %s", nc_name);
            // set name
            netcdf.name = nc_name;
        }
    }

    /* get NetCDF variable */
    status = nc_inq_varid(netcdf.id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    netcdf.name);

    status = nc_get_vara_double(netcdf.id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    netcdf.name);

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
    int              status;
    int              var_id;

    extern nc_struct netcdf;

    if (!netcdf.name) {
        netcdf.name = nc_name;
        status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
        check_nc_status(status, "Error opening %s", netcdf.name);
    }
    else {
        // check if another netcdf has to be opened
        if (strcmp(netcdf.name, nc_name) != false) {
            // close the old netcdf file
            status = nc_close(netcdf.id);
            check_nc_status(status, "Error closing %s", netcdf.name);
            netcdf.id = 0;
            netcdf.name = '\0';

            // open the new netcdf file
            status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
            check_nc_status(status, "Error opening %s", nc_name);

            // set name
            netcdf.name = nc_name;
        }
    }


    /* get NetCDF variable */
    status = nc_inq_varid(netcdf.id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    netcdf.name);

    status = nc_get_vara_float(netcdf.id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    netcdf.name);

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
    int              status;
    int              var_id;

    extern nc_struct netcdf;

    if (!netcdf.name) {
        netcdf.name = nc_name;
        // open netcdf file
        status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
        check_nc_status(status, "Error opening %s", nc_name);
    }
    else {
        // check if another netcdf has to be opened
        if (strcmp(netcdf.name, nc_name) != false) {
            // close the old netcdf file
            status = nc_close(netcdf.id);
            check_nc_status(status, "Error closing %s", netcdf.name);
            netcdf.id = 0;
            netcdf.name = '\0';

            // open the new netcdf file
            status = nc_open(nc_name, NC_NOWRITE, &(netcdf.id));
            check_nc_status(status, "Error opening %s", nc_name);

            // set name
            netcdf.name = nc_name;
        }
    }

    /* get NetCDF variable */
    status = nc_inq_varid(netcdf.id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    netcdf.name);

    status = nc_get_vara_int(netcdf.id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    netcdf.name);

    return status;
}
