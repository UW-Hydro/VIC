/******************************************************************************
 * @section DESCRIPTION
 *
 * Put attribute to netCDF file
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
 * @brief    Put text attribute to netCDF.
 *****************************************************************************/
void
put_nc_attr(int         nc_id,
            int         var_id,
            const char *name,
            const char *value)
{
    int status;

    status = nc_put_att_text(nc_id, var_id, name, strlen(value), value);
    check_nc_status(status, "Error adding %s attribute in ncid %d", name,
                    nc_id);
}
