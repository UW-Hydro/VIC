/******************************************************************************
 * @section DESCRIPTION
 *
 * Routines to compare the global domain to the other VIC input files, such as
 * parameter or state files.
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
 * @brief    Check that the cooridnates and dimensions in a netcdf file matches
             the global domain.
 *****************************************************************************/
void
compare_ncdomain_with_global_domain(nameid_struct *nc_nameid)
{
    extern domain_struct global_domain;

    domain_struct        ncfile_domain;

    size_t               i;

    // read the lat lon coordinates info from ncfile
    // (e.g. parameters file or state file)
    ncfile_domain.locations =
        malloc(global_domain.ncells_total *
               sizeof(*(ncfile_domain.locations)));
    check_alloc_status(ncfile_domain.locations, "Memory allocation error.");
    copy_domain_info(&global_domain, &ncfile_domain);
    get_nc_latlon(nc_nameid, &ncfile_domain);

    // using the ncfile_domain, we can compare the values to the global domain.

    // dimension shapes match (lat/lon)
    if (global_domain.n_nx != ncfile_domain.n_nx) {
        log_err("x dimension in parameters file does not match domain");
    }
    if (global_domain.n_ny != ncfile_domain.n_ny) {
        log_err("y dimension in parameters file does not match domain");
    }

    // loop over all grid cells and check that the two domains are identical
    for (i = 0; i < global_domain.ncells_total; i++) {
        // latitude matches
        if (!assert_close_double(ncfile_domain.locations[i].latitude,
                                 global_domain.locations[i].latitude,
                                 0, 0.01)) {
            log_err("latitude in parameter (%lf) file does not match the "
                    "latitude in the domain file (%lf) for gridcell %zu",
                    ncfile_domain.locations[i].latitude,
                    global_domain.locations[i].latitude, i);
        }
        // longitude matches
        if (!assert_close_double(ncfile_domain.locations[i].longitude,
                                 global_domain.locations[i].longitude,
                                 0, 0.01)) {
            log_err("longitude in parameter (%lf) file does not match the "
                    "longitude in the domain file (%lf) for gridcell %zu",
                    ncfile_domain.locations[i].longitude,
                    global_domain.locations[i].longitude, i);
        }
    }
    free(ncfile_domain.locations);
}
