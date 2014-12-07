/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine handles the startup tasks.
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
 * @brief    Wrapper function for VIC startup tasks.
 *****************************************************************************/
void
vic_start(void)
{
    extern filenames_struct filenames;
    extern filep_struct     filep;
    extern domain_struct    global_domain;
    extern option_struct    options;

    // Initialize Log Destination
    initialize_log();

    // Initialize global structures
    initialize_options();
    initialize_global();
    initialize_parameters();
    initialize_filenames();

    // read global settings
    filep.globalparam = open_file(filenames.global, "r");
    get_global_param(filep.globalparam);

    // Set Log Destination
    setup_logging();

    // set model constants
    if (!strcasecmp(filenames.constants, "MISSING")) {
        filep.constants = open_file(filenames.constants, "r");
        get_parameters(filep.constants);
    }

    // read domain info
    get_global_domain(filenames.domain, &global_domain);

    // decompose the mask

    // get dimensions (number of vegetation types, soil zones, etc)
    options.ROOT_ZONES = get_nc_dimension(filenames.soil, "root_zone");
    options.Nlayer = get_nc_dimension(filenames.soil, "nlayer");
    options.NVEGTYPES = get_nc_dimension(filenames.veg, "veg_class");
    if (options.SNOW_BAND > 1) {
        if (options.SNOW_BAND !=
            get_nc_dimension(filenames.snowband, "snow_band")) {
            log_err("Number of snow bands in global file does not "
                    "match parameter file");
        }
    }
}
