/******************************************************************************
 * @section DESCRIPTION
 *
 * Get global domain data from file.
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
 * @brief    Get global domain information.
 *****************************************************************************/
size_t
get_global_domain(char          *nc_name,
                  domain_struct *global_domain)
{
    int    *run = NULL;
    double *var = NULL;
    double *var_lon = NULL;
    double *var_lat = NULL;
    size_t  i;
    size_t  j;
    size_t  d2count[2];
    size_t  d2start[2];
    size_t  d1count[1];
    size_t  d1start[1];

    global_domain->n_nx = get_nc_dimension(nc_name, global_domain->info.x_dim);
    global_domain->n_ny = get_nc_dimension(nc_name, global_domain->info.y_dim);

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain->n_ny;
    d2count[1] = global_domain->n_nx;

    // get total number of gridcells in domain
    global_domain->ncells_total = global_domain->n_ny * global_domain->n_nx;

    // allocate memory for cells to be run
    run = malloc(global_domain->ncells_total * sizeof(*run));
    if (run == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }

    get_nc_field_int(nc_name, global_domain->info.mask_var, d2start, d2count,
                     run);

    for (i = 0; i < global_domain->ncells_total; i++) {
        if (run[i]) {
            global_domain->ncells_active++;
        }
    }

    // if MASTER_PROC
    global_domain->locations =
        malloc(global_domain->ncells_total * sizeof(*global_domain->locations));
    if (global_domain->locations == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }
    for (i = 0; i < global_domain->ncells_total; i++) {
        initialize_location(&(global_domain->locations[i]));
    }

    // allocate memory for variables
    var = malloc(global_domain->ncells_total * sizeof(*var));
    if (var == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }

    for (i = 0; i < global_domain->ncells_total; i++) {
        if (run[i]) {
            global_domain->locations[i].run = true;
        }
    }

    for (i = 0, j = 0; i < global_domain->ncells_total; i++) {
        if (run[i]) {
            global_domain->locations[i].io_idx = i;
            global_domain->locations[i].global_idx = j;
            j++;
        }
    }

    // Get number of lat/lon dimensions.
    global_domain->info.n_coord_dims = get_nc_varndimensions(nc_name,
                                                             global_domain->info.lon_var);
    if (global_domain->info.n_coord_dims !=
        (size_t) get_nc_varndimensions(nc_name, global_domain->info.lat_var)) {
        log_err("Un even number of dimensions for %s and %s in: %s",
                global_domain->info.lon_var, global_domain->info.lat_var,
                nc_name);
    }

    if (global_domain->info.n_coord_dims == 1) {
        // allocate memory for variables
        var_lon = malloc(global_domain->n_nx * sizeof(*var_lon));
        if (var_lon == NULL) {
            log_err("Memory allocation error in get_global_domain().");
        }
        var_lat = malloc(global_domain->n_ny * sizeof(*var_lat));
        if (var_lat == NULL) {
            log_err("Memory allocation error in get_global_domain().");
        }

        d1start[0] = 0;
        d1count[0] = global_domain->n_nx;

        // get longitude for unmasked grid
        get_nc_field_double(nc_name, global_domain->info.lon_var,
                            d1start, d1count, var_lon);
        for (i = 0; i < global_domain->n_nx; i++) {
            // rescale to [-180., 180]. Note that the if statement is not strictly
            // needed, but it prevents -180 from turning into 180 and vice versa
            if (var_lon[i] < -180.f || var_lon[i] > 180.f) {
                var_lon[i] -= round(var_lon[i] / 360.f) * 360.f;
            }
            global_domain->locations[i].longitude = (double) var_lon[i];
        }

        d1start[0] = 0;
        d1count[0] = global_domain->n_ny;

        // get latitude for unmasked grid
        get_nc_field_double(nc_name, global_domain->info.lat_var,
                            d1start, d1count, var_lat);
        for (i = 0; i < global_domain->n_ny; i++) {
            global_domain->locations[i].latitude = (double) var_lat[i];
        }

        free(var_lon);
        free(var_lat);
    }
    else if (global_domain->info.n_coord_dims == 2) {
        // get longitude for unmasked grid
        get_nc_field_double(nc_name, global_domain->info.lon_var,
                            d2start, d2count, var);
        for (i = 0; i < global_domain->ncells_total; i++) {
            // rescale to [-180., 180]. Note that the if statement is not strictly
            // needed, but it prevents -180 from turning into 180 and vice versa
            if (var[i] < -180.f || var[i] > 180.f) {
                var[i] -= round(var[i] / 360.f) * 360.f;
            }
            global_domain->locations[i].longitude = (double) var[i];
        }

        // get latitude for unmasked grid
        get_nc_field_double(nc_name, global_domain->info.lat_var,
                            d2start, d2count, var);
        for (i = 0; i < global_domain->ncells_total; i++) {
            global_domain->locations[i].latitude = (double) var[i];
        }
    }
    else {
        log_err("Number of dimensions for %s and %s should be 1 or 2 in: %s",
                global_domain->info.lon_var, global_domain->info.lat_var,
                nc_name);
    }

    // get area
    // TBD: read var id from file
    get_nc_field_double(nc_name, global_domain->info.area_var,
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].area = (double) var[i];
    }

    // get fraction
    // TBD: read var id from file
    get_nc_field_double(nc_name, global_domain->info.frac_var,
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].frac = (double) var[i];
    }

    // free memory
    free(var);
    free(run);

    return global_domain->ncells_active;
}

/******************************************************************************
 * @brief    Initialize domain structure.
 *****************************************************************************/
void
initialize_domain(domain_struct *domain)
{
    domain->ncells_total = 0;
    domain->ncells_active = 0;
    domain->n_nx = 0;
    domain->n_ny = 0;
    domain->locations = NULL;

    // Initialize domain info structure
    strcpy(domain->info.lat_var, "MISSING");
    strcpy(domain->info.lon_var, "MISSING");
    strcpy(domain->info.mask_var, "MISSING");
    strcpy(domain->info.area_var, "MISSING");
    strcpy(domain->info.frac_var, "MISSING");
    strcpy(domain->info.y_dim, "MISSING");
    strcpy(domain->info.x_dim, "MISSING");
    domain->info.n_coord_dims = 0;
}

/******************************************************************************
 * @brief    Initialize location structure.
 *****************************************************************************/
void
initialize_location(location_struct *location)
{
    location->run = 0;
    location->latitude = 0;
    location->longitude = 0;
    location->area = 0;
    location->frac = 0;
    location->nveg = 0;
    location->global_idx = 0;
    location->io_idx = 0;
    location->local_idx = 0;
}

/******************************************************************************
 * @brief    Read the number of vegetation type per grid cell from file
 *****************************************************************************/
void
add_nveg_to_global_domain(char          *nc_name,
                          domain_struct *global_domain)
{
    size_t  d2count[2];
    size_t  d2start[2];
    size_t  i;
    double *dvar = NULL;

    dvar = malloc(global_domain->ncells_total * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in add_nveg_to_global_domain().");
    }

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain->n_ny;
    d2count[1] = global_domain->n_nx;
    get_nc_field_double(nc_name, "Nveg", d2start, d2count, dvar);

    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].nveg = (size_t) dvar[i];
    }

    free(dvar);
}

/******************************************************************************
 * @brief    Parse the domain variable types.
 *****************************************************************************/
void
get_domain_type(char *cmdstr)
{
    extern domain_struct global_domain;

    char                 optstr[MAXSTRING];
    char                 ncvarname[MAXSTRING];

    strcpy(ncvarname, "MISSING");

    sscanf(cmdstr, "%*s %s %s", optstr, ncvarname);

    // Lattitude variable name
    if (strcasecmp("LAT", optstr) == 0) {
        strcpy(global_domain.info.lat_var, ncvarname);
    }
    // Longitude variable name
    else if (strcasecmp("LON", optstr) == 0) {
        strcpy(global_domain.info.lon_var, ncvarname);
    }
    // Mask variable name
    else if (strcasecmp("MASK", optstr) == 0) {
        strcpy(global_domain.info.mask_var, ncvarname);
    }
    // Area variable name
    else if (strcasecmp("AREA", optstr) == 0) {
        strcpy(global_domain.info.area_var, ncvarname);
    }
    // Fraction variable name
    else if (strcasecmp("FRAC", optstr) == 0) {
        strcpy(global_domain.info.frac_var, ncvarname);
    }
    // y dimension name
    else if (strcasecmp("YDIM", optstr) == 0) {
        strcpy(global_domain.info.y_dim, ncvarname);
    }
    // x dimension name
    else if (strcasecmp("XDIM", optstr) == 0) {
        strcpy(global_domain.info.x_dim, ncvarname);
    }
    else {
        log_err("Unrecognized domain variable: %s %s", optstr, ncvarname);
    }
}
