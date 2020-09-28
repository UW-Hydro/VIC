/******************************************************************************
 * @section DESCRIPTION
 *
 * Get global domain data from file.
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
 * @brief    Get global domain information.
 *****************************************************************************/
size_t
get_global_domain(nameid_struct *domain_nc_nameid,
                  nameid_struct *param_nc_nameid,
                  domain_struct *global_domain)
{
    int    *run = NULL;
    int    *mask = NULL;
    int     typeid;
    double *var = NULL;
    size_t  i;
    size_t  j;
    size_t  d2count[2];
    size_t  d2start[2];

    global_domain->n_nx = get_nc_dimension(domain_nc_nameid,
                                           global_domain->info.x_dim);
    global_domain->n_ny = get_nc_dimension(domain_nc_nameid,
                                           global_domain->info.y_dim);

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain->n_ny;
    d2count[1] = global_domain->n_nx;

    // get total number of gridcells in domain
    global_domain->ncells_total = global_domain->n_ny * global_domain->n_nx;

    // allocate memory for mask and cells to be run
    mask = malloc(global_domain->ncells_total * sizeof(*mask));
    check_alloc_status(mask, "Memory allocation error.");
    run = malloc(global_domain->ncells_total * sizeof(*run));
    check_alloc_status(run, "Memory allocation error.");

    // Get mask variable from the domain file
    // (check whether mask variable is int type)
    typeid = get_nc_var_type(domain_nc_nameid, global_domain->info.mask_var);
    if (typeid != NC_INT) {
        log_err("Mask variable in the domain file must be integer type.");
    }
    get_nc_field_int(domain_nc_nameid, global_domain->info.mask_var, d2start,
                     d2count,
                     mask);

    // Get run_cell variable from the parameter file
    // (check whether run_cell variable is int type)
    typeid = get_nc_var_type(param_nc_nameid, "run_cell");
    if (typeid != NC_INT) {
        log_err("Run_cell variable in the parameter file must be integer type.");
    }
    get_nc_field_int(param_nc_nameid, "run_cell", d2start, d2count,
                     run);

    // Check whether cells with run_cell == 1 are all within the mask domain
    for (i = 0; i < global_domain->ncells_total; i++) {
        if (run[i] == 1 && mask[i] != 1) {
            log_err("Run_cell = 1 should only appear within the mask of the "
                    "domain file.");
        }
    }

    // Store active cell information into variables
    for (i = 0; i < global_domain->ncells_total; i++) {
        if (run[i] == 1) {
            global_domain->ncells_active++;
        }
    }
    debug("%zu active grid cells found in run_cell in the parameter file.",
          global_domain->ncells_active);

    global_domain->locations =
        malloc(global_domain->ncells_total * sizeof(*global_domain->locations));
    check_alloc_status(global_domain->locations, "Memory allocation error.");
    for (i = 0; i < global_domain->ncells_total; i++) {
        initialize_location(&(global_domain->locations[i]));
    }

    for (i = 0; i < global_domain->ncells_total; i++) {
        if (run[i] == 1) {
            global_domain->locations[i].run = true;
        }
    }

    for (i = 0, j = 0; i < global_domain->ncells_total; i++) {
        if (run[i] == 1) {
            global_domain->locations[i].io_idx = i;
            global_domain->locations[i].global_idx = j;
            j++;
        }
    }

    // allocate memory for variables
    var = malloc(global_domain->ncells_total * sizeof(*var));
    check_alloc_status(var, "Memory allocation error.");

    // get area
    // TBD: read var id from file
    get_nc_field_double(domain_nc_nameid, global_domain->info.area_var,
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].area = var[i];
    }

    // get fraction
    // TBD: read var id from file
    get_nc_field_double(domain_nc_nameid, global_domain->info.frac_var,
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].frac = var[i];
    }

    // get lat and lon coordinates
    get_nc_latlon(domain_nc_nameid, global_domain);

    // check whether lat and lon coordinates in the parameter file match those
    // in the domain file
    compare_ncdomain_with_global_domain(param_nc_nameid);

    // free memory
    free(var);
    free(run);
    free(mask);

    return global_domain->ncells_active;
}

/******************************************************************************
 * @brief    Get lat and lon coordinates information from a netCDF file and
             store in nc_domain structure
 *****************************************************************************/
void
get_nc_latlon(nameid_struct *nc_nameid,
              domain_struct *nc_domain)
{
    double *var = NULL;
    double *var_lon = NULL;
    double *var_lat = NULL;
    size_t  i;
    size_t  j;
    size_t  d2count[2];
    size_t  d2start[2];
    size_t  d1count[1];
    size_t  d1start[1];


    nc_domain->n_nx = get_nc_dimension(nc_nameid,
                                       nc_domain->info.x_dim);
    nc_domain->n_ny = get_nc_dimension(nc_nameid,
                                       nc_domain->info.y_dim);

    // Get number of lat/lon dimensions.
    nc_domain->info.n_coord_dims = get_nc_varndimensions(nc_nameid,
                                                         nc_domain->info.lon_var);
    if (nc_domain->info.n_coord_dims !=
        (size_t) get_nc_varndimensions(nc_nameid, nc_domain->info.lat_var)) {
        log_err("Un even number of dimensions for %s and %s in: %s",
                nc_domain->info.lon_var, nc_domain->info.lat_var,
                nc_nameid->nc_filename);
    }

    if (nc_domain->info.n_coord_dims == 1) {
        // allocate memory for variables
        var_lon = malloc(nc_domain->n_nx * sizeof(*var_lon));
        check_alloc_status(var_lon, "Memory allocation error.");
        var_lat = malloc(nc_domain->n_ny * sizeof(*var_lat));
        check_alloc_status(var_lat, "Memory allocation error.");


        d1start[0] = 0;
        d1count[0] = nc_domain->n_nx;

        // get longitude for unmasked grid
        get_nc_field_double(nc_nameid, nc_domain->info.lon_var,
                            d1start, d1count, var_lon);
        for (j = 0; j < nc_domain->n_ny; j++) {
            for (i = 0; i < nc_domain->n_nx; i++) {
                // rescale to [-180., 180]. Note that the if statement is not strictly
                // needed, but it prevents -180 from turning into 180 and vice versa
                if (var_lon[i] < -180.f || var_lon[i] > 180.f) {
                    var_lon[i] -= round(var_lon[i] / 360.f) * 360.f;
                }
                nc_domain->locations[j * nc_domain->n_nx + i].longitude =
                    (double) var_lon[i];
            }
        }

        d1start[0] = 0;
        d1count[0] = nc_domain->n_ny;

        // get latitude for unmasked grid
        get_nc_field_double(nc_nameid, nc_domain->info.lat_var,
                            d1start, d1count, var_lat);
        for (i = 0; i < nc_domain->n_ny; i++) {
            for (j = 0; j < nc_domain->n_nx; j++) {
                nc_domain->locations[i * nc_domain->n_nx + j].latitude =
                    (double) var_lat[i];
            }
        }

        // free memory
        free(var_lon);
        free(var_lat);
    }
    else if (nc_domain->info.n_coord_dims == 2) {
        // allocate memory for variables
        var = malloc(nc_domain->ncells_total * sizeof(*var));
        check_alloc_status(var, "Memory allocation error.");


        d2start[0] = 0;
        d2start[1] = 0;
        d2count[0] = nc_domain->n_ny;
        d2count[1] = nc_domain->n_nx;

        // get longitude for unmasked grid
        get_nc_field_double(nc_nameid, nc_domain->info.lon_var,
                            d2start, d2count, var);
        for (i = 0; i < nc_domain->ncells_total; i++) {
            // rescale to [-180., 180]. Note that the if statement is not strictly
            // needed, but it prevents -180 from turning into 180 and vice versa
            if (var[i] < -180.f || var[i] > 180.f) {
                var[i] -= round(var[i] / 360.f) * 360.f;
            }
            nc_domain->locations[i].longitude = var[i];
        }

        // get latitude for unmasked grid
        get_nc_field_double(nc_nameid, nc_domain->info.lat_var,
                            d2start, d2count, var);
        for (i = 0; i < nc_domain->ncells_total; i++) {
            nc_domain->locations[i].latitude = var[i];
        }
        // free memory
        free(var);
    }
    else {
        log_err("Number of dimensions for %s and %s should be 1 or 2 in: %s",
                nc_domain->info.lon_var, nc_domain->info.lat_var,
                nc_nameid->nc_filename);
    }
}

/******************************************************************************
 * @brief    Copy domain info from one domain structure to another
 *****************************************************************************/
void
copy_domain_info(domain_struct *domain_from,
                 domain_struct *domain_to)
{
    strcpy(domain_to->info.x_dim, domain_from->info.x_dim);
    strcpy(domain_to->info.y_dim, domain_from->info.y_dim);

    strcpy(domain_to->info.lon_var, domain_from->info.lon_var);
    strcpy(domain_to->info.lat_var, domain_from->info.lat_var);

    domain_to->n_nx = domain_from->n_nx;
    domain_to->n_ny = domain_from->n_ny;

    domain_to->ncells_total = domain_from->ncells_total;
}

/******************************************************************************
 * @brief    Initialize domain structure.
 *****************************************************************************/
void
initialize_domain(domain_struct *domain)
{
    domain->ncells_total = 0;
    domain->ncells_active = 0;
    domain->n_nx = MISSING_USI;
    domain->n_ny = MISSING_USI;
    domain->locations = NULL;

    // Initialize domain info structure
    strcpy(domain->info.lat_var, "MISSING");
    strcpy(domain->info.lon_var, "MISSING");
    strcpy(domain->info.mask_var, "MISSING");
    strcpy(domain->info.area_var, "MISSING");
    strcpy(domain->info.frac_var, "MISSING");
    strcpy(domain->info.y_dim, "MISSING");
    strcpy(domain->info.x_dim, "MISSING");
    domain->info.n_coord_dims = MISSING_USI;
}

/******************************************************************************
 * @brief    Initialize location structure.
 *****************************************************************************/
void
initialize_location(location_struct *location)
{
    location->run = false;
    location->latitude = MISSING;
    location->longitude = MISSING;
    location->area = MISSING;
    location->frac = MISSING;
    location->nveg = MISSING_USI;
    location->global_idx = MISSING_USI;
    location->io_idx = MISSING_USI;
    location->local_idx = MISSING_USI;
}

/******************************************************************************
 * @brief    Read the number of vegetation type per grid cell from file
 *****************************************************************************/
void
add_nveg_to_global_domain(nameid_struct *nc_nameid,
                          domain_struct *global_domain)
{
    size_t d2count[2];
    size_t d2start[2];
    size_t i;
    int   *ivar = NULL;

    ivar = malloc(global_domain->ncells_total * sizeof(*ivar));
    check_alloc_status(ivar, "Memory allocation error.");

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain->n_ny;
    d2count[1] = global_domain->n_nx;
    get_nc_field_int(nc_nameid, "Nveg", d2start, d2count, ivar);

    for (i = 0; i < global_domain->ncells_total; i++) {
        global_domain->locations[i].nveg = (size_t) ivar[i];
    }

    free(ivar);
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
