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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Get global domain information.
 *****************************************************************************/
size_t
get_global_domain(char          *nc_name,
                  domain_struct *global_domain)
{
    int    *run = NULL;
    double *var = NULL;
    size_t  i;
    size_t  j;
    size_t  x;
    size_t  y;
    size_t *idx;
    size_t  d2count[2];
    size_t  d2start[2];

    initialize_domain(global_domain);

    global_domain->n_nx = get_nc_dimension(nc_name, "ni");
    global_domain->n_ny = get_nc_dimension(nc_name, "nj");

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain->n_ny;
    d2count[1] = global_domain->n_nx;

    // allocate memory for cells to be run
    run = (int *) malloc(global_domain->n_ny * global_domain->n_nx *
                         sizeof(int));
    if (run == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }

    get_nc_field_int(nc_name, "run_cell", d2start, d2count, run);

    for (y = 0, i = 0; y < global_domain->n_ny; y++) {
        for (x = 0; x < global_domain->n_nx; x++, i++) {
            if (run[i]) {
                global_domain->ncells_global++;
            }
        }
    }

    // if MASTER_PROC
    global_domain->locations = (location_struct *)
                               malloc(global_domain->ncells_global *
                                      sizeof(location_struct));
    if (global_domain->locations == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }
    for (i = 0; i < global_domain->ncells_global; i++) {
        initialize_location(&(global_domain->locations[i]));
    }

    // allocate memory for variables
    var = (double *) malloc(global_domain->n_ny * global_domain->n_nx *
                            sizeof(double));
    if (var == NULL) {
        log_err("Memory allocation error in get_global_domain().");
    }

    for (y = 0, i = 0, j = 0; y < global_domain->n_ny; y++) {
        for (x = 0; x < global_domain->n_nx; x++, i++) {
            if (run[i]) {
                global_domain->locations[j].global_cell_idx = j;
                global_domain->locations[j].global_x_idx = x;
                global_domain->locations[j].global_y_idx = y;
                j++;
            }
        }
    }

    // get 1D indices used in mapping the netcdf fields to the locations
    idx = (size_t *) malloc(global_domain->ncells_global *
                            sizeof(size_t));
    if (idx == NULL) {
        log_err("Memory allocation error in vic_init().");
    }
    for (i = 0; i < global_domain->ncells_global; i++) {
        idx[i] = get_global_idx(global_domain, i);
    }

    // get longitude -
    // TBD: read var id from file
    get_nc_field_double(nc_name, "xc",
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_global; i++) {
        // rescale to [-180., 180]. Note that the if statement is not strictly
        // needed, but it prevents -180 from turning into 180 and vice versa
        if (var[idx[i]] < -180.f || var[idx[i]] > 180.f) {
            var[idx[i]] -= round(var[idx[i]] / 360.f) * 360.f;
        }
        global_domain->locations[i].longitude = (double) var[idx[i]];
    }

    // get latitude
    // TBD: read var id from file
    get_nc_field_double(nc_name, "yc",
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_global; i++) {
        global_domain->locations[i].latitude = (double) var[idx[i]];
    }

    // get area
    // TBD: read var id from file
    get_nc_field_double(nc_name, "area",
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_global; i++) {
        global_domain->locations[i].area = (double) var[idx[i]];
    }

    // get fraction
    // TBD: read var id from file
    get_nc_field_double(nc_name, "frac",
                        d2start, d2count, var);
    for (i = 0; i < global_domain->ncells_global; i++) {
        global_domain->locations[i].frac = (double) var[idx[i]];
    }

    // free memory
    free(idx);
    free(var);
    free(run);

    // print_domain(global_domain, true);

    return global_domain->ncells_global;
}

/******************************************************************************
 * @brief    Initialize domain structure.
 *****************************************************************************/
void
initialize_domain(domain_struct *domain)
{
    domain->ncells_global = 0;
    domain->n_nx = 0;
    domain->n_ny = 0;
    domain->ncells_local = 0;
    domain->locations = NULL;
}

/******************************************************************************
 * @brief    Initialize location structure.
 *****************************************************************************/
void
initialize_location(location_struct *location)
{
    location->latitude = 0;
    location->longitude = 0;
    location->area = 0;
    location->frac = 0;
    location->global_cell_idx = 0;
    location->global_x_idx = 0;
    location->global_y_idx = 0;
    location->local_cell_idx = 0;
    location->local_x_idx = 0;
    location->local_y_idx = 0;
}

/******************************************************************************
 * @brief    Get global index,
 *****************************************************************************/
size_t
get_global_idx(domain_struct *domain,
               size_t         i)
{
    size_t idx = -1;

    if (i < domain->ncells_global) {
        idx = domain->locations[i].global_y_idx * domain->n_nx +
              domain->locations[i].global_x_idx;
    }

    return idx;
}
