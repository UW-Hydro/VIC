#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(get_global_domain): %s\n", \
                        nc_strerror(e)); }

long int
get_global_domain(char          *ncName,
                  domain_struct *global_domain)
{
    double          *var;
    int              i;
    int              j;
    int             *mask = NULL;
    int              ncid;
    int              nx_dimid;
    int              ny_dimid;
    int              status;
    int              varid;
    int              x;
    int              y;
    location_struct *location;

    initialize_domain(global_domain);

    // open the netcdf file
    status = nc_open(ncName, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get netcdf id for x dimension
    // TBD: read dimension id from file
    status = nc_inq_dimid(ncid, "ni", &nx_dimid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get size of x dimension
    status = nc_inq_dimlen(ncid, nx_dimid, &(global_domain->n_nx));
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get netcdf id for y dimension
    // TBD: read dimension id from file
    status = nc_inq_dimid(ncid, "nj", &ny_dimid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get size of y dimension
    status = nc_inq_dimlen(ncid, ny_dimid, &(global_domain->n_ny));
    if (status != NC_NOERR) {
        ERR(status);
    }

    // allocate memory for mask
    mask = (int *) malloc((size_t) (global_domain->n_ny * global_domain->n_nx) *
                          sizeof(int));
    if (mask == NULL) {
        nrerror("Memory allocation error in get_global_domain().");
    }

    /* get NetCDF variable */
    // TBD: read var id from file
    status = nc_inq_varid(ncid, "mask", &varid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_var_int(ncid, varid, mask);
    if (status != NC_NOERR) {
        ERR(status);
    }

    for (x = 0, i = 0; x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (mask[i]) {
                global_domain->ncells_global++;
            }
        }
    }

    // if MASTER_PROC
    global_domain->locations = (location_struct *)
                               malloc((size_t) global_domain->ncells_global *
                                      sizeof(location_struct));
    if (global_domain->locations == NULL) {
        nrerror("Memory allocation error in get_global_domain().");
    }
    for (i = 0; i < global_domain->ncells_global; i++) {
        initialize_location(&(global_domain->locations[i]));
    }

    // allocate memory for variables
    var = (double *) malloc((size_t) (global_domain->n_ny *
                                      global_domain->n_nx) *
                            sizeof(double));
    if (var == NULL) {
        nrerror("Memory allocation error in get_global_domain().");
    }

    // get longitude
    // TBD: read var id from file
    status = nc_inq_varid(ncid, "xc", &varid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_var_double(ncid, varid, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    for (x = 0, i = 0, j = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (mask[i]) {
                location->longitude = var[i];
                location->global_cell_idx = j;
                location->global_x_idx = x;
                location->global_y_idx = y;
                j++;
                location++;
            }
        }
    }

    // get latitude
    // TBD: read var id from file
    status = nc_inq_varid(ncid, "yc", &varid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_var_double(ncid, varid, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (mask[i]) {
                location->latitude = var[i];
                location++;
            }
        }
    }

    // get area
    // TBD: read var id from file
    status = nc_inq_varid(ncid, "area", &varid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_var_double(ncid, varid, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (mask[i]) {
                location->area = var[i];
                location++;
            }
        }
    }

    // get fraction
    // TBD: read var id from file
    status = nc_inq_varid(ncid, "frac", &varid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_var_double(ncid, varid, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (mask[i]) {
                location->frac = var[i];
                location++;
            }
        }
    }

    // close the netcdf file
    status = nc_close(ncid);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // free memory
    free(var);
    free(mask);

    // print_domain(global_domain, true);

    return global_domain->ncells_global;
}

void
initialize_domain(domain_struct *domain)
{
    domain->ncells_global = 0;
    domain->n_nx = 0;
    domain->n_ny = 0;
    domain->ncells_local = 0;
    domain->locations = NULL;
}

void
initialize_location(location_struct *location)
{
    location->latitude = -1;
    location->longitude = -1;
    location->area = -1;
    location->frac = -1;
    location->global_cell_idx = -1;
    location->global_x_idx = -1;
    location->global_y_idx = -1;
    location->local_cell_idx = -1;
    location->local_x_idx = -1;
    location->local_y_idx = -1;
}

void
print_domain(domain_struct *domain,
             bool           print_loc)
{
    int i;

    printf("domain:\n");
    printf("\tncells_global: %ld\n", domain->ncells_global);
    printf("\tn_nx: %zd\n", domain->n_nx);
    printf("\tn_ny: %zd\n", domain->n_ny);
    printf("\tncells_local: %ld\n", domain->ncells_local);
    printf("\tlocations: %p\n", domain->locations);
    if (print_loc) {
        for (i = 0; i < domain->ncells_global; i++) {
            print_location(&(domain->locations[i]));
        }
    }
}

void
print_location(location_struct *location)
{
    printf("%ld: (%ld, %ld) {%ld: (%ld, %ld)}\n",
           location->global_cell_idx, location->global_x_idx,
           location->global_y_idx, location->local_cell_idx,
           location->local_x_idx, location->local_y_idx);
    printf("\t(%lf, %lf): %lf %lf\n",
           location->longitude, location->latitude,
           location->area, location->frac);
}
