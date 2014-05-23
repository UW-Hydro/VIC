#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(get_global_domain): %s\n", \
                        nc_strerror(e)); }

size_t
get_global_domain(char          *nc_name,
                  domain_struct *global_domain)
{
    int             *run = NULL;
    double          *var = NULL;
    size_t           i;
    size_t           j;
    size_t           x;
    size_t           y;
    size_t           d2count[2];
    size_t           d2start[2];
    int              nc_id;
    int              status;
    int              var_id;
    location_struct *location;

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
        nrerror("Memory allocation error in get_global_domain().");
    }

    get_nc_field_int(nc_name, "run_cell", d2start, d2count, run);

    for (x = 0, i = 0; x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
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
        nrerror("Memory allocation error in get_global_domain().");
    }
    for (i = 0; i < global_domain->ncells_global; i++) {
        initialize_location(&(global_domain->locations[i]));
    }

    // allocate memory for variables
    var = (double *) malloc(global_domain->n_ny * global_domain->n_nx *
                            sizeof(double));
    if (var == NULL) {
        nrerror("Memory allocation error in get_global_domain().");
    }

    // get longitude
    // TBD: read var id from file
    get_nc_field_double(nc_name, "xc", d2start, d2count, var);

    for (x = 0, i = 0, j = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (run[i]) {
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
    get_nc_field_double(nc_name, "yc", d2start, d2count, var);

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (run[i]) {
                location->latitude = var[i];
                location++;
            }
        }
    }

    // get area
    // TBD: read var id from file
    get_nc_field_double(nc_name, "area", d2start, d2count, var);

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (run[i]) {
                location->area = var[i];
                location++;
            }
        }
    }

    // get fraction
    // TBD: read var id from file
    get_nc_field_double(nc_name, "frac", d2start, d2count, var);

    for (x = 0, i = 0, location = global_domain->locations;
         x < global_domain->n_nx; x++) {
        for (y = 0; y < global_domain->n_ny; y++, i++) {
            if (run[i]) {
                location->frac = var[i];
                location++;
            }
        }
    }

    // free memory
    free(var);
    free(run);

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

void
print_domain(domain_struct *domain,
             bool           print_loc)
{
    int i;

    printf("domain:\n");
    printf("\tncells_global: %zd\n", domain->ncells_global);
    printf("\tn_nx: %zd\n", domain->n_nx);
    printf("\tn_ny: %zd\n", domain->n_ny);
    printf("\tncells_local: %zd\n", domain->ncells_local);
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
    printf("%zd: (%zd, %zd) {%zd: (%zd, %zd)}\n",
           location->global_cell_idx, location->global_x_idx,
           location->global_y_idx, location->local_cell_idx,
           location->local_x_idx, location->local_y_idx);
    printf("\t(%lf, %lf): %lf %lf\n",
           location->longitude, location->latitude,
           location->area, location->frac);
}
