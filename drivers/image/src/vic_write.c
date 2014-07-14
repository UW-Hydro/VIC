#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_write(void)
{
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern nc_file_struct      nc_hist_file;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;

    int                        status;
    int                        dimids[MAXDIMS];
    int                        v;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     p;
    size_t                     grid_size;
    size_t                     ndims;
    char                      *cvar = NULL;
    int                       *ivar = NULL;
    double                    *dvar = NULL;
    float                     *fvar = NULL;
    size_t                    *idx = NULL;
    size_t                     d2count[2];
    size_t                     d2start[2];
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];
    size_t                     d5count[5];
    size_t                     d5start[5];
    size_t                     d6count[6];
    size_t                     d6start[6];

    grid_size = global_domain.n_ny * global_domain.n_nx;

    // allocate memory for variables to be stored
    cvar = (char *) malloc(grid_size * sizeof(char));
    if (cvar == NULL) {
        nrerror("Memory allocation error in vic_store().");
    }

    ivar = (int *) malloc(grid_size * sizeof(int));
    if (ivar == NULL) {
        nrerror("Memory allocation error in vic_store().");
    }

    dvar = (double *) malloc(grid_size * sizeof(double));
    if (dvar == NULL) {
        nrerror("Memory allocation error in vic_store().");
    }

    fvar = (float *) malloc(grid_size * sizeof(float));
    if (fvar == NULL) {
        nrerror("Memory allocation error in vic_store().");
    }

    // get 1D indices used in mapping the netcdf fields to the locations
    idx = (size_t *) malloc(global_domain.ncells_global *
                            sizeof(size_t));
    if (idx == NULL) {
        nrerror("Memory allocation error in vic_store().");
    }
    for (i = 0; i < global_domain.ncells_global; i++) {
        idx[i] = get_global_idx(&global_domain, i);
    }

    // initialize starts and counts
    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain.n_ny;
    d2count[1] = global_domain.n_nx;

    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;
    d4count[0] = 1;
    d4count[1] = 1;
    d4count[2] = global_domain.n_ny;
    d4count[3] = global_domain.n_nx;

    d5start[0] = 0;
    d5start[1] = 0;
    d5start[2] = 0;
    d5start[3] = 0;
    d5start[4] = 0;
    d5count[0] = 1;
    d5count[1] = 1;
    d5count[2] = 1;
    d5count[3] = global_domain.n_ny;
    d5count[4] = global_domain.n_nx;

    d6start[0] = 0;
    d6start[1] = 0;
    d6start[2] = 0;
    d6start[3] = 0;
    d6start[4] = 0;
    d6start[5] = 0;
    d6count[0] = 1;
    d6count[1] = 1;
    d6count[2] = 1;
    d6count[3] = 1;
    d6count[4] = global_domain.n_ny;
    d6count[5] = global_domain.n_nx;

    // set missing values
    for (i = 0; i < grid_size; i++) {
        cvar[i] = nc_hist_file.c_fillvalue;
        ivar[i] = nc_hist_file.i_fillvalue;
        dvar[i] = nc_hist_file.d_fillvalue;
        fvar[i] = nc_hist_file.f_fillvalue;
    }

    // initialize dimids to invalid values
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
    }    
}
