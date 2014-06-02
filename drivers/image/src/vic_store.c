#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(vic_store): %s\n", nc_strerror(e)); }

#define MAXDIMS 5

void
vic_store(void)
{
    extern size_t              current;
    extern all_vars_struct    *all_vars;
    extern dmy_struct         *dmy;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern domain_struct       global_domain;
    extern lake_con_struct     lake_con;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    int                        status;
    int                        dimids[MAXDIMS];
    int                        v;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     n;
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

    nc_file_struct             nc_state_file;

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

    // create netcdf file for storing model state - keep the file open
    initialize_state_file(&nc_state_file);

    // set missing values
    for (i = 0; i < grid_size; i++) {
        cvar[i] = nc_state_file.c_fillvalue;
        ivar[i] = nc_state_file.i_fillvalue;
        dvar[i] = nc_state_file.d_fillvalue;
        fvar[i] = nc_state_file.f_fillvalue;
    }

    // initialize dimids to invalid values
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
    }

    // soil thermal node deltas
    ndims = 3;
    dimids[0] = nc_state_file.node_dimid;
    dimids[1] = nc_state_file.nj_dimid;
    dimids[2] = nc_state_file.ni_dimid;
    for (j = 0; j < options.Nnode; j++) {
        d3start[0] = j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            dvar[idx[i]] = (double) soil_con[i].dz_node[j];
        }
        put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                            &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                            dimids, ndims, "dz_node", d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            dvar[idx[i]] = nc_state_file.d_fillvalue;
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }


    // soil thermal node depths
    ndims = 3;
    dimids[0] = nc_state_file.node_dimid;
    dimids[1] = nc_state_file.nj_dimid;
    dimids[2] = nc_state_file.ni_dimid;
    for (j = 0; j < options.Nnode; j++) {
        d3start[0] = j;
        for (i = 0; i < global_domain.ncells_global; i++) {
            dvar[idx[i]] = (double) soil_con[i].Zsum_node[j];
        }
        put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                            &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                            dimids, ndims, "Zsum_node", d3start, d3count, dvar);
        for (i = 0; i < global_domain.ncells_global; i++) {
            dvar[idx[i]] = nc_state_file.d_fillvalue;
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }
    
    // total soil moisture
    ndims = 5;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.layer_dimid;
    dimids[3] = nc_state_file.nj_dimid;
    dimids[4] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d5start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d5start[2] = j;
                for (i = 0; i < global_domain.ncells_global; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[idx[i]] = (double) 
                                       all_vars[i].cell[v][m].layer[j].moist;
                    } else {
                        dvar[idx[i]] = nc_state_file.d_fillvalue;
                    }
                                            
                }
                put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                                    &(nc_state_file.nc_id), 
                                    nc_state_file.d_fillvalue,
                                    dimids, ndims, "Soil_moisture", 
                                    d5start, d5count, dvar);
                for (i = 0; i < global_domain.ncells_global; i++) {
                    dvar[idx[i]] = nc_state_file.d_fillvalue;
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }
    

    // close the netcdf file if it is still open
    if (nc_state_file.open == true) {
        status = nc_close(nc_state_file.nc_id);
        if (status != NC_NOERR) {
            ERR(status);
        }
    }

    free(idx);
    free(cvar);
    free(ivar);
    free(dvar);
    free(fvar);
}

void
initialize_state_file(nc_file_struct *nc)
{
    extern size_t           current;
    extern dmy_struct      *dmy;
    extern filenames_struct filenames;
    extern domain_struct    global_domain;
    extern option_struct    options;

    int                     status;
    int                     old_fill_mode;

    sprintf(nc->fname, "%s.%04d%02d%02d_%02d.nc",
            filenames.statefile, dmy[current].year, dmy[current].month,
            dmy[current].day, dmy[current].hour);

    nc->c_fillvalue = NC_FILL_CHAR;
    nc->i_fillvalue = NC_FILL_INT;
    nc->d_fillvalue = NC_FILL_DOUBLE;
    nc->f_fillvalue = NC_FILL_FLOAT;

    nc->band_size = options.SNOW_BAND;
    nc->frost_size = options.Nfrost;
    nc->layer_size = options.Nlayer;
    nc->ni_size = global_domain.n_nx;
    nc->nj_size = global_domain.n_ny;
    nc->node_size = options.Nnode;
    nc->root_zone_size = options.ROOT_ZONES;
    nc->veg_size = options.NVEGTYPES;

    // open the netcdf file
    status = nc_create(nc->fname, NC_NETCDF4 | NC_CLASSIC_MODEL, &(nc->nc_id));
    if (status != NC_NOERR) {
        ERR(status);
    }
    nc->open = true;

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim snow_band\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim frost_area\n");
        ERR(status);
    }
    
    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim nlayer\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "ni", nc->ni_size, &(nc->ni_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim ni\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "nj", nc->nj_size, &(nc->nj_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim nj\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim node\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim root_zone\n");
        ERR(status);
    }

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim veg_class\n");
        ERR(status);
    }

    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // status = nc_close(nc->nc_id);
    // if (status != NC_NOERR) {
    // ERR(status);
    // }
    // nc->open = false;
}

#undef MAXDIMS
