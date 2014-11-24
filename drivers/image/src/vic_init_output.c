#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(vic_init_output): %s\n", \
                        nc_strerror(e)); }

void
vic_init_output(void)
{
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern domain_struct       global_domain;
    extern filep_struct        filep;
    extern global_param_struct global_param;
    extern nc_file_struct      nc_hist_file;
    extern nc_var_struct       nc_vars[N_OUTVAR_TYPES];
    extern lake_con_struct     lake_con;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    size_t                     i;

    // initialize the output data structures
    for (i = 0; i < global_domain.ncells_global; i++) {
        put_data(&(all_vars[i]), &(atmos[i]), &(soil_con[i]), veg_con[i],
                 veg_lib[i], &lake_con, out_data[i], &(save_data[i]),
                -global_param.nrecs);
    }

    // determine which variables will be written to the history file
    parse_output_info(filep.globalparam, out_data);

    // open the netcdf history file
    initialize_history_file(&nc_hist_file);

    // initialize netcdf info for output variables
    vic_nc_info(&nc_hist_file, out_data, nc_vars);
}

void
initialize_history_file(nc_file_struct *nc)
{
    extern filenames_struct filenames;
    extern domain_struct    global_domain;
    extern option_struct    options;

    int                     status;
    int                     old_fill_mode;

    sprintf(nc->fname, "%s", filenames.result_dir);

    nc->c_fillvalue = NC_FILL_CHAR;
    nc->i_fillvalue = NC_FILL_INT;
    nc->d_fillvalue = NC_FILL_DOUBLE;
    nc->f_fillvalue = NC_FILL_FLOAT;

    nc->band_size = options.SNOW_BAND;
    nc->front_size = MAX_FRONTS;
    nc->frost_size = options.Nfrost;
    nc->layer_size = options.Nlayer;
    nc->ni_size = global_domain.n_nx;
    nc->nj_size = global_domain.n_ny;
    nc->node_size = options.Nnode;
    nc->root_zone_size = options.ROOT_ZONES;
    nc->time_size = NC_UNLIMITED;
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

    status = nc_def_dim(nc->nc_id, "front", nc->front_size,
                        &(nc->front_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim front\n");
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

    status = nc_def_dim(nc->nc_id, "time", nc->time_size,
                        &(nc->time_dimid));
    if (status != NC_NOERR) {
        fprintf(stderr, "nc_def_dim time\n");
        ERR(status);
    }


    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }
}
