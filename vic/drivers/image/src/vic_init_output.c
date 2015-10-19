/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize output structures.
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
 * @brief    Initialzie output structures and determine which variables to
 *           write
 *****************************************************************************/
void
vic_init_output(void)
{
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern domain_struct       local_domain;
    extern filep_struct        filep;
    extern global_param_struct global_param;
    extern MPI_Datatype        mpi_nc_file_struct_type;
    extern int                 mpi_rank;
    extern nc_file_struct      nc_hist_file;
    extern nc_var_struct       nc_vars[N_OUTVAR_TYPES];
    extern lake_con_struct     lake_con;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    int                        status;
    size_t                     i;

    // initialize the output data structures
    for (i = 0; i < local_domain.ncells; i++) {
        put_data(&(all_vars[i]), &(atmos[i]), &(soil_con[i]), veg_con[i],
                 veg_lib[i], &lake_con, out_data[i], &(save_data[i]),
                 -global_param.nrecs);
    }

    if (mpi_rank == 0) {
        // determine which variables will be written to the history file
        parse_output_info(filep.globalparam, out_data);

        // open the netcdf history file
        initialize_history_file(&nc_hist_file);
    }

    // broadcast which variables to write.
    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        status = MPI_Bcast(&out_data[0][i].write, 1, MPI_C_BOOL,
                           0, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
            log_err("MPI error in vic_init_output(): %d\n", status);
        }
    }

    // broadcast history file info. Only the master process will write to it,
    // but the slave processes need some of the information to initialize as
    // well (particularly which variables to write and dimension sizes)
    status = MPI_Bcast(&nc_hist_file, 1, mpi_nc_file_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_init_output(): %d\n", status);
    }

    // initialize netcdf info for output variables
    vic_nc_info(&nc_hist_file, out_data, nc_vars);
}

/******************************************************************************
 * @brief    Initialize history files
 *****************************************************************************/
void
initialize_history_file(nc_file_struct *nc)
{
    extern filenames_struct    filenames;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern global_param_struct global_param;

    int                        status;
    int                        old_fill_mode;

    char                       str[100];
    char                       strUnit[6];

    size_t                     i;
    size_t                     ndims;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    int                        dimids[MAXDIMS];
    int                        time_var_id;
    int                        lon_var_id;
    int                        lat_var_id;
    double                    *dvar;

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
        log_err("Error creating %s", nc->fname);
    }
    nc->open = true;

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    if (status != NC_NOERR) {
        log_err("Error setting fill value in %s", nc->fname);
    }

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining snow_band dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "front", nc->front_size,
                        &(nc->front_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining front dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining frost_area dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining nlayer dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "lon", nc->ni_size, &(nc->ni_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining ni dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "lat", nc->nj_size, &(nc->nj_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining nj dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining node dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining root_zone dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining veg_class dimenension in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "time", nc->time_size,
                        &(nc->time_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining time dimenension in %s", nc->fname);
    }

    // define the netcdf variable time
    status = nc_def_var(nc->nc_id, "time", NC_INT, 1,
                        &(nc->time_dimid), &(time_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining time variable in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, time_var_id, "standard_name",
                             strlen("time"), "time");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    if (global_param.output_steps_per_day == 1) {
        sprintf(strUnit, "days");
    }
    else {
        sprintf(strUnit, "hours");
    }
    sprintf(str, "%s since %i-%02i-%02i 00:00:00", strUnit,
            global_param.startyear, global_param.startmonth,
            global_param.startday);

    status = nc_put_att_text(nc->nc_id, time_var_id, "units",
                             strlen(str), str);
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, time_var_id, "calendar",
                             strlen("standard"), "standard");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }

    ndims = options.COORD_DIMS_OUT;
    dstart[0] = 0;
    dstart[1] = 0;

    if (options.COORD_DIMS_OUT == 1) {
        dimids[0] = nc->ni_dimid;

        dcount[0] = nc->ni_size;
    }
    else if (options.COORD_DIMS_OUT == 2) {
        dimids[0] = nc->nj_dimid;
        dcount[0] = nc->nj_size;

        dimids[1] = nc->ni_dimid;
        dcount[1] = nc->ni_size;
    }
    else {
        log_err("COORD_DIMS_OUT should be 1 or 2");
    }

    // define the netcdf variable longitude
    status = nc_def_var(nc->nc_id, "lon", NC_DOUBLE, ndims,
                        dimids, &(lon_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining lon variable in %s", nc->fname);
    }

    status = nc_put_att_text(nc->nc_id, lon_var_id, "long_name",
                             strlen("longitude"), "longitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, lon_var_id, "units",
                             strlen("degrees_east"), "degrees_east");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, lon_var_id, "standard_name",
                             strlen("longitude"), "longitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }

    if (options.COORD_DIMS_OUT == 1) {
        dimids[0] = nc->nj_dimid;
        dcount[0] = nc->nj_size;
    }

    // define the netcdf variable latitude
    status = nc_def_var(nc->nc_id, "lat", NC_DOUBLE, ndims,
                        dimids, &(lat_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining lat variable in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "long_name",
                             strlen("latitude"), "latitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "units",
                             strlen("degrees_north"), "degrees_north");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "standard_name",
                             strlen("latitude"), "latitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", nc->fname);
    }


    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        log_err("Error leaving define mode for %s", nc->fname);
    }

    // fill the netcdf variables lat/lon
    if (options.COORD_DIMS_OUT == 1) {
        dvar = (double *) malloc(nc->ni_size * sizeof(double));

        dcount[0] = nc->ni_size;
        for (i = 0; i < nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations_grid[i].longitude;
        }
        status = nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", nc->fname);
        }
        free(dvar);

        dvar = (double *) malloc(nc->nj_size * sizeof(double));
        dcount[0] = nc->nj_size;
        for (i = 0; i < nc->nj_size; i++) {
            dvar[i] =
                (double) global_domain.locations_grid[i +
                                                      (i *
                                                       (nc->ni_size -
                                                        1))].latitude;
        }

        status = nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", nc->fname);
        }
        free(dvar);
    }
    else if (options.COORD_DIMS_OUT == 2) {
        dvar = (double *) malloc(nc->nj_size * nc->ni_size * sizeof(double));

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations_grid[i].longitude;
        }
        status = nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", nc->fname);
        }

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations_grid[i].latitude;
        }
        status = nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lat in %s", nc->fname);
        }

        free(dvar);
    }
    else {
        log_err("COORD_DIMS_OUT should be 1 or 2");
    }
}
