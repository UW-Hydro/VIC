/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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
 * @brief    Save model state.
 *****************************************************************************/
void
vic_store(dmy_struct *dmy_current)
{
    extern size_t              current;
    extern dmy_struct         *dmy;
    extern filenames_struct    filenames;
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern int                 mpi_rank;

    int                        status;
    int                        old_fill_mode;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    int                        dimids[MAXDIMS];
    int                        v;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     p;
    size_t                     ndims;
    char                      *cvar = NULL;
    int                       *ivar = NULL;
    double                    *dvar = NULL;
    float                     *fvar = NULL;
    int                        lon_var_id;
    int                        lat_var_id;
    size_t                     d1count[1];
    size_t                     d1start[1];
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

    nc_file_struct             nc_state_file;

    // create netcdf file for storing model state
    sprintf(nc_state_file.fname, "%s.%04d%02d%02d_%05u.nc",
            filenames.statefile, dmy[current].year, dmy[current].month,
            dmy[current].day, dmy[current].dayseconds);

    nc_state_file.c_fillvalue = NC_FILL_CHAR;
    nc_state_file.i_fillvalue = NC_FILL_INT;
    nc_state_file.d_fillvalue = NC_FILL_DOUBLE;
    nc_state_file.f_fillvalue = NC_FILL_FLOAT;

    nc_state_file.ni_size = global_domain.n_nx;
    nc_state_file.nj_size = global_domain.n_ny;
    nc_state_file.veg_size = options.NVEGTYPES;
    nc_state_file.band_size = options.SNOW_BAND;
    nc_state_file.layer_size = options.Nlayer;
    nc_state_file.frost_size = options.Nfrost;
    nc_state_file.node_size = options.Nnode;
    if (options.LAKES) {
        nc_state_file.lake_node_size = options.NLAKENODES;
    }

    // only open and initialize the netcdf file on the first thread
    if (mpi_rank == 0) {
        // open the netcdf file
        status = nc_create(nc_state_file.fname, get_nc_mode(
                               options.STATE_FORMAT), &(nc_state_file.nc_id));
        if (status != NC_NOERR) {
            log_err("Error creating %s", nc_state_file.fname);
        }
        nc_state_file.open = true;

        // Set netcdf file global attributes
        set_global_nc_attributes(nc_state_file.nc_id, NC_STATE_FILE);

        // set the NC_FILL attribute
        status = nc_set_fill(nc_state_file.nc_id, NC_FILL, &old_fill_mode);
        if (status != NC_NOERR) {
            log_err("Error setting fill value in %s", nc_state_file.fname);
        }

        // define netcdf dimensions
        status = nc_def_dim(nc_state_file.nc_id, global_domain.info.x_dim,
                            nc_state_file.ni_size, &(nc_state_file.ni_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining \"%s\" in %s", global_domain.info.x_dim,
                    nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, global_domain.info.y_dim,
                            nc_state_file.nj_size, &(nc_state_file.nj_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining \"%s\" in %s", global_domain.info.y_dim,
                    nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, "veg_class",
                            nc_state_file.veg_size, &(nc_state_file.veg_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining veg_class in %s", nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, "snow_band",
                            nc_state_file.band_size,
                            &(nc_state_file.band_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining snow_band in %s", nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, "nlayer",
                            nc_state_file.layer_size,
                            &(nc_state_file.layer_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining nlayer in %s", nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, "frost_area",
                            nc_state_file.frost_size,
                            &(nc_state_file.frost_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining frost_area in %s", nc_state_file.fname);
        }

        status = nc_def_dim(nc_state_file.nc_id, "soil_node",
                            nc_state_file.node_size,
                            &(nc_state_file.node_dimid));
        if (status != NC_NOERR) {
            log_err("Error defining soil_node in %s", nc_state_file.fname);
        }

        if (options.LAKES) {
            status = nc_def_dim(nc_state_file.nc_id, "lake_node",
                                nc_state_file.lake_node_size,
                                &(nc_state_file.lake_node_dimid));
            if (status != NC_NOERR) {
                log_err("Error defining lake_node in %s", nc_state_file.fname);
            }
        }

        // initialize dimids to invalid values
        for (i = 0; i < MAXDIMS; i++) {
            dimids[i] = -1;
        }

        // write dimension variables

        // Coordinate variables
        ndims = global_domain.info.n_coord_dims;
        dstart[0] = 0;
        dstart[1] = 0;

        if (global_domain.info.n_coord_dims == 1) {
            dimids[0] = nc_state_file.ni_dimid;
            dcount[0] = nc_state_file.ni_size;
        }
        else if (global_domain.info.n_coord_dims == 2) {
            dimids[0] = nc_state_file.nj_dimid;
            dcount[0] = nc_state_file.nj_size;

            dimids[1] = nc_state_file.ni_dimid;
            dcount[1] = nc_state_file.ni_size;
        }
        else {
            log_err("COORD_DIMS_OUT should be 1 or 2");
        }

        // define the netcdf variable longitude
        status = nc_def_var(nc_state_file.nc_id, global_domain.info.lon_var,
                            NC_DOUBLE, ndims, dimids, &(lon_var_id));
        if (status != NC_NOERR) {
            log_err("Error defining lon variable in %s", nc_state_file.fname);
        }

        status = nc_put_att_text(nc_state_file.nc_id, lon_var_id, "long_name", strlen(
                                     "longitude"), "longitude");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }
        status = nc_put_att_text(nc_state_file.nc_id, lon_var_id, "units", strlen(
                                     "degrees_east"), "degrees_east");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }
        status = nc_put_att_text(nc_state_file.nc_id, lon_var_id,
                                 "standard_name", strlen(
                                     "longitude"), "longitude");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }

        if (global_domain.info.n_coord_dims == 1) {
            dimids[0] = nc_state_file.nj_dimid;
            dcount[0] = nc_state_file.nj_size;
        }

        // define the netcdf variable latitude
        status = nc_def_var(nc_state_file.nc_id, global_domain.info.lat_var,
                            NC_DOUBLE, ndims, dimids, &(lat_var_id));
        if (status != NC_NOERR) {
            log_err("Error defining lat variable in %s", nc_state_file.fname);
        }
        status = nc_put_att_text(nc_state_file.nc_id, lat_var_id, "long_name", strlen(
                                     "latitude"), "latitude");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }
        status = nc_put_att_text(nc_state_file.nc_id, lat_var_id, "units", strlen(
                                     "degrees_north"), "degrees_north");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }
        status = nc_put_att_text(nc_state_file.nc_id, lat_var_id,
                                 "standard_name", strlen("latitude"),
                                 "latitude");
        if (status != NC_NOERR) {
            log_err("Error adding attribute in %s", nc_state_file.fname);
        }

        // leave define mode
        status = nc_enddef(nc_state_file.nc_id);
        if (status != NC_NOERR) {
            log_err("Error leaving define mode for %s", nc_state_file.fname);
        }

        // populate lat/lon
        if (global_domain.info.n_coord_dims == 1) {
            dvar = calloc(nc_state_file.ni_size, sizeof(*dvar));
            if (dvar == NULL) {
                log_err("Memory allocation error in vic_store().");
            }

            dcount[0] = nc_state_file.ni_size;
            // implicitly nested loop over ni and nj with j set to 0
            for (i = 0; i < nc_state_file.ni_size; i++) {
                dvar[i] = (double) global_domain.locations[i].longitude;
            }
            status = nc_put_vara_double(nc_state_file.nc_id, lon_var_id, dstart,
                                        dcount, dvar);
            if (status != NC_NOERR) {
                log_err("Error adding data to lon in %s", nc_state_file.fname);
            }
            free(dvar);

            dvar = calloc(nc_state_file.nj_size, sizeof(*dvar));
            if (dvar == NULL) {
                log_err("Memory allocation error in vic_store().");
            }
            dcount[0] = nc_state_file.nj_size;
            // implicitly nested loop over ni and nj with i set to 0;
            // j stride = ni_size
            for (j = 0; j < nc_state_file.nj_size; j++) {
                dvar[j] =
                    (double) global_domain.locations[j *
                                                     nc_state_file.ni_size].
                    latitude;
            }

            status = nc_put_vara_double(nc_state_file.nc_id, lat_var_id, dstart,
                                        dcount, dvar);
            if (status != NC_NOERR) {
                log_err("Error adding data to lon in %s", nc_state_file.fname);
            }
            free(dvar);
        }
        else if (global_domain.info.n_coord_dims == 2) {
            dvar =
                calloc(nc_state_file.nj_size * nc_state_file.ni_size,
                       sizeof(*dvar));
            if (dvar == NULL) {
                log_err("Memory allocation error in vic_store().");
            }

            for (i = 0; i < nc_state_file.nj_size * nc_state_file.ni_size;
                 i++) {
                dvar[i] = (double) global_domain.locations[i].longitude;
            }
            status = nc_put_vara_double(nc_state_file.nc_id, lon_var_id, dstart,
                                        dcount, dvar);
            if (status != NC_NOERR) {
                log_err("Error adding data to lon in %s", nc_state_file.fname);
            }

            for (i = 0; i < nc_state_file.nj_size * nc_state_file.ni_size;
                 i++) {
                dvar[i] = (double) global_domain.locations[i].latitude;
            }
            status = nc_put_vara_double(nc_state_file.nc_id, lat_var_id, dstart,
                                        dcount, dvar);
            if (status != NC_NOERR) {
                log_err("Error adding data to lat in %s", nc_state_file.fname);
            }

            free(dvar);
        }
        else {
            log_err("COORD_DIMS_OUT should be 1 or 2");
        }

        // Variables for other dimensions (all 1-dimensional)
        ndims = 1;
        d1start[0] = 0;

        // vegetation classes
        dimids[0] = nc_state_file.veg_dimid;
        d1count[0] = nc_state_file.veg_size;
        ivar = malloc(nc_state_file.veg_size * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in vic_store().");
        }
        for (j = 0; j < nc_state_file.veg_size; j++) {
            ivar[j] = (int) j + 1;
        }
        put_nc_field_int(nc_state_file.fname, &(nc_state_file.open),
                         &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                         dimids, ndims,
                         "veg_class", d1start, d1count, ivar);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
        free(ivar);

        // snow bands
        dimids[0] = nc_state_file.band_dimid;
        d1count[0] = nc_state_file.band_size;
        ivar = malloc(nc_state_file.band_size * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in vic_store().");
        }
        for (j = 0; j < nc_state_file.band_size; j++) {
            ivar[j] = (int) j;
        }
        put_nc_field_int(nc_state_file.fname, &(nc_state_file.open),
                         &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                         dimids, ndims,
                         "snow_band", d1start, d1count, ivar);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
        free(ivar);

        // soil layers
        dimids[0] = nc_state_file.layer_dimid;
        d1count[0] = nc_state_file.layer_size;
        ivar = malloc(nc_state_file.layer_size * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in vic_store().");
        }
        for (j = 0; j < nc_state_file.layer_size; j++) {
            ivar[j] = (int) j;
        }
        put_nc_field_int(nc_state_file.fname, &(nc_state_file.open),
                         &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                         dimids, ndims,
                         "layer", d1start, d1count, ivar);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
        free(ivar);

        // frost areas
        dimids[0] = nc_state_file.frost_dimid;
        d1count[0] = nc_state_file.frost_size;
        ivar = malloc(nc_state_file.frost_size * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in vic_store().");
        }
        for (j = 0; j < nc_state_file.frost_size; j++) {
            ivar[j] = (int) j;
        }
        put_nc_field_int(nc_state_file.fname, &(nc_state_file.open),
                         &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                         dimids, ndims,
                         "frost_area", d1start, d1count, ivar);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
        free(ivar);

        // soil thermal node deltas
        dimids[0] = nc_state_file.node_dimid;
        d1count[0] = nc_state_file.node_size;
        put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                            &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                            dimids, ndims,
                            "dz_node", d1start, d1count, soil_con[0].dz_node);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // soil thermal node depths
        dimids[0] = nc_state_file.node_dimid;
        d1count[0] = nc_state_file.node_size;
        put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                            &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                            dimids, ndims,
                            "node_depth", d1start, d1count,
                            soil_con[0].Zsum_node);
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        if (options.LAKES) {
            // lake nodes
            dimids[0] = nc_state_file.lake_node_dimid;
            d1count[0] = nc_state_file.lake_node_size;
            ivar = malloc(nc_state_file.lake_node_size * sizeof(*ivar));
            if (ivar == NULL) {
                log_err("Memory allocation error in vic_store().");
            }
            for (j = 0; j < nc_state_file.lake_node_size; j++) {
                ivar[j] = (int) j;
            }
            put_nc_field_int(nc_state_file.fname, &(nc_state_file.open),
                             &(nc_state_file.nc_id), nc_state_file.d_fillvalue,
                             dimids, ndims,
                             "lake_node", d1start, d1count, ivar);
            for (i = 0; i < ndims; i++) {
                dimids[i] = -1;
            }
            free(ivar);
        }
    } // end if (mpi_rank == 0)


    // write state variables

    // allocate memory for variables to be stored
    cvar = malloc(local_domain.ncells_active * sizeof(*cvar));
    if (cvar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    if (ivar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    fvar = malloc(local_domain.ncells_active * sizeof(*fvar));
    if (fvar == NULL) {
        log_err("Memory allocation error in vic_store().");
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
    for (i = 0; i < local_domain.ncells_active; i++) {
        cvar[i] = nc_state_file.c_fillvalue;
        ivar[i] = nc_state_file.i_fillvalue;
        dvar[i] = nc_state_file.d_fillvalue;
        fvar[i] = nc_state_file.f_fillvalue;
    }

    // total soil moisture
    ndims = 5;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.layer_dimid;
    dimids[3] = nc_state_file.nj_dimid;
    dimids[4] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d5start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d5start[2] = j;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][k].layer[j].moist;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Soil_moisture",
                                           d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // ice content
    ndims = 6;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.layer_dimid;
    dimids[3] = nc_state_file.frost_dimid;
    dimids[4] = nc_state_file.nj_dimid;
    dimids[5] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d6start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d6start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d6start[2] = j;
                for (p = 0; p < options.Nfrost; p++) {
                    d6start[3] = p;
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        v = veg_con_map[i].vidx[m];
                        if (v >= 0) {
                            dvar[i] = (double)
                                      all_vars[i].cell[v][k].layer[j].ice[p];
                        }
                        else {
                            dvar[i] = nc_state_file.d_fillvalue;
                        }
                    }
                    gather_put_nc_field_double(nc_state_file.fname,
                                               &(nc_state_file.open),
                                               &(nc_state_file.nc_id),
                                               nc_state_file.d_fillvalue,
                                               dimids, ndims, "Soil_ice",
                                               d6start, d6count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // dew storage: tmpval = veg_var[veg][band].Wdew;
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].veg_var[v][k].Wdew;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Canopy_water",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    if (options.CARBON) {
        // cumulative NPP: tmpval = veg_var[veg][band].AnnualNPP;
        ndims = 4;
        dimids[0] = nc_state_file.veg_dimid;
        dimids[1] = nc_state_file.band_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].veg_var[v][k].AnnualNPP;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "AnnualNPP",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // previous NPP: tmpval = veg_var[veg][band].AnnualNPPPrev;
        ndims = 4;
        dimids[0] = nc_state_file.veg_dimid;
        dimids[1] = nc_state_file.band_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].veg_var[v][k].AnnualNPPPrev;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "AnnualNPPPrev",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // litter carbon: tmpval = cell[veg][band].CLitter;
        ndims = 4;
        dimids[0] = nc_state_file.veg_dimid;
        dimids[1] = nc_state_file.band_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][k].CLitter;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "CLitter",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // intermediate carbon: tmpval = tmpval = cell[veg][band].CInter;
        ndims = 4;
        dimids[0] = nc_state_file.veg_dimid;
        dimids[1] = nc_state_file.band_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][k].CInter;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "CInter",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // slow carbon: tmpval = cell[veg][band].CSlow;
        ndims = 4;
        dimids[0] = nc_state_file.veg_dimid;
        dimids[1] = nc_state_file.band_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][k].CSlow;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "CSlow",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
    }

    // snow age: snow[veg][band].last_snow
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    ivar[i] = (int)
                              all_vars[i].snow[v][k].last_snow;
                }
                else {
                    ivar[i] = nc_state_file.i_fillvalue;
                }
            }
            gather_put_nc_field_int(nc_state_file.fname,
                                    &(nc_state_file.open),
                                    &(nc_state_file.nc_id),
                                    nc_state_file.i_fillvalue,
                                    dimids, ndims, "Snow_age",
                                    d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                ivar[i] = nc_state_file.i_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // melting state: (int)snow[veg][band].MELTING
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    ivar[i] = (int)
                              all_vars[i].snow[v][k].MELTING;
                }
                else {
                    ivar[i] = nc_state_file.i_fillvalue;
                }
            }
            gather_put_nc_field_int(nc_state_file.fname,
                                    &(nc_state_file.open),
                                    &(nc_state_file.nc_id),
                                    nc_state_file.i_fillvalue,
                                    dimids, ndims, "Snow_melt_state",
                                    d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                ivar[i] = nc_state_file.i_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow covered fraction: snow[veg][band].coverage
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].coverage;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_coverage",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow water equivalent: snow[veg][band].swq
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].swq;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_water_equivalent",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow surface temperature: snow[veg][band].surf_temp
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].surf_temp;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_surf_temp",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow surface water: snow[veg][band].surf_water
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].surf_water;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_surf_water",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow pack temperature: snow[veg][band].pack_temp
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].pack_temp;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_pack_temp",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow pack water: snow[veg][band].pack_water
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].pack_water;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_pack_water",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow density: snow[veg][band].density
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].density;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_density",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow cold content: snow[veg][band].coldcontent
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].coldcontent;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_cold_content",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow canopy storage: snow[veg][band].snow_canopy
    ndims = 4;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][k].snow_canopy;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_canopy",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // soil node temperatures: energy[veg][band].T[nidx]
    ndims = 5;
    dimids[0] = nc_state_file.veg_dimid;
    dimids[1] = nc_state_file.band_dimid;
    dimids[2] = nc_state_file.node_dimid;
    dimids[3] = nc_state_file.nj_dimid;
    dimids[4] = nc_state_file.ni_dimid;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d5start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nnode; j++) {
                d5start[2] = j;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].energy[v][k].T[j];
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Soil_node_temp",
                                           d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    if (options.LAKES) {
        // total soil moisture
        ndims = 3;
        dimids[0] = nc_state_file.layer_dimid;
        dimids[1] = nc_state_file.nj_dimid;
        dimids[2] = nc_state_file.ni_dimid;
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.soil.layer[j].moist;
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_soil_moisture",
                                       d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // ice content
        ndims = 4;
        dimids[0] = nc_state_file.layer_dimid;
        dimids[1] = nc_state_file.frost_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (j = 0; j < options.Nlayer; j++) {
            d4start[0] = j;
            for (p = 0; p < options.Nfrost; p++) {
                d4start[1] = p;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = (double)
                              all_vars[i].lake_var.soil.layer[j].ice[p];
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Lake_soil_ice",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        if (options.CARBON) {
            // litter carbon: tmpval = lake_var.soil.CLitter;
            ndims = 2;
            dimids[0] = nc_state_file.nj_dimid;
            dimids[1] = nc_state_file.ni_dimid;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.soil.CLitter;
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_CLitter",
                                       d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
            for (i = 0; i < ndims; i++) {
                dimids[i] = -1;
            }

            // intermediate carbon: tmpval = lake_var.soil.CInter;
            ndims = 2;
            dimids[0] = nc_state_file.nj_dimid;
            dimids[1] = nc_state_file.ni_dimid;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.soil.CInter;
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_CInter",
                                       d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
            for (i = 0; i < ndims; i++) {
                dimids[i] = -1;
            }

            // slow carbon: tmpval = lake_var.soil.CSlow;
            ndims = 2;
            dimids[0] = nc_state_file.nj_dimid;
            dimids[1] = nc_state_file.ni_dimid;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.soil.CSlow;
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_CSlow",
                                       d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
            for (i = 0; i < ndims; i++) {
                dimids[i] = -1;
            }
        }

        // snow age: lake_var.snow.last_snow
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = (int) all_vars[i].lake_var.snow.last_snow;
        }
        gather_put_nc_field_int(nc_state_file.fname,
                                &(nc_state_file.open),
                                &(nc_state_file.nc_id),
                                nc_state_file.i_fillvalue,
                                dimids, ndims, "Lake_snow_age",
                                d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // melting state: (int)lake_var.snow.MELTING
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = (int) all_vars[i].lake_var.snow.MELTING;
        }
        gather_put_nc_field_int(nc_state_file.fname,
                                &(nc_state_file.open),
                                &(nc_state_file.nc_id),
                                nc_state_file.i_fillvalue,
                                dimids, ndims, "Lake_snow_melt_state",
                                d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow covered fraction: lake_var.snow.coverage
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.coverage;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_coverage",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow water equivalent: lake_var.snow.swq
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.swq;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_water_equivalent",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow surface temperature: lake_var.snow.surf_temp
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.surf_temp;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_surf_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow surface water: lake_var.snow.surf_water
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.surf_water;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_surf_water",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow pack temperature: lake_var.snow.pack_temp
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.pack_temp;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_pack_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow pack water: lake_var.snow.pack_water
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.pack_water;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_pack_water",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow density: lake_var.snow.density
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.density;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_density",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow cold content: lake_var.snow.coldcontent
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.coldcontent;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_cold_content",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // snow canopy storage: lake_var.snow.snow_canopy
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.snow.snow_canopy;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_snow_canopy",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // soil node temperatures: lake_var.energy.T[nidx]
        ndims = 3;
        dimids[0] = nc_state_file.node_dimid;
        dimids[1] = nc_state_file.nj_dimid;
        dimids[2] = nc_state_file.ni_dimid;
        for (j = 0; j < options.Nnode; j++) {
            d3start[0] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.soil.layer[j].moist;
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_soil_node_temp",
                                       d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake active layers: lake_var.activenod
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = (int) all_vars[i].lake_var.activenod;
        }
        gather_put_nc_field_int(nc_state_file.fname,
                                &(nc_state_file.open),
                                &(nc_state_file.nc_id),
                                nc_state_file.i_fillvalue,
                                dimids, ndims, "Lake_active_layers",
                                d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            ivar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake layer thickness: lake_var.dz
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.dz;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_layer_dz",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake surface layer thickness: lake_var.surfdz
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.surfdz;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_surf_layer_dz",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake depth: lake_var.ldepth
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.ldepth;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_depth",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake layer surface areas: lake_var.surface[ndix]
        ndims = 3;
        dimids[0] = nc_state_file.lake_node_dimid;
        dimids[1] = nc_state_file.nj_dimid;
        dimids[2] = nc_state_file.ni_dimid;
        for (j = 0; j < options.NLAKENODES; j++) {
            d3start[0] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.surface[j];
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_layer_surf_area",
                                       d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake surface area: lake_var.sarea
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.sarea;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_surf_area",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake volume: lake_var.volume
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.volume;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_volume",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake layer temperatures: lake_var.temp[nidx]
        ndims = 3;
        dimids[0] = nc_state_file.lake_node_dimid;
        dimids[1] = nc_state_file.nj_dimid;
        dimids[2] = nc_state_file.ni_dimid;
        for (j = 0; j < options.NLAKENODES; j++) {
            d3start[0] = j;
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) all_vars[i].lake_var.temp[j];
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Lake_layer_temp",
                                       d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // vertical average lake temperature: lake_var.tempavg
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.tempavg;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_average_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice area fraction: lake_var.areai
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.areai;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_area_frac",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // new lake ice area fraction: lake_var.new_ice_area
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.new_ice_area;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_area_frac_new",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice water equivalent: lake_var.ice_water_eq
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.ice_water_eq;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_water_equivalent",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice height: lake_var.hice
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.hice;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_height",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice temperature: lake_var.tempi
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.tempi;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow water equivalent: lake_var.swe
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.swe;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims,
                                   "Lake_ice_snow_water_equivalen",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow surface temperature: lake_var.surf_temp
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.surf_temp;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_surf_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow pack temperature: lake_var.pack_temp
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.pack_temp;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_pack_temp",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow coldcontent: lake_var.coldcontent
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.coldcontent;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_cold_content",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow surface water: lake_var.surf_water
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.surf_water;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_surf_water",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow pack water: lake_var.pack_water
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.pack_water;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_pack_water",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow albedo: lake_var.SAlbedo
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.SAlbedo;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_albedo",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // lake ice snow depth: lake_var.sdepth
        ndims = 2;
        dimids[0] = nc_state_file.nj_dimid;
        dimids[1] = nc_state_file.ni_dimid;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) all_vars[i].lake_var.sdepth;
        }
        gather_put_nc_field_double(nc_state_file.fname,
                                   &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Lake_ice_snow_depth",
                                   d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.i_fillvalue;
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
    }

    // close the netcdf file if it is still open
    if (nc_state_file.open == true) {
        status = nc_close(nc_state_file.nc_id);
        if (status != NC_NOERR) {
            log_err("Error closing %s", nc_state_file.fname);
        }
    }

    free(cvar);
    free(ivar);
    free(dvar);
    free(fvar);
}
