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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Initialzie output structures and determine which variables to
 *           write
 *****************************************************************************/
void
vic_init_output(dmy_struct *dmy_current)
{
    extern option_struct   options;
    extern domain_struct   local_domain;
    extern filep_struct    filep;
    extern MPI_Comm        MPI_COMM_VIC;
    extern int             mpi_rank;
    extern double       ***out_data;
    extern stream_struct **output_streams;
    extern nc_file_struct *nc_hist_files;

    size_t                 i;

    // initialize the output data structures
    set_output_met_data_info();

    out_data = create_outdata(local_domain.ncells_active);

    if (mpi_rank == 0) {
        // determine which variables will be written to the history file
        parse_output_info(filep.globalparam, output_streams);
    }

    // for (i = 0; i < options.Noutstreams; i++) {
    // status = MPI_Bcast(&(output_streams[i]), 1, mpi_stream_struct_type, 0, MPI_COMM_VIC);
    // if (status != MPI_SUCCESS) {
    // log_err("MPI error in vic_init_output(): %d\n", status);
    // }
    // }

    // allocation of output structures
    output_streams =
        calloc(local_domain.ncells_active, sizeof(*output_streams));
    if (output_streams == NULL) {
        log_err("Memory allocation error in vic_init_output().");
    }

    if (mpi_rank == 0) {
        for (i = 0; i < options.Noutstreams; i++) {
            nc_hist_files = calloc(options.Noutstreams, sizeof(*nc_hist_files));
            if (nc_hist_files == NULL) {
                log_err("Memory allocation error in vic_init_output().");
            }
        }

        for (i = 0; i < options.Noutstreams; i++) {
            // open the netcdf history file
            initialize_history_file(&(nc_hist_files[i]), &(*output_streams[i]),
                                    dmy_current);
        }
    }
}

/******************************************************************************
 * @brief    Initialize history files
 *****************************************************************************/
void
initialize_history_file(nc_file_struct *nc,
                        stream_struct  *stream,
                        dmy_struct     *dmy_current)
{
    extern filenames_struct    filenames;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern global_param_struct global_param;
    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    int                        status;
    int                        old_fill_mode;

    char                       str[MAXSTRING];
    char                       unit_str[MAXSTRING];
    char                       calendar_str[MAXSTRING];
    char                       cell_method[MAXSTRING];

    size_t                     i;
    size_t                     j;
    size_t                     ndims;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    int                        dimids[MAXDIMS];
    int                        time_var_id;
    int                        lon_var_id;
    int                        lat_var_id;
    unsigned int               varid;
    double                    *dvar;


    // This could be further refined but for now, I've choosen a file naming
    // Convention that goes like this:
    // If FREQ_NDAYS -- filename = result_dir/prefix.YYYY-MM-DD.nc
    if (stream->agg_alarm.freq == FREQ_NDAYS) {
        sprintf(stream->filename, "%s/%s.%04d-%02d-%02d.nc",
                filenames.result_dir,
                stream->prefix, dmy_current->year, dmy_current->month,
                dmy_current->day);
    }
    // If FREQ_NMONTHS -- filename = result_dir/prefix.YYYY-MM.nc
    else if (stream->agg_alarm.freq == FREQ_NMONTHS) {
        sprintf(stream->filename, "%s/%s.%04d-%02d.nc", filenames.result_dir,
                stream->prefix, dmy_current->year, dmy_current->month);
    }
    // If FREQ_NYEARS -- filename = result_dir/prefix.YYYY.nc
    else if (stream->agg_alarm.freq == FREQ_NYEARS) {
        sprintf(stream->filename, "%s/%s.%04d.nc", filenames.result_dir,
                stream->prefix, dmy_current->year);
    }
    // For all other cases -- filename = result_dir/prefix.YYYY-MM-DD-SSSSS.nc
    else {
        sprintf(stream->filename, "%s/%s.%04d-%02d-%02d-%05u.nc",
                filenames.result_dir,
                stream->prefix, dmy_current->year, dmy_current->month,
                dmy_current->day, dmy_current->dayseconds);
    }

    initialize_nc_file(nc, stream->nvars, stream->varid);

    // open the netcdf file
    status =
        nc_create(stream->filename, get_nc_mode(stream->file_format),
                  &(nc->nc_id));
    if (status != NC_NOERR) {
        log_err("Error creating %s", stream->filename);
    }
    nc->open = true;

    // Set netcdf file global attributes
    set_global_nc_attributes(nc->nc_id, NC_HISTORY_FILE);

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    if (status != NC_NOERR) {
        log_err("Error setting fill value in %s", stream->filename);
    }

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining snow_band dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "front", nc->front_size,
                        &(nc->front_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining front dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining frost_area dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining nlayer dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, global_domain.info.x_dim, nc->ni_size,
                        &(nc->ni_dimid));

    if (status != NC_NOERR) {
        log_err("Error defining x dimenension in %s", stream->filename);
    }
    status = nc_def_dim(nc->nc_id, global_domain.info.y_dim, nc->nj_size,
                        &(nc->nj_dimid));

    if (status != NC_NOERR) {
        log_err("Error defining y dimenension in %s", stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining node dimenension in %s", stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining root_zone dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining veg_class dimenension in %s",
                stream->filename);
    }

    status = nc_def_dim(nc->nc_id, "time", nc->time_size,
                        &(nc->time_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining time dimenension in %s", stream->filename);
    }

    // define the netcdf variable time
    status = nc_def_var(nc->nc_id, "time", NC_DOUBLE, 1,
                        &(nc->time_dimid), &(time_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining time variable in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, time_var_id, "standard_name",
                             strlen("time"), "time");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }

    // adding units attribute to time variable
    str_from_time_units(global_param.time_units, unit_str);

    sprintf(str, "%s since %s", unit_str, global_param.time_origin_str);

    status = nc_put_att_text(nc->nc_id, time_var_id, "units",
                             strlen(str), str);
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }

    // adding calendar attribute to time variable
    str_from_calendar(global_param.calendar, calendar_str);

    status = nc_put_att_text(nc->nc_id, time_var_id, "calendar",
                             strlen(calendar_str), calendar_str);
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }

    ndims = global_domain.info.n_coord_dims;
    dstart[0] = 0;
    dstart[1] = 0;

    if (global_domain.info.n_coord_dims == 1) {
        dimids[0] = nc->ni_dimid;
        dcount[0] = nc->ni_size;
    }
    else if (global_domain.info.n_coord_dims == 2) {
        dimids[0] = nc->nj_dimid;
        dcount[0] = nc->nj_size;

        dimids[1] = nc->ni_dimid;
        dcount[1] = nc->ni_size;
    }
    else {
        log_err("n_coord_dims should be 1 or 2");
    }

    // define the netcdf variable longitude
    status =
        nc_def_var(nc->nc_id, global_domain.info.lon_var, NC_DOUBLE, ndims,
                   dimids, &(lon_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining lon variable in %s", stream->filename);
    }

    status = nc_put_att_text(nc->nc_id, lon_var_id, "long_name",
                             strlen("longitude"), "longitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, lon_var_id, "units",
                             strlen("degrees_east"), "degrees_east");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, lon_var_id, "standard_name",
                             strlen("longitude"), "longitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }

    if (global_domain.info.n_coord_dims == 1) {
        dimids[0] = nc->nj_dimid;
        dcount[0] = nc->nj_size;
    }

    // define the netcdf variable latitude
    status = nc_def_var(nc->nc_id, global_domain.info.lat_var, NC_DOUBLE, ndims,
                        dimids, &(lat_var_id));
    if (status != NC_NOERR) {
        log_err("Error defining lat variable in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "long_name",
                             strlen("latitude"), "latitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "units",
                             strlen("degrees_north"), "degrees_north");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }
    status = nc_put_att_text(nc->nc_id, lat_var_id, "standard_name",
                             strlen("latitude"), "latitude");
    if (status != NC_NOERR) {
        log_err("Error adding attribute in %s", stream->filename);
    }

    // create output variables
    for (j = 0; j < stream->nvars; j++) {
        varid = stream->varid[j];

        // define the variable
        status = nc_def_var(nc->nc_id,
                            nc->nc_vars[j].nc_var_name,
                            nc->nc_vars[j].nc_type,
                            nc->nc_vars[j].nc_dims,
                            nc->nc_vars[j].nc_dimids,
                            &(nc->nc_vars[j].nc_varid));
        if (status != NC_NOERR) {
            log_err("Error defining variable %s in %s",
                    nc->nc_vars[j].nc_var_name, stream->filename);
        }

        // Add compression (only works for netCDF4 filetype)
        if (stream->compress) {
            status = nc_def_var_deflate(nc->nc_id, nc->nc_vars[j].nc_varid,
                                        true, true, stream->compress);
            if (status != NC_NOERR) {
                log_err(
                    "Error setting compression level in %s for variable: %s",
                    stream->filename, nc->nc_vars[j].nc_var_name);
            }
        }

        // set the fill value attribute
        if (nc->nc_vars[j].nc_type == NC_DOUBLE) {
            status = nc_put_att_double(nc->nc_id, nc->nc_vars[j].nc_varid,
                                       "_FillValue", NC_DOUBLE, 1,
                                       &(nc->d_fillvalue));
        }
        else if (nc->nc_vars[j].nc_type == NC_FLOAT) {
            status = nc_put_att_float(nc->nc_id, nc->nc_vars[j].nc_varid,
                                      "_FillValue", NC_FLOAT, 1,
                                      &(nc->f_fillvalue));
        }
        else if (nc->nc_vars[j].nc_type == NC_INT) {
            status = nc_put_att_int(nc->nc_id, nc->nc_vars[j].nc_varid,
                                    "_FillValue", NC_INT, 1,
                                    &(nc->i_fillvalue));
        }
        else {
            // TODO: Add NC_SHORT, NC_CHAR, and NC_BYTE
            log_err("NC_TYPE %d not supported at this time",
                    nc->nc_vars[j].nc_type);
        }
        if (status != NC_NOERR) {
            log_err("Error putting _FillValue attribute to %s in %s",
                    nc->nc_vars[j].nc_var_name, stream->filename);
        }

        put_nc_attr(nc->nc_id, nc->nc_vars[j].nc_varid, "long_name",
                    out_metadata[varid].long_name);
        put_nc_attr(nc->nc_id, nc->nc_vars[j].nc_varid, "standard_name",
                    out_metadata[varid].standard_name);
        put_nc_attr(nc->nc_id, nc->nc_vars[j].nc_varid, "units",
                    out_metadata[varid].units);
        put_nc_attr(nc->nc_id, nc->nc_vars[j].nc_varid, "description",
                    out_metadata[varid].description);

        if (cell_method_from_agg_type(stream->aggtype[j], cell_method)) {
            put_nc_attr(nc->nc_id, nc->nc_vars[j].nc_varid, "cell_methods",
                        cell_method);
            // NOTE: if cell_methods == variance, units should be ^2
        }
    }

    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        log_err("Error leaving define mode for %s", stream->filename);
    }

    // fill the netcdf variables lat/lon
    if (global_domain.info.n_coord_dims == 1) {
        dvar = calloc(nc->ni_size, sizeof(*dvar));
        if (dvar == NULL) {
            log_err("Memory allocation error in vic_init_output().");
        }

        dcount[0] = nc->ni_size;
        for (i = 0; i < nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].longitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", stream->filename);
        }
        free(dvar);

        dvar = calloc(nc->nj_size, sizeof(*dvar));
        if (dvar == NULL) {
            log_err("Memory allocation error in vic_init_output().");
        }
        dcount[0] = nc->nj_size;
        for (i = 0; i < nc->nj_size; i++) {
            dvar[i] =
                (double) global_domain.locations[i * nc->ni_size].latitude;
        }

        status =
            nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", stream->filename);
        }
        free(dvar);
    }
    else if (global_domain.info.n_coord_dims == 2) {
        dvar = calloc(nc->nj_size * nc->ni_size, sizeof(*dvar));
        if (dvar == NULL) {
            log_err("Memory allocation error in vic_init_output().");
        }

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].longitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lon in %s", stream->filename);
        }

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].latitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        if (status != NC_NOERR) {
            log_err("Error adding data to lat in %s", stream->filename);
        }

        free(dvar);
    }
    else {
        log_err("n_coord_dims should be 1 or 2");
    }
}

/******************************************************************************
 * @brief    Determine the netCDF file format
 *****************************************************************************/
int
get_nc_mode(unsigned short int format)
{
    if (format == NETCDF3_CLASSIC) {
        return NC_CLASSIC_MODEL;
    }
    else if (format == NETCDF3_64BIT_OFFSET) {
        return NC_64BIT_OFFSET;
    }
    else if (format == NETCDF4_CLASSIC) {
        return (NC_NETCDF4 | NC_CLASSIC_MODEL);
    }
    else if (format == NETCDF4) {
        return NC_NETCDF4;
    }
    else {
        log_err("Unrecognized netCDF file format");
    }
}

/******************************************************************************
 * @brief    Set global netcdf attributes (either history or state file)
 *****************************************************************************/
void
set_global_nc_attributes(int ncid,
                         unsigned short int
                         file_type)
{
    char       tmpstr[MAXSTRING];
    char       userstr[MAXSTRING];
    char       hoststr[MAXSTRING];
    char       mpistr[MPI_MAX_LIBRARY_VERSION_STRING];
    int        len;
    int        status;
    time_t     curr_date_time;
    struct tm *timeinfo;

    // datestr
    curr_date_time = time(NULL);
    if (curr_date_time == -1) {
        log_err("Something went wrong getting the current time!");
    }
    timeinfo = localtime(&curr_date_time);

    // username
    if (getlogin_r(userstr, MAXSTRING) != 0) {
        log_err("Error getting username");
    }
    // hostname
    if (gethostname(hoststr, MAXSTRING) != 0) {
        log_err("Error getting hostname");
    }

    // Set global attributes
    if (file_type == NC_HISTORY_FILE) {
        put_nc_attr(ncid, NC_GLOBAL, "title", "VIC History File");
    }
    else if (file_type == NC_STATE_FILE) {
        put_nc_attr(ncid, NC_GLOBAL, "title", "VIC State File");
    }
    else {
        put_nc_attr(ncid, NC_GLOBAL, "title", "Unknown");
    }

    // TODO: pass in driver as an argmument to this function
    put_nc_attr(ncid, NC_GLOBAL, "source", "VIC Image Driver");
    sprintf(tmpstr, "Created by %s on %s on %s",
            userstr, hoststr, asctime(timeinfo));
    put_nc_attr(ncid, NC_GLOBAL, "history", tmpstr);
    put_nc_attr(ncid, NC_GLOBAL, "references",
                "Primary Historical Reference for VIC: Liang, X., D. P. "
                "Lettenmaier, E. F. Wood, and S. J. Burges, 1994: A Simple "
                "hydrologically Based Model of Land Surface Water and Energy "
                "Fluxes for GSMs, J. Geophys. Res., 99(D7), 14,415-14,428.");
    put_nc_attr(ncid, NC_GLOBAL, "comment",
                "Output from the Variable Infiltration Capacity (VIC)"
                "Macroscale Hydrologic Model");
    put_nc_attr(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
    put_nc_attr(ncid, NC_GLOBAL, "netcdf_lib_version", nc_inq_libvers());
    status = MPI_Get_library_version(mpistr, &len);
    if (status == MPI_SUCCESS) {
        put_nc_attr(ncid, NC_GLOBAL, "mpi_lib_version", mpistr);
    }

    // Useful attributes from VIC
    put_nc_attr(ncid, NC_GLOBAL, "VIC_Model_Version", VERSION);
    // TODO: pass in driver as an argmument to this function
    put_nc_attr(ncid, NC_GLOBAL, "VIC_Driver", "Image");
}

/******************************************************************************
 * @brief    Set global netcdf attributes (either history or state file)
 *****************************************************************************/
void
initialize_nc_file(nc_file_struct *nc_file,
                   size_t          nvars,
                   unsigned int   *varids)
{
    extern option_struct options;
    extern domain_struct global_domain;

    size_t               i;

    // Set fill values
    nc_file->c_fillvalue = NC_FILL_CHAR;
    nc_file->i_fillvalue = NC_FILL_INT;
    nc_file->d_fillvalue = NC_FILL_DOUBLE;
    nc_file->f_fillvalue = NC_FILL_FLOAT;

    // set ids to MISSING
    nc_file->nc_id = MISSING;
    nc_file->band_dimid = MISSING;
    nc_file->front_dimid     MISSING;
    nc_file->frost_dimid     MISSING;
    nc_file->lake_node_dimid MISSING;
    nc_file->layer_dimid     MISSING;
    nc_file->ni_dimid = MISSING;
    nc_file->nj_dimid = MISSING;
    nc_file->node_dimid = MISSING;
    nc_file->root_zone_dimid MISSING;
    nc_file->time_dimid = MISSING;
    nc_file->veg_dimid = MISSING;

    // Set dimension sizes
    nc_file->band_size = options.SNOW_BAND;
    nc_file->front_size = MAX_FRONTS;
    nc_file->frost_size = options.Nfrost;
    nc_file->layer_size = options.Nlayer;
    nc_file->ni_size = global_domain.n_nx;
    nc_file->nj_size = global_domain.n_ny;
    nc_file->node_size = options.Nnode;
    nc_file->root_zone_size = options.ROOT_ZONES;
    nc_file->time_size = NC_UNLIMITED;
    nc_file->veg_size = options.NVEGTYPES;

    // allocate memory for nc_vars
    nc_file->nc_vars = calloc(nvars, sizeof(*(nc_file->nc_vars)));
    if (nc_file->nc_vars == NULL) {
        log_err("Memory allocation error");
    }

    for (i = 0; i < nvars; i++) {
        set_nc_var_info(varids[i], nc_file, &(nc_file->nc_vars[i]));
    }
}
