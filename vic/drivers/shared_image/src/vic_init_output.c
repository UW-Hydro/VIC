/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize output structures.
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
 * @brief    Initialzie output structures and determine which variables to
 *           write
 *****************************************************************************/
void
vic_init_output(dmy_struct *dmy_current)
{
    extern all_vars_struct   *all_vars;
    extern force_data_struct *force;
    extern domain_struct      local_domain;
    extern filep_struct       filep;
    extern MPI_Comm           MPI_COMM_VIC;
    extern int                mpi_rank;
    extern nc_file_struct    *nc_hist_files;
    extern lake_con_struct    lake_con;
    extern double          ***out_data;
    extern save_data_struct  *save_data;
    extern soil_con_struct   *soil_con;
    extern veg_con_struct   **veg_con;
    extern veg_lib_struct   **veg_lib;
    extern option_struct      options;
    extern MPI_Datatype       mpi_alarm_struct_type;
    extern stream_struct     *output_streams;

    int                       status;
    size_t                    i;
    size_t                    streamnum;
    size_t                    nstream_vars[MAX_OUTPUT_STREAMS];
    bool                      default_outputs = false;
    timer_struct              timer;

    // initialize the output data structures
    set_output_met_data_info();

    // allocate out_data
    alloc_out_data(local_domain.ncells_active, out_data);

    // initialize the save data structures
    for (i = 0; i < local_domain.ncells_active; i++) {
        initialize_save_data(&(all_vars[i]), &(force[i]), &(soil_con[i]),
                             veg_con[i], veg_lib[i], &lake_con, out_data[i],
                             &(save_data[i]), &timer);
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        // count the number of streams and variables in the global parameter file
        count_nstreams_nvars(filep.globalparam, &(options.Noutstreams),
                             nstream_vars);

        // If there weren't any output streams specified, get the defaults
        if (options.Noutstreams == 0) {
            default_outputs = true;
            get_default_nstreams_nvars(&(options.Noutstreams), nstream_vars);
        }
    }

    // broadcast Noutstreams and nstream_vars
    status = MPI_Bcast(&(options.Noutstreams), 1, MPI_AINT, VIC_MPI_ROOT,
                       MPI_COMM_VIC);
    check_mpi_status(status, "MPI Error.");
    status = MPI_Bcast(&(nstream_vars), MAX_OUTPUT_STREAMS, MPI_AINT,
                       VIC_MPI_ROOT,
                       MPI_COMM_VIC);
    check_mpi_status(status, "MPI Error.");

    // allocate output streams
    output_streams = calloc(options.Noutstreams, sizeof(*output_streams));
    check_alloc_status(output_streams, "Memory allocation error.");

    // allocate netcdf history files array
    nc_hist_files = calloc(options.Noutstreams, sizeof(*nc_hist_files));
    check_alloc_status(nc_hist_files, "Memory allocation error.");

    // allocate memory for streams, initialize to default/missing values
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        setup_stream(&(output_streams[streamnum]), nstream_vars[streamnum],
                     local_domain.ncells_active);
    }

    if (mpi_rank == VIC_MPI_ROOT) {
        if (default_outputs) {
            // determine which variables will be written to the history file
            set_output_defaults(&output_streams, dmy_current, NETCDF4_CLASSIC);
        }
        else {
            // set output defaults
            parse_output_info(filep.globalparam, &output_streams, dmy_current);
        }
    }

    // Now broadcast the arrays of shape nvars
    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        // prefix
        status = MPI_Bcast(output_streams[streamnum].prefix,
                           MAXSTRING, MPI_CHAR, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // filename
        status = MPI_Bcast(output_streams[streamnum].filename,
                           MAXSTRING, MPI_CHAR, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // skip fh

        // file_format
        status = MPI_Bcast(&(output_streams[streamnum].file_format),
                           1, MPI_UNSIGNED_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // compress
        status = MPI_Bcast(&(output_streams[streamnum].compress),
                           1, MPI_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // type
        status = MPI_Bcast(output_streams[streamnum].type,
                           output_streams[streamnum].nvars,
                           MPI_UNSIGNED_SHORT, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // mult
        status = MPI_Bcast(output_streams[streamnum].mult,
                           output_streams[streamnum].nvars,
                           MPI_DOUBLE, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // format
        // skip broadcast

        // varid
        status = MPI_Bcast(output_streams[streamnum].varid,
                           output_streams[streamnum].nvars,
                           MPI_UNSIGNED, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // aggtype
        status = MPI_Bcast(output_streams[streamnum].aggtype,
                           output_streams[streamnum].nvars, MPI_UNSIGNED_SHORT,
                           VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // skip agg data

        // Now brodcast the alarms
        status = MPI_Bcast(&(output_streams[streamnum].agg_alarm), 1,
                           mpi_alarm_struct_type, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");
        status = MPI_Bcast(&(output_streams[streamnum].write_alarm), 1,
                           mpi_alarm_struct_type, VIC_MPI_ROOT, MPI_COMM_VIC);
        check_mpi_status(status, "MPI error.");

        // allocate agg data
        alloc_aggdata(&(output_streams[streamnum]));

        // setup netcdf files
        initialize_nc_file(&(nc_hist_files[streamnum]),
                           output_streams[streamnum].nvars,
                           output_streams[streamnum].varid,
                           output_streams[streamnum].type);
    }
    // validate streams
    validate_streams(&output_streams);
}

/******************************************************************************
 * @brief    Initialize history file
 *****************************************************************************/
void
initialize_history_file(nc_file_struct *nc,
                        stream_struct  *stream)
{
    extern filenames_struct    filenames;
    extern domain_struct       global_domain;
    extern option_struct       options;
    extern global_param_struct global_param;
    extern metadata_struct     out_metadata[N_OUTVAR_TYPES];

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
    int                        lon_var_id;
    int                        lat_var_id;
    unsigned int               varid;
    double                    *dvar;


    // This could be further refined but for now, I've chosen a file naming
    // Convention that goes like this:
    switch (stream->agg_alarm.freq) {
    // If FREQ_NDAYS -- filename = result_dir/prefix.YYYY-MM-DD.nc
    case FREQ_NDAYS:
        sprintf(stream->filename, "%s/%s.%04d-%02d-%02d.nc",
                filenames.result_dir,
                stream->prefix, stream->time_bounds[0].year,
                stream->time_bounds[0].month,
                stream->time_bounds[0].day);
        break;
    case FREQ_NMONTHS:
        // If FREQ_NMONTHS -- filename = result_dir/prefix.YYYY-MM.nc
        sprintf(stream->filename, "%s/%s.%04d-%02d.nc", filenames.result_dir,
                stream->prefix, stream->time_bounds[0].year,
                stream->time_bounds[0].month);
        break;
    case FREQ_NYEARS:
        // If FREQ_NYEARS -- filename = result_dir/prefix.YYYY.nc
        sprintf(stream->filename, "%s/%s.%04d.nc", filenames.result_dir,
                stream->prefix, stream->time_bounds[0].year);
        break;
    default:
        // For all other cases -- filename = result_dir/prefix.YYYY-MM-DD-SSSSS.nc
        sprintf(stream->filename, "%s/%s.%04d-%02d-%02d-%05u.nc",
                filenames.result_dir,
                stream->prefix, stream->time_bounds[0].year,
                stream->time_bounds[0].month,
                stream->time_bounds[0].day,
                stream->time_bounds[0].dayseconds);
    }

    // open the netcdf file
    status = nc_create(stream->filename,
                       get_nc_mode(stream->file_format),
                       &(nc->nc_id));
    check_nc_status(status, "Error creating %s", stream->filename);
    nc->open = true;

    // Set netcdf file global attributes
    set_global_nc_attributes(nc->nc_id, NC_HISTORY_FILE);

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    check_nc_status(status, "Error setting fill value in %s", stream->filename);

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    check_nc_status(status, "Error defining snow_band dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "front", nc->front_size,
                        &(nc->front_dimid));
    check_nc_status(status, "Error defining front dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    check_nc_status(status, "Error defining frost_area dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    check_nc_status(status, "Error defining nlayer dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, global_domain.info.x_dim, nc->ni_size,
                        &(nc->ni_dimid));

    check_nc_status(status, "Error defining x dimension in %s",
                    stream->filename);
    status = nc_def_dim(nc->nc_id, global_domain.info.y_dim, nc->nj_size,
                        &(nc->nj_dimid));

    check_nc_status(status, "Error defining y dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    check_nc_status(status, "Error defining node dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    check_nc_status(status, "Error defining root_zone dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    check_nc_status(status, "Error defining veg_class dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "time", nc->time_size,
                        &(nc->time_dimid));
    check_nc_status(status, "Error defining time dimension in %s",
                    stream->filename);

    status = nc_def_dim(nc->nc_id, "nv", 2, &(nc->time_bounds_dimid));
    check_nc_status(status, "Error defining time bounds dimension in %s",
                    stream->filename);

    if (options.LAKES) {
        status = nc_def_dim(nc->nc_id, "lake_node", nc->lake_node_size,
                            &(nc->lake_node_dimid));
        check_nc_status(status, "Error defining lake_node dimension in %s",
                        stream->filename);
    }

    // define the netcdf variable time
    status = nc_def_var(nc->nc_id, "time", NC_DOUBLE, 1,
                        &(nc->time_dimid), &(nc->time_varid));
    check_nc_status(status, "Error defining time variable in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, nc->time_varid, "standard_name",
                             strlen("time"), "time");
    check_nc_status(status, "Error adding attribute in %s", stream->filename);

    // adding units attribute to time variable
    str_from_time_units(global_param.time_units, unit_str);

    sprintf(str, "%s since %s", unit_str, global_param.time_origin_str);

    status = nc_put_att_text(nc->nc_id, nc->time_varid, "units",
                             strlen(str), str);
    check_nc_status(status, "Error adding attribute in %s", stream->filename);

    // adding calendar attribute to time variable
    str_from_calendar(global_param.calendar, calendar_str);

    status = nc_put_att_text(nc->nc_id, nc->time_varid, "calendar",
                             strlen(calendar_str), calendar_str);
    check_nc_status(status, "Error adding attribute in %s", stream->filename);

    // adding bounds attribute to time variable
    status = nc_put_att_text(nc->nc_id, nc->time_varid, "bounds",
                             strlen("time_bnds"), "time_bnds");
    check_nc_status(status, "Error adding attribute in %s", stream->filename);

    // define the netcdf variable time_bnds
    dimids[0] = nc->time_dimid;
    dimids[1] = nc->time_bounds_dimid;
    status = nc_def_var(nc->nc_id, "time_bnds", NC_DOUBLE, 2,
                        dimids, &(nc->time_bounds_varid));
    check_nc_status(status, "Error defining time bounds variable in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, nc->time_bounds_varid, "standard_name",
                             strlen("time_bounds"), "time_bounds");
    check_nc_status(status, "Error adding attribute in %s", stream->filename);
    status = nc_put_att_text(nc->nc_id, nc->time_bounds_varid, "units",
                             strlen(str), str);
    check_nc_status(status, "Error adding attribute in %s", stream->filename);
    status = nc_put_att_text(nc->nc_id, nc->time_bounds_varid, "calendar",
                             strlen(calendar_str), calendar_str);
    check_nc_status(status, "Error adding attribute in %s", stream->filename);

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
    check_nc_status(status, "Error defining lon variable in %s",
                    stream->filename);

    status = nc_put_att_text(nc->nc_id, lon_var_id, "long_name",
                             strlen("longitude"), "longitude");
    check_nc_status(status, "Error adding longitude long_name attribute in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, lon_var_id, "units",
                             strlen("degrees_east"), "degrees_east");
    check_nc_status(status, "Error adding longitude units attribute in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, lon_var_id, "standard_name",
                             strlen("longitude"), "longitude");
    check_nc_status(status,
                    "Error adding longitude standard_name attribute in %s",
                    stream->filename);

    if (global_domain.info.n_coord_dims == 1) {
        dimids[0] = nc->nj_dimid;
        dcount[0] = nc->nj_size;
    }

    // define the netcdf variable latitude
    status = nc_def_var(nc->nc_id, global_domain.info.lat_var, NC_DOUBLE, ndims,
                        dimids, &(lat_var_id));
    check_nc_status(status, "Error defining lat variable in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, lat_var_id, "long_name",
                             strlen("latitude"), "latitude");
    check_nc_status(status, "Error adding latitude long_name attribute in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, lat_var_id, "units",
                             strlen("degrees_north"), "degrees_north");
    check_nc_status(status, "Error adding latitude units attribute in %s",
                    stream->filename);
    status = nc_put_att_text(nc->nc_id, lat_var_id, "standard_name",
                             strlen("latitude"), "latitude");
    check_nc_status(status,
                    "Error adding latitude standard_name attribute in %s",
                    stream->filename);

    // create output variables
    for (j = 0; j < stream->nvars; j++) {
        varid = stream->varid[j];

        set_nc_var_dimids(varid, nc, &(nc->nc_vars[j]));

        // define the variable
        status = nc_def_var(nc->nc_id,
                            out_metadata[varid].varname,
                            nc->nc_vars[j].nc_type,
                            nc->nc_vars[j].nc_dims,
                            nc->nc_vars[j].nc_dimids,
                            &(nc->nc_vars[j].nc_varid));
        check_nc_status(status, "Error defining variable %s in %s.  Status: %d",
                        out_metadata[varid].varname, stream->filename, status);

        // Add compression (only works for netCDF4 filetype)
        if (stream->compress) {
            status = nc_def_var_deflate(nc->nc_id, nc->nc_vars[j].nc_varid,
                                        true, true, stream->compress);
            check_nc_status(
                status,
                "Error setting compression level in %s for variable: %s",
                stream->filename, out_metadata[varid].varname);
        }

        // set the fill value attribute
        switch (nc->nc_vars[j].nc_type) {
        case NC_DOUBLE:
            status = nc_put_att_double(nc->nc_id, nc->nc_vars[j].nc_varid,
                                       "_FillValue", NC_DOUBLE, 1,
                                       &(nc->d_fillvalue));
            break;
        case NC_FLOAT:
            status = nc_put_att_float(nc->nc_id, nc->nc_vars[j].nc_varid,
                                      "_FillValue", NC_FLOAT, 1,
                                      &(nc->f_fillvalue));
            break;
        case NC_INT:
            status = nc_put_att_int(nc->nc_id, nc->nc_vars[j].nc_varid,
                                    "_FillValue", NC_INT, 1,
                                    &(nc->i_fillvalue));
            break;
        case NC_SHORT:
            log_err("NC_SHORT not supported yet");
            break;
        case NC_CHAR:
            log_err("NC_CHAR not supported yet");
            break;
        case NC_BYTE:
            log_err("NC_BYTE not supported yet");
            break;
        default:
            log_err("NC_TYPE %d not supported at this time",
                    nc->nc_vars[j].nc_type);
        }
        // check status of fill value setting
        check_nc_status(status,
                        "Error (%d) putting _FillValue attribute to %s in %s",
                        status, out_metadata[varid].varname, stream->filename);

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
    check_nc_status(status, "Error leaving define mode for %s",
                    stream->filename);

    // fill the netcdf variables lat/lon
    if (global_domain.info.n_coord_dims == 1) {
        dvar = calloc(nc->ni_size, sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");

        dcount[0] = nc->ni_size;
        for (i = 0; i < nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].longitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        check_nc_status(status, "Error adding data to lon in %s",
                        stream->filename);
        free(dvar);

        dvar = calloc(nc->nj_size, sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");
        dcount[0] = nc->nj_size;
        for (i = 0; i < nc->nj_size; i++) {
            dvar[i] =
                (double) global_domain.locations[i * nc->ni_size].latitude;
        }

        status =
            nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        check_nc_status(status, "Error adding data to lon in %s",
                        stream->filename);
        free(dvar);
    }
    else if (global_domain.info.n_coord_dims == 2) {
        dvar = calloc(nc->nj_size * nc->ni_size, sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].longitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lon_var_id, dstart, dcount, dvar);
        check_nc_status(status, "Error adding data to lon in %s",
                        stream->filename);

        for (i = 0; i < nc->nj_size * nc->ni_size; i++) {
            dvar[i] = (double) global_domain.locations[i].latitude;
        }
        status =
            nc_put_vara_double(nc->nc_id, lat_var_id, dstart, dcount, dvar);
        check_nc_status(status, "Error adding data to lat in %s",
                        stream->filename);

        free(dvar);
    }
    else {
        log_err("n_coord_dims should be 1 or 2");
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
    char           tmpstr[MAXSTRING];
    char           userstr[MAXSTRING];
    char           hoststr[MAXSTRING];
    char           mpistr[MPI_MAX_LIBRARY_VERSION_STRING];
    int            len;
    int            status;
    time_t         curr_date_time;
    struct tm     *timeinfo;
    uid_t          uid;
    struct passwd *pw;

    // datestr
    curr_date_time = time(NULL);
    if (curr_date_time == -1) {
        log_err("Something went wrong getting the current time!");
    }
    timeinfo = localtime(&curr_date_time);

    // username
    uid = geteuid();
    pw = getpwuid(uid);

    if (pw) {
        strcpy(userstr, pw->pw_name);
    }
    else {
        strcpy(userstr, "unknown");
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
                "Output from the Variable Infiltration Capacity (VIC) "
                "Macroscale Hydrologic Model");
    put_nc_attr(ncid, NC_GLOBAL, "Conventions", "CF-1.6");
    put_nc_attr(ncid, NC_GLOBAL, "netcdf_lib_version", nc_inq_libvers());
    status = MPI_Get_library_version(mpistr, &len);
    if (status == MPI_SUCCESS) {
        put_nc_attr(ncid, NC_GLOBAL, "mpi_lib_version", mpistr);
    }

    // Useful attributes from VIC
    put_nc_attr(ncid, NC_GLOBAL, "VIC_Model_Version", VERSION);
    put_nc_attr(ncid, NC_GLOBAL, "VIC_GIT_VERSION", GIT_VERSION);
    // TODO: pass in driver as an argmument to this function
    put_nc_attr(ncid, NC_GLOBAL, "VIC_Driver", "Image");
}

/******************************************************************************
 * @brief    Set global netcdf attributes (either history or state file)
 *****************************************************************************/
void
initialize_nc_file(nc_file_struct     *nc_file,
                   size_t              nvars,
                   unsigned int       *varids,
                   unsigned short int *dtypes)
{
    extern option_struct options;
    extern domain_struct global_domain;

    size_t               i;

    nc_file->open = false;

    // Set fill values
    nc_file->c_fillvalue = NC_FILL_CHAR;
    nc_file->s_fillvalue = NC_FILL_SHORT;
    nc_file->i_fillvalue = NC_FILL_INT;
    nc_file->d_fillvalue = NC_FILL_DOUBLE;
    nc_file->f_fillvalue = NC_FILL_FLOAT;

    // set ids to MISSING
    nc_file->nc_id = MISSING;
    nc_file->band_dimid = MISSING;
    nc_file->front_dimid = MISSING;
    nc_file->frost_dimid = MISSING;
    nc_file->lake_node_dimid = MISSING;
    nc_file->layer_dimid = MISSING;
    nc_file->ni_dimid = MISSING;
    nc_file->nj_dimid = MISSING;
    nc_file->node_dimid = MISSING;
    nc_file->root_zone_dimid = MISSING;
    nc_file->time_dimid = MISSING;
    nc_file->veg_dimid = MISSING;

    // Set dimension sizes
    nc_file->band_size = options.SNOW_BAND;
    nc_file->front_size = MAX_FRONTS;
    nc_file->frost_size = options.Nfrost;
    nc_file->lake_node_size = MAX_LAKE_NODES;
    nc_file->layer_size = options.Nlayer;
    nc_file->ni_size = global_domain.n_nx;
    nc_file->nj_size = global_domain.n_ny;
    nc_file->node_size = options.Nnode;
    nc_file->root_zone_size = options.ROOT_ZONES;
    nc_file->time_size = NC_UNLIMITED;
    nc_file->veg_size = options.NVEGTYPES;

    // allocate memory for nc_vars
    nc_file->nc_vars = calloc(nvars, sizeof(*(nc_file->nc_vars)));
    check_alloc_status(nc_file->nc_vars, "Memory allocation error.");

    for (i = 0; i < nvars; i++) {
        set_nc_var_info(varids[i], dtypes[i], nc_file, &(nc_file->nc_vars[i]));
    }
}
