/******************************************************************************
 * @section DESCRIPTION
 *
 * Write output to netcdf file
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
 * @brief    Write output data and convert units if necessary.
 *****************************************************************************/
void
vic_write_output(dmy_struct *dmy)
{
    extern option_struct   options;
    extern stream_struct  *output_streams;
    extern nc_file_struct *nc_hist_files;

    size_t                 stream_idx;

    // Write data
    for (stream_idx = 0; stream_idx < options.Noutstreams; stream_idx++) {
        if (raise_alarm(&(output_streams[stream_idx].agg_alarm), dmy)) {
            debug("raised alarm for stream %zu", stream_idx);
            vic_write(&(output_streams[stream_idx]),
                      &(nc_hist_files[stream_idx]), dmy);
            reset_stream(&(output_streams[stream_idx]), dmy);
        }
    }
}

/******************************************************************************
 * @brief    Write output to netcdf file. Currently everything is cast to
 *           double
 *****************************************************************************/
void
vic_write(stream_struct  *stream,
          nc_file_struct *nc_hist_file,
          dmy_struct     *dmy_current)
{
    extern global_param_struct global_param;
    extern domain_struct       local_domain;
    extern int                 mpi_rank;
    extern metadata_struct     out_metadata[N_OUTVAR_TYPES];

    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     ndims;
    double                     dtime;
    double                    *dvar = NULL;
    float                     *fvar = NULL;
    int                       *ivar = NULL;
    short int                 *svar = NULL;
    char                      *cvar = NULL;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    unsigned int               varid;
    int                        status;
    double                     offset;
    double                     bounds[2];

    if (mpi_rank == VIC_MPI_ROOT) {
        // If the output file is not open, initialize the history file now.
        if (nc_hist_file->open == false) {
            // open the netcdf history file
            initialize_history_file(nc_hist_file, stream, dmy_current);
        }
    }

    // initialize dimids to invalid values - helps debugging
    for (i = 0; i < MAXDIMS; i++) {
        dstart[i] = -1;
        dcount[i] = 0;
    }

    for (k = 0; k < stream->nvars; k++) {
        varid = stream->varid[k];

        if (nc_hist_file->nc_vars[k].nc_type == NC_DOUBLE) {
            if (dvar == NULL) {
                // allocate memory for variables to be stored
                dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
                check_alloc_status(dvar, "Memory allocation error");
            }
        }
        else if (nc_hist_file->nc_vars[k].nc_type == NC_FLOAT) {
            if (fvar == NULL) {
                // allocate memory for variables to be stored
                fvar = malloc(local_domain.ncells_active * sizeof(*fvar));
                check_alloc_status(fvar, "Memory allocation error");
            }
        }
        else if (nc_hist_file->nc_vars[k].nc_type == NC_INT) {
            if (ivar == NULL) {
                // allocate memory for variables to be stored
                ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
                check_alloc_status(ivar, "Memory allocation error");
            }
        }
        else if (nc_hist_file->nc_vars[k].nc_type == NC_SHORT) {
            if (svar == NULL) {
                // allocate memory for variables to be stored
                svar = malloc(local_domain.ncells_active * sizeof(*svar));
                check_alloc_status(svar, "Memory allocation error");
            }
        }
        else if (nc_hist_file->nc_vars[k].nc_type == NC_CHAR) {
            if (cvar == NULL) {
                // allocate memory for variables to be stored
                cvar = malloc(local_domain.ncells_active * sizeof(*cvar));
                check_alloc_status(cvar, "Memory allocation error");
            }
        }
        else {
            log_err("Unsupported nc_type encountered");
        }

        ndims = nc_hist_file->nc_vars[k].nc_dims;
        for (j = 0; j < ndims; j++) {
            dstart[j] = 0;
            dcount[j] = 1;
        }
        // The size of the last two dimensions are the grid size; files are
        // written one slice at a time, so all counts are 1, except the last
        // two
        for (j = ndims - 2; j < ndims; j++) {
            dcount[j] = nc_hist_file->nc_vars[k].nc_counts[j];
        }
        dstart[0] = stream->write_alarm.count;  // Position in the time dimensions

        for (j = 0; j < out_metadata[varid].nelem; j++) {
            // if there is more than one layer, then dstart needs to advance
            dstart[1] = j;
            if (nc_hist_file->nc_vars[k].nc_type == NC_DOUBLE) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = (double) stream->aggdata[i][k][j][0];
                }
                gather_put_nc_field_double(nc_hist_file->nc_id,
                                           nc_hist_file->nc_vars[k].nc_varid,
                                           nc_hist_file->d_fillvalue,
                                           dstart, dcount, dvar);
            }
            else if (nc_hist_file->nc_vars[k].nc_type == NC_FLOAT) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    fvar[i] = (float) stream->aggdata[i][k][j][0];
                }
                gather_put_nc_field_float(nc_hist_file->nc_id,
                                          nc_hist_file->nc_vars[k].nc_varid,
                                          nc_hist_file->f_fillvalue,
                                          dstart, dcount, fvar);
            }
            else if (nc_hist_file->nc_vars[k].nc_type == NC_INT) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    ivar[i] = (int) stream->aggdata[i][k][j][0];
                }
                gather_put_nc_field_int(nc_hist_file->nc_id,
                                        nc_hist_file->nc_vars[k].nc_varid,
                                        nc_hist_file->i_fillvalue,
                                        dstart, dcount, ivar);
            }
            else if (nc_hist_file->nc_vars[k].nc_type == NC_SHORT) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    svar[i] = (short int) stream->aggdata[i][k][j][0];
                }
                gather_put_nc_field_short(nc_hist_file->nc_id,
                                          nc_hist_file->nc_vars[k].nc_varid,
                                          nc_hist_file->s_fillvalue,
                                          dstart, dcount, svar);
            }
            else if (nc_hist_file->nc_vars[k].nc_type == NC_CHAR) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    cvar[i] = (char) stream->aggdata[i][k][j][0];
                }
                gather_put_nc_field_schar(nc_hist_file->nc_id,
                                          nc_hist_file->nc_vars[k].nc_varid,
                                          nc_hist_file->d_fillvalue,
                                          dstart, dcount, cvar);
            }
            else {
                log_err("Unsupported nc_type encountered");
            }
        }

        // reset dimids to invalid values - helps debugging
        for (j = 0; j < MAXDIMS; j++) {
            dstart[j] = -1;
            dcount[j] = -1;
        }
    }

    // write to file
    if (mpi_rank == VIC_MPI_ROOT) {
        // Add time variable
        dstart[0] = stream->write_alarm.count;

        // timestamp is the beginning of the aggregation window
        dtime = date2num(global_param.time_origin_num,
                         &(stream->time_bounds[0]), 0.,
                         global_param.calendar, global_param.time_units);

        status = nc_put_var1_double(nc_hist_file->nc_id,
                                    nc_hist_file->time_varid,
                                    dstart, &dtime);
        check_nc_status(status, "Error writing time variable");

        // Add time bounds variable
        dstart[1] = 0;
        dcount[0] = 1;
        dcount[1] = 2;
        bounds[0] = dtime;
        dt_seconds_to_time_units(global_param.time_units, global_param.dt,
                                 &offset);
        bounds[1] = offset + date2num(global_param.time_origin_num,
                                      &(stream->time_bounds[1]), 0.,
                                      global_param.calendar,
                                      global_param.time_units);

        status = nc_put_vara_double(nc_hist_file->nc_id,
                                    nc_hist_file->time_bounds_varid,
                                    dstart, dcount, bounds);
        check_nc_status(status, "Error writing time bounds variable");
    }

    // Advance the position in the history file
    stream->write_alarm.count++;
    if (raise_alarm(&(stream->write_alarm), dmy_current)) {
        // close this history file
        if (mpi_rank == VIC_MPI_ROOT) {
            status = nc_close(nc_hist_file->nc_id);
            check_nc_status(status, "Error closing history file");
            nc_hist_file->open = false;
        }
        reset_alarm(&(stream->write_alarm), dmy_current);
    }
    else {
        // Force sync with disk (GH:#596)
        if (mpi_rank == VIC_MPI_ROOT) {
            status = nc_sync(nc_hist_file->nc_id);
            check_nc_status(status, "Error syncing netCDF file %s",
                            stream->filename);
        }
    }

    // free memory
    if (dvar != NULL) {
        free(dvar);
    }
    if (fvar != NULL) {
        free(fvar);
    }
    if (ivar != NULL) {
        free(ivar);
    }
    if (svar != NULL) {
        free(svar);
    }
    if (cvar != NULL) {
        free(cvar);
    }
}
