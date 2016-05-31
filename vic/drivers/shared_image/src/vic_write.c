/******************************************************************************
 * @section DESCRIPTION
 *
 * Write output to netcdf file
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
            // reset_nc_file(&(nc_hist_files[stream_idx]));
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
    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    int                        dimids[MAXDIMS];
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     ndims;
    double                    *dvar = NULL;
    size_t                     dcount[MAXDIMS];
    size_t                     dstart[MAXDIMS];
    unsigned int               varid;

    // Advance the position in the history file
    stream->write_alarm.count++;

    // allocate memory for variables to be stored
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_write().");
    }

    // set missing values
    for (i = 0; i < local_domain.ncells_active; i++) {
        dvar[i] = nc_hist_file->d_fillvalue;
    }

    // initialize dimids to invalid values - helps debugging
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
        dstart[i] = -1;
        dcount[i] = -1;
    }

    for (k = 0; k < stream->nvars; k++) {
        varid = stream->varid[k];

        ndims = nc_hist_file->nc_vars[k].nc_dims;
        for (j = 0; j < ndims; j++) {
            dimids[j] = nc_hist_file->nc_vars[k].nc_dimids[j];
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
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = (double) stream->aggdata[i][k][j][0];
            }
            gather_put_nc_field_double(stream->filename, &(nc_hist_file->open),
                                       &(nc_hist_file->nc_id),
                                       nc_hist_file->d_fillvalue,
                                       dimids, ndims,
                                       out_metadata[varid].varname,
                                       dstart, dcount, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_hist_file->d_fillvalue;
            }
        }

        // reset dimids to invalid values - helps debugging
        for (j = 0; j < MAXDIMS; j++) {
            dimids[j] = -1;
            dstart[j] = -1;
            dcount[j] = -1;
        }
    }

    // write to file
    if (mpi_rank == 0) {
        // ADD Time variable
        dimids[0] = nc_hist_file->time_dimid;
        dstart[0] = stream->write_alarm.count;
        dcount[0] = 1;

        dvar[0] = date2num(global_param.time_origin_num, dmy_current, 0.,
                           global_param.calendar, global_param.time_units);


        put_nc_field_double(stream->filename, &(nc_hist_file->open),
                            &(nc_hist_file->nc_id),
                            nc_hist_file->d_fillvalue,
                            dimids, 1, "time",
                            dstart, dcount, dvar);
    }

    // free memory
    free(dvar);
}
