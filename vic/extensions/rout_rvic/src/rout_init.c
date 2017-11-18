/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize routing model parameters
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

#include <rout.h>

/******************************************************************************
 * @brief    Initialize routing model parameters
 *****************************************************************************/
void
rout_init(void)
{
    extern int              mpi_rank;
    extern rout_struct      rout;
    extern domain_struct    global_domain;
    extern filenames_struct filenames;
    int                     status;

    if (mpi_rank == VIC_MPI_ROOT) {
        int    *ivar = NULL;
        double *dvar = NULL;

        size_t  i;
        size_t  j;
        size_t  i1start;
        size_t  d3count[3];
        size_t  d3start[3];

        i1start = 0;

        d3start[0] = 0;
        d3start[1] = 0;
        d3start[2] = 0;
        d3count[0] = rout.rout_param.n_timesteps;
        d3count[1] = rout.rout_param.n_sources;
        d3count[2] = 1; // tracers dimension

        // allocate memory for variables to be read
        ivar = malloc(rout.rout_param.n_sources * sizeof(*ivar));
        check_alloc_status(ivar, "Memory allocation error.");

        // allocate memory for variables to be read
        dvar = malloc(
            rout.rout_param.n_timesteps * rout.rout_param.n_sources *
            sizeof(*dvar));
        check_alloc_status(dvar, "Memory allocation error.");

        // The Ring
        for (j = 0; j < rout.rout_param.n_outlets; j++) {
            for (i = 0; i < rout.rout_param.n_timesteps; i++) {
                rout.ring[j * rout.rout_param.n_timesteps + i] = 0.0;
            }
        }

        // discharge
        for (j = 0; j < global_domain.ncells_active; j++) {
            rout.discharge[j] = 0.0;
        }

        // source2outlet_ind: source to outlet index mapping
        get_nc_field_int(&(filenames.rout_params),
                         "source2outlet_ind",
                         &i1start, &rout.rout_param.n_sources, ivar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source2outlet_ind[i] = (int) ivar[i];
        }

        // source_time_offset: Number of leading timesteps ommited
        get_nc_field_int(&(filenames.rout_params),
                         "source_time_offset",
                         &i1start, &rout.rout_param.n_sources, ivar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source_time_offset[i] = (int) ivar[i];
        }

        // source_x_ind: x grid coordinate of source grid cell
        get_nc_field_int(&(filenames.rout_params),
                         "source_x_ind",
                         &i1start, &rout.rout_param.n_sources, ivar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source_x_ind[i] = (int) ivar[i];
        }

        // source_y_ind: y grid coordinate of source grid cell
        get_nc_field_int(&(filenames.rout_params),
                         "source_y_ind",
                         &i1start, &rout.rout_param.n_sources, ivar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source_y_ind[i] = (int) ivar[i];
        }

        // source_lat: Latitude coordinate of source grid cell
        get_nc_field_double(&(filenames.rout_params),
                            "source_lat",
                            &i1start, &rout.rout_param.n_sources, dvar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source_lat[i] = (double) dvar[i];
        }

        // source_lon: Longitude coordinate of source grid cell
        get_nc_field_double(&(filenames.rout_params),
                            "source_lon",
                            &i1start, &rout.rout_param.n_sources, dvar);
        for (i = 0; i < rout.rout_param.n_sources; i++) {
            rout.rout_param.source_lon[i] = (double) dvar[i];
        }

        // outlet_lat: Latitude coordinate of source grid cell
        get_nc_field_double(&(filenames.rout_params),
                            "outlet_lat",
                            &i1start, &rout.rout_param.n_outlets, dvar);
        for (i = 0; i < rout.rout_param.n_outlets; i++) {
            rout.rout_param.outlet_lat[i] = (double) dvar[i];
        }

        // outlet_lon: Longitude coordinate of source grid cell
        get_nc_field_double(&(filenames.rout_params),
                            "outlet_lon",
                            &i1start, &rout.rout_param.n_outlets, dvar);
        for (i = 0; i < rout.rout_param.n_outlets; i++) {
            rout.rout_param.outlet_lon[i] = (double) dvar[i];
        }

        // Unit Hydrograph:
        get_nc_field_double(&(filenames.rout_params),
                            "unit_hydrograph",
                            d3start, d3count, dvar);
        for (i = 0;
             i < (rout.rout_param.n_timesteps * rout.rout_param.n_sources);
             i++) {
            rout.rout_param.unit_hydrograph[i] = (double) dvar[i];
        }

        // TODO: add check: what to do in case no VIC gridcell exists for a rout source?
        // Mapping: Let the routing-source index numbers correspond to the VIC index numbers
        size_t i_source;
        for (i_source = 0; i_source < rout.rout_param.n_sources; i_source++) {
            for (i = 0; i < global_domain.ncells_total; i++) {
                if (rout.rout_param.source_lat[i_source] ==
                    global_domain.locations[i].latitude &&
                    rout.rout_param.source_lon[i_source] ==
                    global_domain.locations[i].longitude) {
                    rout.rout_param.source_VIC_index[i_source] = i;
                }
            }
        }

        // Check source index of VIC gridcell
        for (i_source = 0; i_source < rout.rout_param.n_sources; i_source++) {
            if ((size_t)rout.rout_param.source_VIC_index[i_source] < 0 ||
                (size_t)rout.rout_param.source_VIC_index[i_source] >
                global_domain.ncells_total) {
                log_err("invalid source, index of VIC gridcell");
            }
        }

        // Mapping: Let i_outlet the routing-outlet index numbers correspond to the VIC index numbers
        size_t i_outlet;
        for (i_outlet = 0; i_outlet < rout.rout_param.n_outlets; i_outlet++) {
            for (i = 0; i < global_domain.ncells_total; i++) {
                if (rout.rout_param.outlet_lat[i_outlet] ==
                    global_domain.locations[i].latitude &&
                    rout.rout_param.outlet_lon[i_outlet] ==
                    global_domain.locations[i].longitude) {
                    rout.rout_param.outlet_VIC_index[i_outlet] = i;
                }
            }
        }

        // Check outlet index of VIC gridcell
        for (i_outlet = 0; i_outlet < rout.rout_param.n_outlets; i_outlet++) {
            if ((size_t)rout.rout_param.outlet_VIC_index[i_outlet] < 0 ||
                (size_t)rout.rout_param.outlet_VIC_index[i_outlet] >
                global_domain.ncells_total) {
                log_err("invalid outlet, index of VIC gridcell");
            }
        }

        // close parameter file
        status = nc_close(filenames.rout_params.nc_id);
        check_nc_status(status, "Error closing %s",
                        filenames.rout_params.nc_filename);

        // cleanup
        free(ivar);
        free(dvar);
    }
}
