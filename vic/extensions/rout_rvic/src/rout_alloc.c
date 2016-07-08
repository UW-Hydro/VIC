/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for Routing structures.
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
#include <rout.h>

/******************************************************************************
 * @brief    Allocate memory for Routing structures.
 *****************************************************************************/
void
rout_alloc(void)
{
    extern int mpi_rank;
    if (mpi_rank == VIC_MPI_ROOT) {
        extern domain_struct global_domain;
        extern rout_struct   rout;
        int                  ivar;
        size_t               d1count[1];
        size_t               d1start[1];

        d1count[0] = 0;
        d1start[0] = 1;

        // Get some values and dimensions
        get_nc_field_int(rout.param_filename,
                         "full_time_length",
                         d1start,
                         d1count,
                         &ivar);
        rout.rout_param.full_time_length = (int) ivar;

        rout.rout_param.nTimesteps = get_nc_dimension(rout.param_filename,
                                                      "timesteps");
        rout.rout_param.nOutlets = get_nc_dimension(rout.param_filename,
                                                    "outlets");
        rout.rout_param.nSources = get_nc_dimension(rout.param_filename,
                                                    "sources");

        // Allocate memory in rout param_struct
        rout.rout_param.source2outlet_ind = malloc(
            rout.rout_param.nSources *
            sizeof(*rout.rout_param.source2outlet_ind));
        if (rout.rout_param.source2outlet_ind == NULL) {
            log_err(
                "Memory allocation error in rout.rout_param.source2outlet_ind().");
        }
        rout.rout_param.source_time_offset = malloc(
            rout.rout_param.nSources *
            sizeof(*rout.rout_param.source_time_offset));
        if (rout.rout_param.source_time_offset == NULL) {
            log_err(
                "Memory allocation error in rout.rout_param.source_time_offset().");
        }
        rout.rout_param.source_x_ind = malloc(
            rout.rout_param.nSources * sizeof(*rout.rout_param.source_x_ind));
        if (rout.rout_param.source_x_ind == NULL) {
            log_err("Memory allocation error in rout.rout_param.source_x_ind().");
        }
        rout.rout_param.source_y_ind = malloc(
            rout.rout_param.nSources * sizeof(*rout.rout_param.source_y_ind));
        if (rout.rout_param.source_y_ind == NULL) {
            log_err("Memory allocation error in rout.rout_param.source_y_ind().");
        }
        rout.rout_param.source_lat = malloc(
            rout.rout_param.nSources * sizeof(*rout.rout_param.source_lat));
        if (rout.rout_param.source_lat == NULL) {
            log_err("Memory allocation error in rout.rout_param.source_lat().");
        }
        rout.rout_param.source_lon = malloc(
            rout.rout_param.nSources * sizeof(*rout.rout_param.source_lon));
        if (rout.rout_param.source_lon == NULL) {
            log_err("Memory allocation error in rout.rout_param.source_lon().");
        }
        rout.rout_param.source_VIC_index = malloc(
            rout.rout_param.nSources *
            sizeof(*rout.rout_param.source_VIC_index));
        if (rout.rout_param.source_VIC_index == NULL) {
            log_err(
                "Memory allocation error in rout.rout_param.source_VIC_index().");
        }
        rout.rout_param.outlet_lat = malloc(
            rout.rout_param.nOutlets * sizeof(*rout.rout_param.outlet_lat));
        if (rout.rout_param.outlet_lat == NULL) {
            log_err("Memory allocation error in rout.rout_param.outlet_lat().");
        }
        rout.rout_param.outlet_lon = malloc(
            rout.rout_param.nOutlets * sizeof(*rout.rout_param.outlet_lon));
        if (rout.rout_param.outlet_lon == NULL) {
            log_err("Memory allocation error in rout.rout_param.outlet_lon().");
        }
        rout.rout_param.outlet_VIC_index = malloc(
            rout.rout_param.nOutlets *
            sizeof(rout.rout_param.outlet_VIC_index));
        if (rout.rout_param.outlet_VIC_index == NULL) {
            log_err(
                "Memory allocation error in rout.rout_param.outlet_VIC_index().");
        }
        rout.rout_param.unit_hydrograph = malloc(
            rout.rout_param.nSources * rout.rout_param.nTimesteps *
            sizeof(*rout.rout_param.unit_hydrograph));
        if (rout.rout_param.unit_hydrograph == NULL) {
            log_err(
                "Memory allocation error in rout.rout_param.unit_hydrograph().");
        }
        rout.rout_param.aggrunin =
            malloc(
                global_domain.ncells_total * sizeof(*rout.rout_param.aggrunin));
        if (rout.rout_param.aggrunin == NULL) {
            log_err("Memory allocation error in rout.rout_param.aggrunin().");
        }
        rout.discharge =
            malloc(global_domain.ncells_total * sizeof(*rout.discharge));
        if (rout.discharge == NULL) {
            log_err("Memory allocation error in rout.rout_param.discharge().");
        }

        // Allocate memory for the ring
        rout.ring = malloc(
            rout.rout_param.full_time_length * rout.rout_param.nOutlets *
            sizeof(*rout.ring));
        if (rout.ring == NULL) {
            log_err("Memory allocation error in rout.ring().");
        }
    }
}
