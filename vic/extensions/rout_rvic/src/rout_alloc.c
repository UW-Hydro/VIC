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
    extern domain_struct local_domain;
    extern rout_struct   rout;
    // char       *nc_name="./input/stehekin_parameters_01.rvic.prm.Stehekin.20150727.nc";

    // Get some dimensions
    rout.rout_param.iSubsetLength = get_nc_dimension(rout.param_filename,
                                                     "timesteps");
    rout.rout_param.iOutlets = get_nc_dimension(rout.param_filename, "outlets");
    rout.rout_param.iSources = get_nc_dimension(rout.param_filename, "sources");

    // Allocate memory in rout param_struct
    rout.rout_param.source2outlet_ind = (size_t *) malloc(
        rout.rout_param.iSources * sizeof(size_t));
    if (rout.rout_param.source2outlet_ind == NULL) {
        log_err(
            "Memory allocation error in rout.rout_param.source2outlet_ind().");
    }
    rout.rout_param.source_time_offset = (int *) malloc(
        rout.rout_param.iSources * sizeof(int));
    if (rout.rout_param.source_time_offset == NULL) {
        log_err(
            "Memory allocation error in rout.rout_param.source_time_offset().");
    }
    rout.rout_param.source_x_ind = (int *) malloc(
        rout.rout_param.iSources * sizeof(int));
    if (rout.rout_param.source_x_ind == NULL) {
        log_err("Memory allocation error in rout.rout_param.source_x_ind().");
    }
    rout.rout_param.source_y_ind = (int *) malloc(
        rout.rout_param.iSources * sizeof(int));
    if (rout.rout_param.source_y_ind == NULL) {
        log_err("Memory allocation error in rout.rout_param.source_y_ind().");
    }
    rout.rout_param.source_lat = (double *) malloc(
        rout.rout_param.iSources * sizeof(double));
    if (rout.rout_param.source_lat == NULL) {
        log_err("Memory allocation error in rout.rout_param.source_lat().");
    }
    rout.rout_param.source_lon = (double *) malloc(
        rout.rout_param.iSources * sizeof(double));
    if (rout.rout_param.source_lon == NULL) {
        log_err("Memory allocation error in rout.rout_param.source_lon().");
    }
    rout.rout_param.source_VIC_index = (int *) malloc(
        rout.rout_param.iSources * sizeof(int));
    if (rout.rout_param.source_VIC_index == NULL) {
        log_err(
            "Memory allocation error in rout.rout_param.source_VIC_index().");
    }
    rout.rout_param.outlet_lat = (double *) malloc(
        rout.rout_param.iOutlets * sizeof(double));
    if (rout.rout_param.outlet_lat == NULL) {
        log_err("Memory allocation error in rout.rout_param.outlet_lat().");
    }
    rout.rout_param.outlet_lon = (double *) malloc(
        rout.rout_param.iOutlets * sizeof(double));
    if (rout.rout_param.outlet_lon == NULL) {
        log_err("Memory allocation error in rout.rout_param.outlet_lon().");
    }
    rout.rout_param.outlet_VIC_index = (int *) malloc(
        rout.rout_param.iOutlets * sizeof(int));
    if (rout.rout_param.outlet_VIC_index == NULL) {
        log_err(
            "Memory allocation error in rout.rout_param.outlet_VIC_index().");
    }
    rout.rout_param.unit_hydrograph = (double *) malloc(
        rout.rout_param.iSources * rout.rout_param.iSubsetLength *
        sizeof(double));
    if (rout.rout_param.unit_hydrograph == NULL) {
        log_err("Memory allocation error in rout.rout_param.unit_hydrograph().");
    }
    rout.rout_param.aggrunin =
        (double *)malloc(local_domain.ncells_active * sizeof(double));
    if (rout.rout_param.aggrunin == NULL) {
        log_err("Memory allocation error in rout.rout_param.aggrunin().");
    }
    rout.discharge = (double *)malloc(local_domain.ncells_active * sizeof(double));
    if (rout.discharge == NULL) {
        log_err("Memory allocation error in rout.rout_param.discharge().");
    }

    // Allocate memory for the ring
    rout.ring = (double *)malloc(
        rout.rout_param.iSubsetLength * rout.rout_param.iOutlets *
        sizeof(double));
    if (rout.ring == NULL) {
        log_err("Memory allocation error in rout.ring().");
    }
}
