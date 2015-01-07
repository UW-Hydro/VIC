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
    extern int                 mpi_rank;
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

        // initialize netcdf info for output variables
        vic_nc_info(&nc_hist_file, out_data, nc_vars);
    }
}

/******************************************************************************
 * @brief    Initialize history files
 *****************************************************************************/
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
        log_ncerr(status);
    }
    nc->open = true;

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "front", nc->front_size,
                        &(nc->front_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "ni", nc->ni_size, &(nc->ni_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "nj", nc->nj_size, &(nc->nj_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }

    status = nc_def_dim(nc->nc_id, "time", nc->time_size,
                        &(nc->time_dimid));
    if (status != NC_NOERR) {
        log_ncerr(status);
    }


    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        log_ncerr(status);
    }
}
