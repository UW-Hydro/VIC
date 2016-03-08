/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine handles the startup tasks.
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

#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Wrapper function for VIC startup tasks.
 *****************************************************************************/
void
vic_start(void)
{
    int                        local_ncells_active;
    int                        status;
    location_struct           *mapped_locations = NULL;
    location_struct           *active_locations = NULL;
    size_t                     i;
    extern size_t             *filter_active_cells;
    extern size_t             *mpi_map_mapping_array;
    extern filenames_struct    filenames;
    extern filep_struct        filep;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern global_param_struct global_param;
    extern MPI_Datatype        mpi_global_struct_type;
    extern MPI_Datatype        mpi_filenames_struct_type;
    extern MPI_Datatype        mpi_location_struct_type;
    extern MPI_Datatype        mpi_option_struct_type;
    extern MPI_Datatype        mpi_param_struct_type;
    extern int                *mpi_map_local_array_sizes;
    extern int                *mpi_map_global_array_offsets;
    extern int                 mpi_rank;
    extern int                 mpi_size;
    extern option_struct       options;
    extern parameters_struct   param;
    size_t                     j;

    // Initialize global structures
    initialize_options();
    initialize_global();
    initialize_parameters();
    initialize_filenames();

    if (mpi_rank == 0) {
        // Initialize the global domain
        initialize_domain(&global_domain);

        // read global settings
        filep.globalparam = open_file(filenames.global, "r");
        get_global_param(filep.globalparam);
    }

    status = MPI_Bcast(&filenames, 1, mpi_filenames_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    // Set Log Destination
    setup_logging(mpi_rank);

    if (mpi_rank == 0) {
        // set model constants
        if (strcasecmp(filenames.constants, "MISSING")) {
            filep.constants = open_file(filenames.constants, "r");
            get_parameters(filep.constants);
        }

        // read domain info
        get_global_domain(filenames.domain, &global_domain);

        // add the number of vegetation type to the location info in the
        // global domain struct. This just makes life easier
        add_nveg_to_global_domain(filenames.veglib, &global_domain);

        // decompose the mask
        mpi_map_decomp_domain(global_domain.ncells_active, mpi_size,
                              &mpi_map_local_array_sizes,
                              &mpi_map_global_array_offsets,
                              &mpi_map_mapping_array);

        // get the indices for the active cells (used in reading and writing)
        filter_active_cells = malloc(global_domain.ncells_active *
                                     sizeof(*filter_active_cells));

        j = 0;
        for (i = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                filter_active_cells[j] = global_domain.locations[i].io_idx;
                j++;
            }
        }

        // get dimensions (number of vegetation types, soil zones, etc)
        options.ROOT_ZONES = get_nc_dimension(filenames.soil, "root_zone");
        options.Nlayer = get_nc_dimension(filenames.soil, "nlayer");
        options.NVEGTYPES = get_nc_dimension(filenames.veg, "veg_class");
        if (options.SNOW_BAND > 1) {
            if (options.SNOW_BAND !=
                get_nc_dimension(filenames.snowband, "snow_band")) {
                log_err("Number of snow bands in global file does not "
                        "match parameter file");
            }
        }
        if (options.LAKES) {
            options.NLAKENODES = get_nc_dimension(filenames.lakeparam,
                                                  "lake_node");
        }

        // Check that model parameters are valid
        validate_parameters();
    }

    // broadcast global, option, param structures as well as global valies
    // such as NF and NR
    status = MPI_Bcast(&NF, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    status = MPI_Bcast(&NR, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    status = MPI_Bcast(&global_param, 1, mpi_global_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    status = MPI_Bcast(&options, 1, mpi_option_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    status = MPI_Bcast(&param, 1, mpi_param_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    // setup the local domain_structs

    // First scatter the array sizes
    status = MPI_Scatter(mpi_map_local_array_sizes, 1, MPI_INT,
                         &local_ncells_active, 1, MPI_INT, 0, MPI_COMM_WORLD);
    local_domain.ncells_active = (size_t) local_ncells_active;
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }

    // Allocate memory for the local locations
    local_domain.locations = malloc(local_domain.ncells_active *
                                    sizeof(*local_domain.locations));
    if (local_domain.locations == NULL) {
        log_err("malloc error in vic_start()\n");
    }
    for (i = 0; i < local_domain.ncells_active; i++) {
        initialize_location(&(local_domain.locations[i]));
    }

    // map the location vector to a temporary array so they can be scattered
    if (mpi_rank == 0) {
        mapped_locations = malloc(global_domain.ncells_active *
                                  sizeof(*mapped_locations));
        if (mapped_locations == NULL) {
            log_err("malloc error in vic_start()\n");
        }
        for (i = 0; i < global_domain.ncells_active; i++) {
            initialize_location(&(mapped_locations[i]));
        }

        active_locations = (location_struct *) malloc(
            global_domain.ncells_active * sizeof(location_struct));
        if (active_locations == NULL) {
            log_err("malloc error in vic_start()\n");
        }
        for (i = 0; i < global_domain.ncells_active; i++) {
            initialize_location(&(active_locations[i]));
        }

        for (i = 0, j = 0; i < global_domain.ncells_total; i++) {
            if (global_domain.locations[i].run) {
                active_locations[j] = global_domain.locations[i];
                j++;
            }
        }

        map(sizeof(location_struct), global_domain.ncells_active,
            mpi_map_mapping_array, NULL, active_locations,
            mapped_locations);
    }

    // Scatter the locations
    status = MPI_Scatterv(mapped_locations, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets,
                          mpi_location_struct_type,
                          local_domain.locations, local_domain.ncells_active,
                          mpi_location_struct_type,
                          0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in vic_start(): %d\n", status);
    }
    // Set the local index value
    for (i = 0; i < (size_t) local_domain.ncells_active; i++) {
        local_domain.locations[i].local_idx = i;
    }

    // cleanup
    if (mpi_rank == 0) {
        free(mapped_locations);
        free(active_locations);
    }
}
