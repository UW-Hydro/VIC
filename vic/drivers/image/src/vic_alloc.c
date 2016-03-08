/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
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
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void
vic_alloc(void)
{
    extern all_vars_struct    *all_vars;
    extern atmos_data_struct  *atmos;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern out_data_struct   **out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;
    extern lake_con_struct    *lake_con;
    size_t                     i;
    size_t                     j;

    // allocate memory for atmos structure
    atmos = malloc(local_domain.ncells_active * sizeof(*atmos));
    if (atmos == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for veg_hist structure
    veg_hist = malloc(local_domain.ncells_active * sizeof(*veg_hist));
    if (veg_hist == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for soil structure
    soil_con = malloc(local_domain.ncells_active * sizeof(*soil_con));
    if (soil_con == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation mapping structure
    veg_con_map = malloc(local_domain.ncells_active * sizeof(*veg_con_map));
    if (veg_con_map == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_con = malloc(local_domain.ncells_active * sizeof(*veg_con));
    if (veg_con == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_lib = malloc(local_domain.ncells_active * sizeof(*veg_lib));
    if (veg_lib == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    if (options.LAKES) {
        // allocate memory for lake structure
        lake_con = malloc(local_domain.ncells_active * sizeof(*lake_con));
        if (lake_con == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
    }

    // all_vars allocation
    all_vars = malloc(local_domain.ncells_active * sizeof(*all_vars));
    if (all_vars == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // out_data allocation
    out_data = malloc(local_domain.ncells_active * sizeof(*out_data));
    if (out_data == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // save_data allocation
    save_data = malloc(local_domain.ncells_active * sizeof(*save_data));
    if (save_data == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for individual grid cells
    for (i = 0; i < local_domain.ncells_active; i++) {
        // atmos allocation - allocate enough memory for NR+1 steps
        alloc_atmos(&(atmos[i]));

        // snow band allocation
        soil_con[i].AreaFract = calloc(options.SNOW_BAND,
                                       sizeof(*(soil_con[i].AreaFract)));
        if (soil_con[i].AreaFract == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].BandElev = calloc(options.SNOW_BAND,
                                      sizeof(*(soil_con[i].BandElev)));
        if (soil_con[i].BandElev == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Tfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Tfactor)));
        if (soil_con[i].Tfactor == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Pfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Pfactor)));
        if (soil_con[i].Pfactor == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].AboveTreeLine = calloc(options.SNOW_BAND,
                                           sizeof(*(soil_con[i].AboveTreeLine)));
        if (soil_con[i].AboveTreeLine == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        initialize_soil_con(&(soil_con[i]));

        // vegetation tile allocation

        veg_con_map[i].nv_types = options.NVEGTYPES;

        veg_con_map[i].vidx = calloc(veg_con_map[i].nv_types,
                                     sizeof(*(veg_con_map[i].vidx)));
        if (veg_con_map[i].vidx == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        veg_con_map[i].Cv = calloc(veg_con_map[i].nv_types,
                                   sizeof(*(veg_con_map[i].Cv)));
        if (veg_con_map[i].Cv == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        veg_con_map[i].nv_active = (size_t) local_domain.locations[i].nveg + 1;
        if (options.AboveTreelineVeg >= 0) {
            veg_con_map[i].nv_active += 1;
        }

        veg_con[i] = malloc((veg_con_map[i].nv_active) * sizeof(*(veg_con[i])));
        if (veg_con[i] == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].zone_depth = calloc(options.ROOT_ZONES,
                                              sizeof(*(veg_con[i][j].zone_depth)));
            if (veg_con[i][j].zone_depth == NULL) {
                log_err("Memory allocation error in vic_alloc().");
            }
            veg_con[i][j].zone_fract = calloc(options.ROOT_ZONES,
                                              sizeof(*(veg_con[i][j].zone_fract)));
            if (veg_con[i][j].zone_fract == NULL) {
                log_err("Memory allocation error in vic_alloc().");
            }
            if (options.CARBON) {
                veg_con[i][j].CanopLayerBnd = calloc(options.Ncanopy,
                                                     sizeof(*(veg_con[i][j].
                                                              CanopLayerBnd)));
                if (veg_con[i][j].CanopLayerBnd == NULL) {
                    log_err("Memory allocation error in vic_alloc().");
                }
            }
            initialize_veg_con(&(veg_con[i][j]));
        }

        // vegetation library allocation - there is a veg library for each
        // active grid cell

        veg_lib[i] = calloc(options.NVEGTYPES, sizeof(*(veg_lib[i])));
        if (veg_lib[i] == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        all_vars[i] = make_all_vars(veg_con_map[i].nv_active);

        out_data[i] = create_output_list();

        // allocate memory for veg_hist
        veg_hist[i] = calloc(veg_con_map[i].nv_active, sizeof(*(veg_hist[i])));
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            alloc_veg_hist(&(veg_hist[i][j]));
        }

    }
}
