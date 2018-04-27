/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
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
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void
vic_alloc(void)
{
    extern all_vars_struct    *all_vars;
    extern force_data_struct  *force;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern double           ***out_data;
    extern save_data_struct   *save_data;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;
    extern lake_con_struct    *lake_con;
    size_t                     i;
    size_t                     j;

    // allocate memory for force structure
    force = malloc(local_domain.ncells_active * sizeof(*force));
    check_alloc_status(force, "Memory allocation error.");

    // allocate memory for veg_hist structure
    veg_hist = malloc(local_domain.ncells_active * sizeof(*veg_hist));
    check_alloc_status(veg_hist, "Memory allocation error.");

    // allocate memory for soil structure
    soil_con = malloc(local_domain.ncells_active * sizeof(*soil_con));
    check_alloc_status(soil_con, "Memory allocation error.");

    // allocate memory for vegetation mapping structure
    veg_con_map = malloc(local_domain.ncells_active * sizeof(*veg_con_map));
    check_alloc_status(veg_con_map, "Memory allocation error.");

    // allocate memory for vegetation structure
    veg_con = malloc(local_domain.ncells_active * sizeof(*veg_con));
    check_alloc_status(veg_con, "Memory allocation error.");

    // allocate memory for vegetation structure
    veg_lib = malloc(local_domain.ncells_active * sizeof(*veg_lib));
    check_alloc_status(veg_lib, "Memory allocation error.");

    if (options.LAKES) {
        // allocate memory for lake structure
        lake_con = malloc(local_domain.ncells_active * sizeof(*lake_con));
        check_alloc_status(lake_con, "Memory allocation error.");
    }

    // all_vars allocation
    all_vars = malloc(local_domain.ncells_active * sizeof(*all_vars));
    check_alloc_status(all_vars, "Memory allocation error.");

    // out_data allocation
    out_data = malloc(local_domain.ncells_active * sizeof(*out_data));
    check_alloc_status(out_data, "Memory allocation error.");

    // save_data allocation
    save_data = calloc(local_domain.ncells_active, sizeof(*save_data));
    check_alloc_status(save_data, "Memory allocation error.");

    // allocate memory for individual grid cells
    for (i = 0; i < local_domain.ncells_active; i++) {
        // force allocation - allocate enough memory for NR+1 steps
        alloc_force(&(force[i]));

        // snow band allocation
        soil_con[i].AreaFract = calloc(options.SNOW_BAND,
                                       sizeof(*(soil_con[i].AreaFract)));
        check_alloc_status(soil_con[i].AreaFract, "Memory allocation error.");
        soil_con[i].BandElev = calloc(options.SNOW_BAND,
                                      sizeof(*(soil_con[i].BandElev)));
        check_alloc_status(soil_con[i].BandElev, "Memory allocation error.");
        soil_con[i].Tfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Tfactor)));
        check_alloc_status(soil_con[i].Tfactor, "Memory allocation error.");
        soil_con[i].Pfactor = calloc(options.SNOW_BAND,
                                     sizeof(*(soil_con[i].Pfactor)));
        check_alloc_status(soil_con[i].Pfactor, "Memory allocation error.");
        soil_con[i].AboveTreeLine = calloc(options.SNOW_BAND,
                                           sizeof(*(soil_con[i].AboveTreeLine)));
        check_alloc_status(soil_con[i].AboveTreeLine,
                           "Memory allocation error.");

        initialize_soil_con(&(soil_con[i]));

        // vegetation tile allocation

        veg_con_map[i].nv_types = options.NVEGTYPES;

        veg_con_map[i].vidx = calloc(veg_con_map[i].nv_types,
                                     sizeof(*(veg_con_map[i].vidx)));
        check_alloc_status(veg_con_map[i].vidx, "Memory allocation error.");
        veg_con_map[i].Cv = calloc(veg_con_map[i].nv_types,
                                   sizeof(*(veg_con_map[i].Cv)));
        check_alloc_status(veg_con_map[i].Cv, "Memory allocation error.");

        veg_con_map[i].nv_active = (size_t) local_domain.locations[i].nveg + 1;
        if (options.AboveTreelineVeg >= 0) {
            veg_con_map[i].nv_active += 1;
        }

        veg_con[i] = malloc((veg_con_map[i].nv_active) * sizeof(*(veg_con[i])));
        check_alloc_status(veg_con[i], "Memory allocation error.");

        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].zone_depth = calloc(options.ROOT_ZONES,
                                              sizeof(*(veg_con[i][j].zone_depth)));
            check_alloc_status(veg_con[i][j].zone_depth,
                               "Memory allocation error.");
            veg_con[i][j].zone_fract = calloc(options.ROOT_ZONES,
                                              sizeof(*(veg_con[i][j].zone_fract)));
            check_alloc_status(veg_con[i][j].zone_fract,
                               "Memory allocation error.");
            if (options.CARBON) {
                veg_con[i][j].CanopLayerBnd = calloc(options.Ncanopy,
                                                     sizeof(*(veg_con[i][j].
                                                              CanopLayerBnd)));
                check_alloc_status(veg_con[i][j].CanopLayerBnd,
                                   "Memory allocation error.");
            }
            initialize_veg_con(&(veg_con[i][j]));
        }

        // vegetation library allocation - there is a veg library for each
        // active grid cell

        veg_lib[i] = calloc(options.NVEGTYPES, sizeof(*(veg_lib[i])));
        check_alloc_status(veg_lib[i], "Memory allocation error.");

        all_vars[i] = make_all_vars(veg_con_map[i].nv_active - 1);

        // allocate memory for veg_hist
        veg_hist[i] = calloc(veg_con_map[i].nv_active, sizeof(*(veg_hist[i])));
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            alloc_veg_hist(&(veg_hist[i][j]));
        }
    }
}
