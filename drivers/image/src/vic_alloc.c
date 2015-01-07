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
    size_t                     i;
    size_t                     j;

    // allocate memory for atmos structure
    atmos = (atmos_data_struct *)
            malloc((size_t) local_domain.ncells *
                   sizeof(atmos_data_struct));
    if (atmos == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for veg_hist structure
    veg_hist = (veg_hist_struct **)
               malloc((size_t) local_domain.ncells *
                      sizeof(veg_hist_struct *));
    if (veg_hist == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for soil structure
    soil_con = (soil_con_struct *)
               malloc((size_t) local_domain.ncells *
                      sizeof(soil_con_struct));
    if (soil_con == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation mapping structure
    veg_con_map = (veg_con_map_struct *)
                  malloc((size_t) local_domain.ncells *
                         sizeof(veg_con_map_struct));
    if (veg_con_map == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_con = (veg_con_struct **)
              malloc((size_t) local_domain.ncells *
                     sizeof(veg_con_struct *));
    if (veg_con == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_lib = (veg_lib_struct **)
              malloc((size_t) local_domain.ncells *
                     sizeof(veg_lib_struct *));
    if (veg_lib == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // all_vars allocation
    all_vars = (all_vars_struct *)
               malloc((size_t) local_domain.ncells *
                      sizeof(all_vars_struct));
    if (all_vars == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // out_data allocation
    out_data = (out_data_struct **)
               malloc((size_t) local_domain.ncells *
                      sizeof(out_data_struct *));
    if (out_data == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // save_data allocation
    save_data = (save_data_struct *)
                malloc((size_t) local_domain.ncells *
                       sizeof(save_data_struct));
    if (save_data == NULL) {
        log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for individual grid cells
    for (i = 0; i < local_domain.ncells; i++) {
        // atmos allocation - allocate enough memory for NR+1 steps
        alloc_atmos(&(atmos[i]));

        // snow band allocation
        soil_con[i].AreaFract = (double *) calloc(options.SNOW_BAND,
                                                  sizeof(double));
        if (soil_con[i].AreaFract == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].BandElev = (double *) calloc(options.SNOW_BAND,
                                                 sizeof(double));
        if (soil_con[i].BandElev == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Tfactor = (double *) calloc(options.SNOW_BAND,
                                                sizeof(double));
        if (soil_con[i].Tfactor == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Pfactor = (double *) calloc(options.SNOW_BAND,
                                                sizeof(double));
        if (soil_con[i].Pfactor == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        soil_con[i].AboveTreeLine = (bool *) calloc(options.SNOW_BAND,
                                                    sizeof(bool));
        if (soil_con[i].AboveTreeLine == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        initialize_soil_con(&(soil_con[i]));

        // vegetation tile allocation

        veg_con_map[i].nv_types = options.NVEGTYPES;

        veg_con_map[i].vidx = (int *) calloc(veg_con_map[i].nv_types,
                                             sizeof(int));
        if (veg_con_map[i].vidx == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }
        veg_con_map[i].Cv = (double *) calloc(veg_con_map[i].nv_types,
                                              sizeof(double));
        if (veg_con_map[i].Cv == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        veg_con_map[i].nv_active = (size_t) local_domain.locations[i].nveg + 1;
        if (options.AboveTreelineVeg >= 0) {
            veg_con_map[i].nv_active += 1;
        }

        veg_con[i] = (veg_con_struct *)
                     malloc((size_t) (veg_con_map[i].nv_active) *
                            sizeof(veg_con_struct));
        if (veg_con[i] == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].zone_depth = calloc(options.ROOT_ZONES,
                                              sizeof(float));
            if (veg_con[i][j].zone_depth == NULL) {
                log_err("Memory allocation error in vic_alloc().");
            }
            veg_con[i][j].zone_fract = calloc(options.ROOT_ZONES,
                                              sizeof(float));
            if (veg_con[i][j].zone_fract == NULL) {
                log_err("Memory allocation error in vic_alloc().");
            }
            if (options.CARBON) {
                veg_con[i][j].CanopLayerBnd = calloc(options.Ncanopy,
                                                     sizeof(double));
                if (veg_con[i][j].CanopLayerBnd == NULL) {
                    log_err("Memory allocation error in vic_alloc().");
                }
            }
            initialize_veg_con(&(veg_con[i][j]));
        }

        // vegetation library allocation - there is a veg library for each
        // active grid cell

        veg_lib[i] = (veg_lib_struct *) calloc(options.NVEGTYPES,
                                               sizeof(veg_lib_struct));
        if (veg_lib[i] == NULL) {
            log_err("Memory allocation error in vic_alloc().");
        }

        all_vars[i] = make_all_vars(veg_con_map[i].nv_active);

        out_data[i] = create_output_list();

        // allocate memory for veg_hist
        veg_hist[i] = (veg_hist_struct *) calloc(veg_con_map[i].nv_active,
                                                 sizeof(veg_hist_struct));
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            alloc_veg_hist(&(veg_hist[i][j]));
        }
    }
}
