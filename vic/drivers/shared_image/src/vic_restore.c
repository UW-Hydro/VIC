/******************************************************************************
 * @section DESCRIPTION
 *
 * Read initial model state.
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
 * @brief    Read initial model state.
 *****************************************************************************/
void
vic_restore(void)
{
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern veg_con_map_struct *veg_con_map;
    extern filenames_struct    filenames;

    int                        v;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     p;
    char                      *cvar = NULL;
    int                       *ivar = NULL;
    double                    *dvar = NULL;
    float                     *fvar = NULL;
    size_t                     d2count[2];
    size_t                     d2start[2];
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];
    size_t                     d5count[5];
    size_t                     d5start[5];
    size_t                     d6count[6];
    size_t                     d6start[6];

    // validate state file dimensions and coordinate variables
    check_init_state_file();

    // read state variables

    // allocate memory for variables to be stored
    cvar = malloc(local_domain.ncells_active * sizeof(*cvar));
    if (cvar == NULL) {
        log_err("Memory allocation error in vic_restore().");
    }

    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    if (ivar == NULL) {
        log_err("Memory allocation error in vic_restore().");
    }

    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_restore().");
    }

    fvar = malloc(local_domain.ncells_active * sizeof(*fvar));
    if (fvar == NULL) {
        log_err("Memory allocation error in vic_restore().");
    }

    // initialize starts and counts
    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain.n_ny;
    d2count[1] = global_domain.n_nx;

    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;
    d4count[0] = 1;
    d4count[1] = 1;
    d4count[2] = global_domain.n_ny;
    d4count[3] = global_domain.n_nx;

    d5start[0] = 0;
    d5start[1] = 0;
    d5start[2] = 0;
    d5start[3] = 0;
    d5start[4] = 0;
    d5count[0] = 1;
    d5count[1] = 1;
    d5count[2] = 1;
    d5count[3] = global_domain.n_ny;
    d5count[4] = global_domain.n_nx;

    d6start[0] = 0;
    d6start[1] = 0;
    d6start[2] = 0;
    d6start[3] = 0;
    d6start[4] = 0;
    d6start[5] = 0;
    d6count[0] = 1;
    d6count[1] = 1;
    d6count[2] = 1;
    d6count[3] = 1;
    d6count[4] = global_domain.n_ny;
    d6count[5] = global_domain.n_nx;

    // total soil moisture
    for (m = 0; m < options.NVEGTYPES; m++) {
        d5start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d5start[2] = j;
                get_scatter_nc_field_double(filenames.init_state,
                                            "Soil_moisture",
                                            d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].cell[v][k].layer[j].moist = dvar[i];
                    }
                }
            }
        }
    }

    // ice content
    for (m = 0; m < options.NVEGTYPES; m++) {
        d6start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d6start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d6start[2] = j;
                for (p = 0; p < options.Nfrost; p++) {
                    d6start[3] = p;
                    get_scatter_nc_field_double(filenames.init_state,
                                                "Soil_ice",
                                                d6start, d6count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        v = veg_con_map[i].vidx[m];
                        if (v >= 0) {
                            all_vars[i].cell[v][k].layer[j].ice[p] = dvar[i];
                        }
                    }
                }
            }
        }
    }

    // dew storage: tmpval = veg_var[veg][band].Wdew;
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Canopy_water",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].veg_var[v][k].Wdew = dvar[i];
                }
            }
        }
    }

    if (options.CARBON) {
        // cumulative NPP: tmpval = veg_var[veg][band].AnnualNPP;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.init_state,
                                            "AnnualNPP",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].veg_var[v][k].AnnualNPP = dvar[i];
                    }
                }
            }
        }

        // previous NPP: tmpval = veg_var[veg][band].AnnualNPPPrev;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.init_state,
                                            "AnnualNPPPrev",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].veg_var[v][k].AnnualNPPPrev = dvar[i];
                    }
                }
            }
        }

        // litter carbon: tmpval = cell[veg][band].CLitter;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.init_state,
                                            "CLitter",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].cell[v][k].CLitter = dvar[i];
                    }
                }
            }
        }

        // intermediate carbon: tmpval = tmpval = cell[veg][band].CInter;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.init_state,
                                            "CInter",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].cell[v][k].CInter = dvar[i];
                    }
                }
            }
        }

        // slow carbon: tmpval = cell[veg][band].CSlow;
        for (m = 0; m < options.NVEGTYPES; m++) {
            d4start[0] = m;
            for (k = 0; k < options.SNOW_BAND; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.init_state,
                                            "CSlow",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].cell[v][k].CSlow = dvar[i];
                    }
                }
            }
        }
    }

    // snow age: snow[veg][band].last_snow
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_int(filenames.init_state,
                                     "Snow_age",
                                     d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].last_snow = ivar[i];
                }
            }
        }
    }

    // melting state: (int)snow[veg][band].MELTING
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_int(filenames.init_state,
                                     "Snow_melt_state",
                                     d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].MELTING = ivar[i];
                }
            }
        }
    }

    // snow covered fraction: snow[veg][band].coverage
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_coverage",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].coverage = dvar[i];
                }
            }
        }
    }

    // snow water equivalent: snow[veg][band].swq
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_water_equivalent",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].swq = dvar[i];
                }
            }
        }
    }

    // snow surface temperature: snow[veg][band].surf_temp
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_surf_temp",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].surf_temp = dvar[i];
                }
            }
        }
    }

    // snow surface water: snow[veg][band].surf_water
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_surf_water",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].surf_water = dvar[i];
                }
            }
        }
    }

    // snow pack temperature: snow[veg][band].pack_temp
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_pack_temp",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].pack_temp = dvar[i];
                }
            }
        }
    }

    // snow pack water: snow[veg][band].pack_water
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_pack_water",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].pack_water = dvar[i];
                }
            }
        }
    }

    // snow density: snow[veg][band].density
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_density",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].density = dvar[i];
                }
            }
        }
    }

    // snow cold content: snow[veg][band].coldcontent
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_cold_content",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].coldcontent = dvar[i];
                }
            }
        }
    }

    // snow canopy storage: snow[veg][band].snow_canopy
    for (m = 0; m < options.NVEGTYPES; m++) {
        d4start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Snow_canopy",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[m];
                if (v >= 0) {
                    all_vars[i].snow[v][k].snow_canopy = dvar[i];
                }
            }
        }
    }

    // soil node temperatures: energy[veg][band].T[nidx]
    for (m = 0; m < options.NVEGTYPES; m++) {
        d5start[0] = m;
        for (k = 0; k < options.SNOW_BAND; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nnode; j++) {
                d5start[2] = j;
                get_scatter_nc_field_double(filenames.init_state,
                                            "Soil_node_temp",
                                            d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[m];
                    if (v >= 0) {
                        all_vars[i].energy[v][k].T[j] = dvar[i];
                    }
                }
            }
        }
    }

    if (options.LAKES) {
        // total soil moisture
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_soil_moisture",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.soil.layer[j].moist = dvar[i];
            }
        }

        // ice content
        for (j = 0; j < options.Nlayer; j++) {
            d4start[0] = j;
            for (p = 0; p < options.Nfrost; p++) {
                d4start[1] = p;
                get_scatter_nc_field_double(filenames.init_state,
                                            "Lake_soil_ice",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    all_vars[i].lake_var.soil.layer[j].ice[p] = dvar[i];
                }
            }
        }

        if (options.CARBON) {
            // litter carbon: tmpval = lake_var.soil.CLitter;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_CLitter",
                                        d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.soil.CLitter = dvar[i];
            }

            // intermediate carbon: tmpval = lake_var.soil.CInter;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_CInter",
                                        d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.soil.CInter = dvar[i];
            }

            // slow carbon: tmpval = lake_var.soil.CSlow;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_CSlow",
                                        d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.soil.CSlow = dvar[i];
            }
        }

        // snow age: lake_var.snow.last_snow
        get_scatter_nc_field_int(filenames.init_state,
                                 "Lake_snow_age",
                                 d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.last_snow = ivar[i];
        }

        // melting state: (int)lake_var.snow.MELTING
        get_scatter_nc_field_int(filenames.init_state,
                                 "Lake_snow_melt_state",
                                 d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.MELTING = ivar[i];
        }

        // snow covered fraction: lake_var.snow.coverage
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_coverage",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.coverage = dvar[i];
        }

        // snow water equivalent: lake_var.snow.swq
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_water_equivalent",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.swq = dvar[i];
        }

        // snow surface temperature: lake_var.snow.surf_temp
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_surf_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.surf_temp = dvar[i];
        }

        // snow surface water: lake_var.snow.surf_water
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_surf_water",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.surf_water = dvar[i];
        }

        // snow pack temperature: lake_var.snow.pack_temp
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_pack_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.pack_temp = dvar[i];
        }

        // snow pack water: lake_var.snow.pack_water
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_pack_water",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.pack_water = dvar[i];
        }

        // snow density: lake_var.snow.density
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_density",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.density = dvar[i];
        }

        // snow cold content: lake_var.snow.coldcontent
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_cold_content",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.coldcontent = dvar[i];
        }

        // snow canopy storage: lake_var.snow.snow_canopy
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_snow_canopy",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.snow.snow_canopy = dvar[i];
        }

        // soil node temperatures: lake_var.energy.T[nidx]
        for (j = 0; j < options.Nnode; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_soil_node_temp",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.soil.layer[j].moist = dvar[i];
            }
        }

        // lake active layers: lake_var.activenod
        get_scatter_nc_field_int(filenames.init_state,
                                 "Lake_active_layers",
                                 d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.activenod = ivar[i];
        }

        // lake layer thickness: lake_var.dz
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_layer_dz",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.dz = dvar[i];
        }

        // lake surface layer thickness: lake_var.surfdz
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_surf_layer_dz",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.surfdz = dvar[i];
        }

        // lake depth: lake_var.ldepth
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_depth",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.ldepth = dvar[i];
        }

        // lake layer surface areas: lake_var.surface[ndix]
        for (j = 0; j < options.NLAKENODES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_layer_surf_area",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.surface[j] = dvar[i];
            }
        }

        // lake surface area: lake_var.sarea
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_surf_area",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.sarea = dvar[i];
        }

        // lake volume: lake_var.volume
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_volume",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.volume = dvar[i];
        }

        // lake layer temperatures: lake_var.temp[nidx]
        for (j = 0; j < options.NLAKENODES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.init_state,
                                        "Lake_layer_temp",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                all_vars[i].lake_var.temp[j] = dvar[i];
            }
        }

        // vertical average lake temperature: lake_var.tempavg
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_average_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.tempavg = dvar[i];
        }

        // lake ice area fraction: lake_var.areai
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_area_frac",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.areai = dvar[i];
        }

        // new lake ice area fraction: lake_var.new_ice_area
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_area_frac_new",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.new_ice_area = dvar[i];
        }

        // lake ice water equivalent: lake_var.ice_water_eq
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_water_equivalent",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.ice_water_eq = dvar[i];
        }

        // lake ice height: lake_var.hice
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_height",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.hice = dvar[i];
        }

        // lake ice temperature: lake_var.tempi
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.tempi = dvar[i];
        }

        // lake ice snow water equivalent: lake_var.swe
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_water_equivalen",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.swe = dvar[i];
        }

        // lake ice snow surface temperature: lake_var.surf_temp
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_surf_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.surf_temp = dvar[i];
        }

        // lake ice snow pack temperature: lake_var.pack_temp
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_pack_temp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.pack_temp = dvar[i];
        }

        // lake ice snow coldcontent: lake_var.coldcontent
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_cold_content",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.coldcontent = dvar[i];
        }

        // lake ice snow surface water: lake_var.surf_water
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_surf_water",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.surf_water = dvar[i];
        }

        // lake ice snow pack water: lake_var.pack_water
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_pack_water",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.pack_water = dvar[i];
        }

        // lake ice snow albedo: lake_var.SAlbedo
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_albedo",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.SAlbedo = dvar[i];
        }

        // lake ice snow depth: lake_var.sdepth
        get_scatter_nc_field_double(filenames.init_state,
                                    "Lake_ice_snow_depth",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            all_vars[i].lake_var.sdepth = dvar[i];
        }
    }

    free(cvar);
    free(ivar);
    free(dvar);
    free(fvar);
}

void
check_init_state_file(void)
{
    extern filenames_struct filenames;
    extern domain_struct    global_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;

    int                     status;
    size_t                  dimlen;
    size_t                  i;
    size_t                  d1count[1];
    size_t                  d1start[1];
    size_t                  d2count[2];
    size_t                  d2start[2];
    int                     lon_var_id;
    int                     lat_var_id;
    double                 *dvar;
    double                  rtol = 0.0; // maybe move this to a .h file
    double                  abs_tol = 0.0001; // maybe move this to a .h file
    nc_file_struct          nc;

    // open the netcdf file
    status = nc_open(filenames.init_state, NC_SHARE, &(nc.nc_id));
    if (status != NC_NOERR) {
        log_err("Error opening %s", filenames.init_state);
    }

    // read and validate dimension lengths
    dimlen = get_nc_dimension(filenames.init_state, global_domain.info.x_dim);
    if (dimlen != global_domain.n_nx) {
        log_err("Number of grid columns in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, global_domain.info.y_dim);
    if (dimlen != global_domain.n_ny) {
        log_err("Number of grid rows in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, "veg_class");
    if (dimlen != options.NVEGTYPES) {
        log_err("Number of veg classes in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, "snow_band");
    if (dimlen != options.SNOW_BAND) {
        log_err("Number of snow bands in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, "nlayer");
    if (dimlen != options.Nlayer) {
        log_err("Number of soil layers in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, "frost_area");
    if (dimlen != options.Nfrost) {
        log_err("Number of frost areas in state file does not "
                "match parameter file");
    }
    dimlen = get_nc_dimension(filenames.init_state, "soil_node");
    if (dimlen != options.Nnode) {
        log_err("Number of soil nodes in state file does not "
                "match parameter file");
    }
    if (options.LAKES) {
        dimlen = get_nc_dimension(filenames.init_state, "lake_node");
        if (dimlen != options.NLAKENODES) {
            log_err("Number of lake nodes in state file does not "
                    "match parameter file");
        }
    }

    // read dimension variables

    // lat/lon
    status = nc_inq_varid(nc.nc_id, global_domain.info.lon_var, &lon_var_id);
    if (status != NC_NOERR) {
        log_err("Unable to find variable \"%s\" in %s",
                global_domain.info.lon_var, filenames.init_state);
    }
    status = nc_inq_varid(nc.nc_id, global_domain.info.lat_var, &lat_var_id);
    if (status != NC_NOERR) {
        log_err("Unable to find variable \"%s\" in %s",
                global_domain.info.lat_var, filenames.init_state);
    }
    if (global_domain.info.n_coord_dims == 1) {
        d1start[0] = 0;
        dvar = calloc(global_domain.n_nx, sizeof(*dvar));
        d1count[0] = global_domain.n_nx;
        status = nc_get_vara_double(nc.nc_id, lon_var_id,
                                    d1start, d1count, dvar);
        if (status != NC_NOERR) {
            log_err("Error reading data from \"%s\" in %s",
                    global_domain.info.lon_var, filenames.init_state);
        }
        for (i = 0; i < global_domain.n_nx; i++) {
            if (!assert_close_double(dvar[i], global_domain.locations[i].longitude, rtol, abs_tol)) {
                log_err("Longitudes in initial state file do not "
                        "match parameter file");
            }
        }
        free(dvar);

        dvar = calloc(global_domain.n_ny, sizeof(*dvar));
        d1count[0] = global_domain.n_ny;
        status = nc_get_vara_double(nc.nc_id, lat_var_id,
                                    d1start, d1count, dvar);
        if (status != NC_NOERR) {
            log_err("Error reading data from \"%s\" in %s",
                    global_domain.info.lat_var, filenames.init_state);
        }
        for (i = 0; i < global_domain.n_ny; i++) {
            if (!assert_close_double(dvar[i], global_domain.locations[i].latitude, rtol, abs_tol)) {
                log_err("Latitudes in initial state file do not "
                        "match parameter file");
            }
        }
        free(dvar);
    }
    else if (global_domain.info.n_coord_dims == 2) {
        d2start[0] = 0;
        d2start[1] = 0;
        dvar = calloc(global_domain.n_ny * global_domain.n_nx, sizeof(*dvar));
        d2count[0] = global_domain.n_ny;
        d2count[1] = global_domain.n_nx;
        status = nc_get_vara_double(nc.nc_id, lon_var_id,
                                    d2start, d2count, dvar);
        if (status != NC_NOERR) {
            log_err("Error reading data from \"%s\" in %s",
                    global_domain.info.lon_var, filenames.init_state);
        }
        for (i = 0; i < global_domain.n_ny * global_domain.n_nx; i++) {
            if (dvar[i] != (double) global_domain.locations[i].longitude) {
                log_err("Longitudes in initial state file do not "
                        "match parameter file");
            }
        }
        status = nc_get_vara_double(nc.nc_id, lat_var_id,
                                    d2start, d2count, dvar);
        if (status != NC_NOERR) {
            log_err("Error reading data from \"%s\" in %s",
                    global_domain.info.lat_var, filenames.init_state);
        }
        for (i = 0; i < global_domain.n_ny * global_domain.n_nx; i++) {
            if (dvar[i] != (double) global_domain.locations[i].latitude) {
                log_err("Latitudes in initial state file do not "
                        "match parameter file");
            }
        }
        free(dvar);
    }
    else {
        log_err("global_domain.info.n_coord_dims should be 1 or 2");
    }

    // Variables for other dimensions
    d1start[0] = 0;
    d1count[0] = options.Nnode;

    // soil thermal node deltas
    dvar = calloc(options.Nnode, sizeof(*dvar));
    get_nc_field_double(filenames.init_state, "dz_node",
                        d1start, d1count, dvar);
    for (i = 0; i < options.Nnode; i++) {
        if (dvar[i] != soil_con[0].dz_node[i]) {
            log_err("Soil node intervals in state file do not match "
                    "those computed by VIC");
        }
    }
    free(dvar);

    // soil thermal node depths
    dvar = calloc(options.Nnode, sizeof(*dvar));
    get_nc_field_double(filenames.init_state, "node_depth",
                        d1start, d1count, dvar);
    for (i = 0; i < options.Nnode; i++) {
        if (dvar[i] != soil_con[0].Zsum_node[i]) {
            log_err("Soil node depths in state file do not match "
                    "those computed by VIC");
        }
    }
    free(dvar);
}
