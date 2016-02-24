/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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
 * @brief    Save model state.
 *****************************************************************************/
void
vic_store(void)
{
    extern all_vars_struct    *all_vars;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;

    int                        status;
    int                        dimids[MAXDIMS];
    int                        v;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     p;
    size_t                     ndims;
    char                      *cvar = NULL;
    int                       *ivar = NULL;
    double                    *dvar = NULL;
    float                     *fvar = NULL;
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];
    size_t                     d5count[5];
    size_t                     d5start[5];
    size_t                     d6count[6];
    size_t                     d6start[6];

    nc_file_struct             nc_state_file;

    // allocate memory for variables to be stored
    cvar = malloc(local_domain.ncells_active * sizeof(*cvar));
    if (cvar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    if (ivar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    fvar = malloc(local_domain.ncells_active * sizeof(*fvar));
    if (fvar == NULL) {
        log_err("Memory allocation error in vic_store().");
    }

    // initialize starts and counts
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

    // create netcdf file for storing model state - keep the file open
    initialize_state_file(&nc_state_file);

    // set missing values
    for (i = 0; i < local_domain.ncells_active; i++) {
        cvar[i] = nc_state_file.c_fillvalue;
        ivar[i] = nc_state_file.i_fillvalue;
        dvar[i] = nc_state_file.d_fillvalue;
        fvar[i] = nc_state_file.f_fillvalue;
    }

    // initialize dimids to invalid values
    for (i = 0; i < MAXDIMS; i++) {
        dimids[i] = -1;
    }

    // soil thermal node deltas
    ndims = 3;
    dimids[0] = nc_state_file.node_dimid;
    dimids[1] = nc_state_file.nj_dimid;
    dimids[2] = nc_state_file.ni_dimid;
    for (j = 0; j < options.Nnode; j++) {
        d3start[0] = j;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) soil_con[i].dz_node[j];
        }
        gather_put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "dz_node", d3start, d3count,
                                   dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }


    // soil thermal node depths
    ndims = 3;
    dimids[0] = nc_state_file.node_dimid;
    dimids[1] = nc_state_file.nj_dimid;
    dimids[2] = nc_state_file.ni_dimid;
    for (j = 0; j < options.Nnode; j++) {
        d3start[0] = j;
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = (double) soil_con[i].Zsum_node[j];
        }
        gather_put_nc_field_double(nc_state_file.fname, &(nc_state_file.open),
                                   &(nc_state_file.nc_id),
                                   nc_state_file.d_fillvalue,
                                   dimids, ndims, "Zsum_node", d3start, d3count,
                                   dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            dvar[i] = nc_state_file.d_fillvalue;
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // total soil moisture
    ndims = 5;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.layer_dimid;
    dimids[3] = nc_state_file.nj_dimid;
    dimids[4] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d5start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d5start[2] = j;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][m].layer[j].moist;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Soil_moisture",
                                           d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // total soil moisture
    ndims = 6;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.layer_dimid;
    dimids[3] = nc_state_file.frost_dimid;
    dimids[4] = nc_state_file.nj_dimid;
    dimids[5] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d6start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d6start[1] = k;
            for (j = 0; j < options.Nlayer; j++) {
                d6start[2] = j;
                for (p = 0; p < options.Nfrost; p++) {
                    d6start[3] = p;
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        v = veg_con_map[i].vidx[k];
                        if (v >= 0) {
                            dvar[i] = (double)
                                      all_vars[i].cell[v][m].layer[j].ice[p];
                        }
                        else {
                            dvar[i] = nc_state_file.d_fillvalue;
                        }
                    }
                    gather_put_nc_field_double(nc_state_file.fname,
                                               &(nc_state_file.open),
                                               &(nc_state_file.nc_id),
                                               nc_state_file.d_fillvalue,
                                               dimids, ndims, "Ice_content",
                                               d6start, d6count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // dew storage: tmpval = veg_var[veg][band].Wdew;
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].veg_var[v][m].Wdew;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Wdew",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    if (options.CARBON) {
        // cumulative NPP: tmpval = veg_var[veg][band].AnnualNPP;
        ndims = 4;
        dimids[0] = nc_state_file.band_dimid;
        dimids[1] = nc_state_file.veg_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.SNOW_BAND; m++) {
            d4start[0] = m;
            for (k = 0; k < options.NVEGTYPES; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].veg_var[v][m].AnnualNPP;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "AnnualNPP",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // previous NPP: tmpval = veg_var[veg][band].AnnualNPPPrev;
        ndims = 4;
        dimids[0] = nc_state_file.band_dimid;
        dimids[1] = nc_state_file.veg_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.SNOW_BAND; m++) {
            d4start[0] = m;
            for (k = 0; k < options.NVEGTYPES; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].veg_var[v][m].AnnualNPPPrev;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "AnnualNPPPrev",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // litter carbon: tmpval = cell[veg][band].CLitter;
        ndims = 4;
        dimids[0] = nc_state_file.band_dimid;
        dimids[1] = nc_state_file.veg_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.SNOW_BAND; m++) {
            d4start[0] = m;
            for (k = 0; k < options.NVEGTYPES; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][m].CLitter;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "CLitter",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // intermediate carbon: tmpval = tmpval = cell[veg][band].CInter;
        ndims = 4;
        dimids[0] = nc_state_file.band_dimid;
        dimids[1] = nc_state_file.veg_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.SNOW_BAND; m++) {
            d4start[0] = m;
            for (k = 0; k < options.NVEGTYPES; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][m].CInter;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Cinter",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }

        // slow carbon: tmpval = cell[veg][band].CSlow;
        ndims = 4;
        dimids[0] = nc_state_file.band_dimid;
        dimids[1] = nc_state_file.veg_dimid;
        dimids[2] = nc_state_file.nj_dimid;
        dimids[3] = nc_state_file.ni_dimid;
        for (m = 0; m < options.SNOW_BAND; m++) {
            d4start[0] = m;
            for (k = 0; k < options.NVEGTYPES; k++) {
                d4start[1] = k;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].cell[v][m].CSlow;
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "CSlow",
                                           d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
        for (i = 0; i < ndims; i++) {
            dimids[i] = -1;
        }
    }

    // last day of snow: snow[veg][band].last_snow
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    ivar[i] = (int)
                              all_vars[i].snow[v][m].last_snow;
                }
                else {
                    ivar[i] = nc_state_file.i_fillvalue;
                }
            }
            gather_put_nc_field_int(nc_state_file.fname,
                                    &(nc_state_file.open),
                                    &(nc_state_file.nc_id),
                                    nc_state_file.i_fillvalue,
                                    dimids, ndims, "Last_snow",
                                    d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                ivar[i] = nc_state_file.i_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // melting state: (int)snow[veg][band].MELTING
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    ivar[i] = (int)
                              all_vars[i].snow[v][m].MELTING;
                }
                else {
                    ivar[i] = nc_state_file.i_fillvalue;
                }
            }
            gather_put_nc_field_int(nc_state_file.fname,
                                    &(nc_state_file.open),
                                    &(nc_state_file.nc_id),
                                    nc_state_file.i_fillvalue,
                                    dimids, ndims, "Melt_state",
                                    d4start, d4count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                ivar[i] = nc_state_file.i_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow covered fraction: snow[veg][band].coverage
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].coverage;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_coverage",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow water equivalent: snow[veg][band].swq
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].swq;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_water_equivalent",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow surface temperature: snow[veg][band].surf_temp
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].surf_temp;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_surf_temp",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow surface water: snow[veg][band].surf_water
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].surf_water;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_surf_water",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow pack temperature: snow[veg][band].pack_temp
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].pack_temp;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_pack_temp",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow pack water: snow[veg][band].pack_water
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].pack_water;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_pack_water",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow density: snow[veg][band].density
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].density;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_density",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow cold content: snow[veg][band].coldcontent
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].coldcontent;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_cold_content",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // snow canopy storage: snow[veg][band].snow_canopy
    ndims = 4;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.nj_dimid;
    dimids[3] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d4start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d4start[1] = k;
            for (i = 0; i < local_domain.ncells_active; i++) {
                v = veg_con_map[i].vidx[k];
                if (v >= 0) {
                    dvar[i] = (double)
                              all_vars[i].snow[v][m].snow_canopy;
                }
                else {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
            gather_put_nc_field_double(nc_state_file.fname,
                                       &(nc_state_file.open),
                                       &(nc_state_file.nc_id),
                                       nc_state_file.d_fillvalue,
                                       dimids, ndims, "Snow_canopy",
                                       d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                dvar[i] = nc_state_file.d_fillvalue;
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    // soil node temperatures: energy[veg][band].T[nidx]
    ndims = 5;
    dimids[0] = nc_state_file.band_dimid;
    dimids[1] = nc_state_file.veg_dimid;
    dimids[2] = nc_state_file.node_dimid;
    dimids[3] = nc_state_file.nj_dimid;
    dimids[4] = nc_state_file.ni_dimid;
    for (m = 0; m < options.SNOW_BAND; m++) {
        d5start[0] = m;
        for (k = 0; k < options.NVEGTYPES; k++) {
            d5start[1] = k;
            for (j = 0; j < options.Nnode; j++) {
                d5start[2] = j;
                for (i = 0; i < local_domain.ncells_active; i++) {
                    v = veg_con_map[i].vidx[k];
                    if (v >= 0) {
                        dvar[i] = (double)
                                  all_vars[i].energy[v][m].T[j];
                    }
                    else {
                        dvar[i] = nc_state_file.d_fillvalue;
                    }
                }
                gather_put_nc_field_double(nc_state_file.fname,
                                           &(nc_state_file.open),
                                           &(nc_state_file.nc_id),
                                           nc_state_file.d_fillvalue,
                                           dimids, ndims, "Node_temperature",
                                           d5start, d5count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    dvar[i] = nc_state_file.d_fillvalue;
                }
            }
        }
    }
    for (i = 0; i < ndims; i++) {
        dimids[i] = -1;
    }

    if (options.LAKES) {
        log_err("LAKES option not yet implemented in vic_store()");
    }

    // close the netcdf file if it is still open
    if (nc_state_file.open == true) {
        status = nc_close(nc_state_file.nc_id);
        if (status != NC_NOERR) {
            log_err("Error closing %s", nc_state_file.fname);
        }
    }

    free(cvar);
    free(ivar);
    free(dvar);
    free(fvar);
}

void
initialize_state_file(nc_file_struct *nc)
{
    extern size_t           current;
    extern dmy_struct      *dmy;
    extern filenames_struct filenames;
    extern domain_struct    global_domain;
    extern option_struct    options;

    int                     status;
    int                     old_fill_mode;

    sprintf(nc->fname, "%s.%04d%02d%02d_%05u.nc",
            filenames.statefile, dmy[current].year, dmy[current].month,
            dmy[current].day, dmy[current].dayseconds);

    nc->c_fillvalue = NC_FILL_CHAR;
    nc->i_fillvalue = NC_FILL_INT;
    nc->d_fillvalue = NC_FILL_DOUBLE;
    nc->f_fillvalue = NC_FILL_FLOAT;

    nc->band_size = options.SNOW_BAND;
    nc->frost_size = options.Nfrost;
    nc->layer_size = options.Nlayer;
    nc->ni_size = global_domain.n_nx;
    nc->nj_size = global_domain.n_ny;
    nc->node_size = options.Nnode;
    nc->root_zone_size = options.ROOT_ZONES;
    nc->veg_size = options.NVEGTYPES;

    // open the netcdf file
    status = nc_create(nc->fname, NC_NETCDF4 | NC_CLASSIC_MODEL, &(nc->nc_id));
    if (status != NC_NOERR) {
        log_err("Error creating %s", nc->fname);
    }
    nc->open = true;

    // set the NC_FILL attribute
    status = nc_set_fill(nc->nc_id, NC_FILL, &old_fill_mode);
    if (status != NC_NOERR) {
        log_err("Error setting fill value in %s", nc->fname);
    }

    // define netcdf dimensions
    status = nc_def_dim(nc->nc_id, "snow_band", nc->band_size,
                        &(nc->band_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining snow_band in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "frost_area", nc->frost_size,
                        &(nc->frost_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining frost_area in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "nlayer", nc->layer_size,
                        &(nc->layer_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining nlayer in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "ni", nc->ni_size, &(nc->ni_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining ni in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "nj", nc->nj_size, &(nc->nj_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining nj in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "node", nc->node_size, &(nc->node_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining node in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "root_zone", nc->root_zone_size,
                        &(nc->root_zone_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining root_zone in %s", nc->fname);
    }

    status = nc_def_dim(nc->nc_id, "veg_class", nc->veg_size,
                        &(nc->veg_dimid));
    if (status != NC_NOERR) {
        log_err("Error defining veg_class in %s", nc->fname);
    }

    // leave define mode
    status = nc_enddef(nc->nc_id);
    if (status != NC_NOERR) {
        log_err("Error leaving define mode for %s", nc->fname);
    }
}
