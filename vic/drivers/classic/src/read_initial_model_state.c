/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initializes the model state at hour 0 of the date defined in
 * the given state file.
 *
 * Soil moisture, soil thermal, and snowpack variables  are stored for each
 * vegetation type and snow band.  However moisture variables from the
 * distributed precipitation model are averaged so that the model is restarted
 * with mu = 1.
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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This subroutine initializes the model state at hour 0 of the date
 *           defined in the given state file.
 *****************************************************************************/
void
read_initial_model_state(FILE            *init_state,
                         all_vars_struct *all_vars,
                         int              Nveg,
                         int              Nbands,
                         int              cellnum,
                         soil_con_struct *soil_con,
                         lake_con_struct  lake_con)
{
    extern option_struct options;

    char                 tmpstr[MAXSTRING];
    char                 tmpchar;
    int                  veg, iveg;
    int                  band, iband;
    size_t               lidx;
    size_t               nidx;
    int                  tmp_cellnum;
    int                  tmp_Nveg;
    int                  tmp_Nband;
    int                  tmp_char;
    int                  byte, Nbytes;
    int                  node;
    size_t               frost_area;

    cell_data_struct   **cell;
    snow_data_struct   **snow;
    energy_bal_struct  **energy;
    veg_var_struct     **veg_var;
    lake_var_struct     *lake_var;

    cell = all_vars->cell;
    veg_var = all_vars->veg_var;
    snow = all_vars->snow;
    energy = all_vars->energy;
    lake_var = &all_vars->lake_var;

    /* read cell information */
    if (options.STATE_FORMAT == BINARY) {
        fread(&tmp_cellnum, sizeof(int), 1, init_state);
        fread(&tmp_Nveg, sizeof(int), 1, init_state);
        fread(&tmp_Nband, sizeof(int), 1, init_state);
        fread(&Nbytes, sizeof(int), 1, init_state);
    }
    else {
        fscanf(init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
    }
    // Skip over unused cell information
    while (tmp_cellnum != cellnum && !feof(init_state)) {
        if (options.STATE_FORMAT == BINARY) {
            // skip rest of current cells info
            for (byte = 0; byte < Nbytes; byte++) {
                fread(&tmpchar, 1, 1, init_state);
            }
            // read info for next cell
            fread(&tmp_cellnum, sizeof(int), 1, init_state);
            fread(&tmp_Nveg, sizeof(int), 1, init_state);
            fread(&tmp_Nband, sizeof(int), 1, init_state);
            fread(&Nbytes, sizeof(int), 1, init_state);
        }
        else {
            // skip rest of current cells info
            fgets(tmpstr, MAXSTRING, init_state); // skip rest of general cell info
            for (veg = 0; veg <= tmp_Nveg; veg++) {
                for (band = 0; band < tmp_Nband; band++) {
                    fgets(tmpstr, MAXSTRING, init_state); // skip snowband info
                }
            }
            if (options.LAKES) {
                fgets(tmpstr, MAXSTRING, init_state); // skip lake info
            }
            // read info for next cell
            fscanf(init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
        } // end if
    } // end while

    if (feof(init_state)) {
        log_err("Requested grid cell (%d) is not in the model state file.",
                cellnum);
    }

    if (tmp_Nveg != Nveg) {
        log_err("The number of vegetation types in cell %d (%d) does not equal"
                "that defined in vegetation parameter file (%d).  Check your "
                "input files.", cellnum, tmp_Nveg, Nveg);
    }
    if (tmp_Nband != Nbands) {
        log_err("The number of snow bands in cell %d (%d) does not equal that"
                "defined in the snow band file (%d).  Check your input files.",
                cellnum, tmp_Nband, Nbands);
    }

    /* Read soil thermal node deltas */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fread(&soil_con->dz_node[nidx], sizeof(double), 1, init_state);
        }
        else {
            fscanf(init_state, "%lf", &soil_con->dz_node[nidx]);
        }
    }
    if (options.Nnode == 1) {
        soil_con->dz_node[0] = 0;
    }

    /* Read soil thermal node depths */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fread(&soil_con->Zsum_node[nidx], sizeof(double), 1, init_state);
        }
        else {
            fscanf(init_state, "%lf", &soil_con->Zsum_node[nidx]);
        }
    }
    if (options.Nnode == 1) {
        soil_con->Zsum_node[0] = 0;
    }
    if (soil_con->Zsum_node[options.Nnode - 1] - soil_con->dp > DBL_EPSILON) {
        log_warn("Sum of soil nodes (%f) exceeds defined damping depth (%f)."
                 "Resetting damping depth.",
                 soil_con->Zsum_node[options.Nnode - 1], soil_con->dp);
        soil_con->dp = soil_con->Zsum_node[options.Nnode - 1];
    }

    /* Input for all vegetation types */
    for (veg = 0; veg <= Nveg; veg++) {
        /* Input for all snow bands */
        for (band = 0; band < Nbands; band++) {
            /* Read cell identification information */
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&iveg, sizeof(int), 1, init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&iband, sizeof(int), 1, init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, "%d %d", &iveg, &iband) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            if (iveg != veg || iband != band) {
                log_err("The vegetation and snow band indices in the model "
                        "state file (veg = %d, band = %d) do not match those "
                        "currently requested (veg = %d , band = %d).  Model "
                        "state file must be stored with variables for all "
                        "vegetation indexed by variables for all snow bands.",
                        iveg, iband, veg, band);
            }

            /* Read total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                if (options.STATE_FORMAT == BINARY) {
                    if (fread(&cell[veg][band].layer[lidx].moist,
                              sizeof(double), 1, init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
                else {
                    if (fscanf(init_state, " %lf",
                               &cell[veg][band].layer[lidx].moist) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
            }

            /* Read average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    if (options.STATE_FORMAT == BINARY) {
                        if (fread(&cell[veg][band].layer[lidx].ice[frost_area],
                                  sizeof(double), 1, init_state) != 1) {
                            log_err("End of model state file found"
                                    "unexpectedly");
                        }
                    }
                    else {
                        if (fscanf(init_state, " %lf",
                                   &cell[veg][band].layer[lidx].ice[frost_area])
                            ==
                            EOF) {
                            log_err("End of model state file found"
                                    "unexpectedly");
                        }
                    }
                }
            }

            if (veg < Nveg) {
                /* Read dew storage */
                if (options.STATE_FORMAT == BINARY) {
                    if (fread(&veg_var[veg][band].Wdew, sizeof(double), 1,
                              init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
                else {
                    if (fscanf(init_state, " %lf",
                               &veg_var[veg][band].Wdew) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                }

                if (options.CARBON) {
                    if (options.STATE_FORMAT == BINARY) {
                        /* Read cumulative annual NPP */
                        if (fread(&(veg_var[veg][band].AnnualNPP),
                                  sizeof(double), 1, init_state) != 1) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fread(&(veg_var[veg][band].AnnualNPPPrev),
                                  sizeof(double), 1, init_state) != 1) {
                            log_err("End of model state file found unexpectedly");
                        }
                        /* Read Soil Carbon Storage */
                        if (fread(&(cell[veg][band].CLitter), sizeof(double), 1,
                                  init_state) != 1) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fread(&(cell[veg][band].CInter), sizeof(double), 1,
                                  init_state) != 1) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fread(&(cell[veg][band].CSlow), sizeof(double), 1,
                                  init_state) != 1) {
                            log_err("End of model state file found unexpectedly");
                        }
                    }
                    else {
                        /* Read cumulative annual NPP */
                        if (fscanf(init_state, " %lf",
                                   &veg_var[veg][band].AnnualNPP) == EOF) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fscanf(init_state, " %lf",
                                   &veg_var[veg][band].AnnualNPPPrev) == EOF) {
                            log_err("End of model state file found unexpectedly");
                        }
                        /* Read Soil Carbon Storage */
                        if (fscanf(init_state, " %lf",
                                   &cell[veg][band].CLitter) == EOF) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fscanf(init_state, " %lf",
                                   &cell[veg][band].CInter) == EOF) {
                            log_err("End of model state file found unexpectedly");
                        }
                        if (fscanf(init_state, " %lf",
                                   &cell[veg][band].CSlow) == EOF) {
                            log_err("End of model state file found unexpectedly");
                        }
                    }
                }
            }

            /* Read snow data */
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg][band].last_snow, sizeof(int), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].MELTING, sizeof(char), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].coverage, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].swq, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].surf_temp, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].surf_water, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].pack_temp, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].pack_water, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].density, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].coldcontent, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&snow[veg][band].snow_canopy, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state,
                           " %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                           &snow[veg][band].last_snow, &tmp_char,
                           &snow[veg][band].coverage, &snow[veg][band].swq,
                           &snow[veg][band].surf_temp,
                           &snow[veg][band].surf_water,
                           &snow[veg][band].pack_temp,
                           &snow[veg][band].pack_water,
                           &snow[veg][band].density,
                           &snow[veg][band].coldcontent,
                           &snow[veg][band].snow_canopy) ==
                    EOF) {
                    log_err("End of model state file found unexpectedly");
                }
                snow[veg][band].MELTING = (char)tmp_char;
            }

            /* Read soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                if (options.STATE_FORMAT == BINARY) {
                    if (fread(&energy[veg][band].T[nidx], sizeof(double), 1,
                              init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
                else {
                    if (fscanf(init_state, " %lf",
                               &energy[veg][band].T[nidx]) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
            }

            /* Read foliage temperature*/
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&energy[veg][band].Tfoliage, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &energy[veg][band].Tfoliage) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read outgoing longwave from understory */
            /* TO-DO: this is a flux. Saving it to the state file is a temporary solution! */
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&energy[veg][band].LongUnderOut, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &energy[veg][band].LongUnderOut) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read thermal flux through the snow pack */
            /* TO-DO: this is a flux. Saving it to the state file is a temporary solution! */
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&energy[veg][band].snow_flux, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &energy[veg][band].snow_flux) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }
    }
    if (options.LAKES) {
        if (options.STATE_FORMAT == BINARY) {
            /* Read total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                if (fread(&lake_var->soil.layer[lidx].moist, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    if (fread(&lake_var->soil.layer[lidx].ice[frost_area],
                              sizeof(double), 1, init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
            }

            if (options.CARBON) {
                /* Read Soil Carbon Storage */
                if (fread(&(lake_var->soil.CLitter), sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&(lake_var->soil.CInter), sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fread(&(lake_var->soil.CSlow), sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read snow data */
            if (fread(&lake_var->snow.last_snow, sizeof(int), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.MELTING, sizeof(char), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.coverage, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.swq, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.surf_temp, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.surf_water, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.pack_temp, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.pack_water, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.density, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.coldcontent, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->snow.snow_canopy, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (lake_var->snow.density > 0.) {
                lake_var->snow.depth = CONST_RHOFW * lake_var->snow.swq /
                                       lake_var->snow.density;
            }

            /* Read soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                if (fread(&lake_var->energy.T[nidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read lake-specific variables */
            if (fread(&lake_var->activenod, sizeof(int), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->dz, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->surfdz, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->ldepth, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            for (node = 0; node <= lake_var->activenod; node++) {
                if (fread(&lake_var->surface[node], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            if (fread(&lake_var->sarea, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->volume, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            for (node = 0; node < lake_var->activenod; node++) {
                if (fread(&lake_var->temp[node], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            if (fread(&lake_var->tempavg, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->areai, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->new_ice_area, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->ice_water_eq, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->hice, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->tempi, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->swe, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->surf_temp, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->pack_temp, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->coldcontent, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->surf_water, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->pack_water, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->SAlbedo, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&lake_var->sdepth, sizeof(double), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
        }
        else {
            /* Read total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                if (fscanf(init_state, " %lf",
                           &lake_var->soil.layer[lidx].moist) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    if (fscanf(init_state, " %lf",
                               &lake_var->soil.layer[lidx].ice[frost_area]) ==
                        EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
            }

            if (options.CARBON) {
                /* Read Soil Carbon Storage */
                if (fscanf(init_state, " %lf",
                           &lake_var->soil.CLitter) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fscanf(init_state, " %lf", &lake_var->soil.CInter) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
                if (fscanf(init_state, " %lf", &lake_var->soil.CSlow) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read snow data */
            if (fscanf(init_state, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &lake_var->snow.last_snow, &tmp_char,
                       &lake_var->snow.coverage, &lake_var->snow.swq,
                       &lake_var->snow.surf_temp, &lake_var->snow.surf_water,
                       &lake_var->snow.pack_temp, &lake_var->snow.pack_water,
                       &lake_var->snow.density, &lake_var->snow.coldcontent,
                       &lake_var->snow.snow_canopy) ==
                EOF) {
                log_err("End of model state file found unexpectedly");
            }
            lake_var->snow.MELTING = (char)tmp_char;
            if (lake_var->snow.density > 0.) {
                lake_var->snow.depth = CONST_RHOFW * lake_var->snow.swq /
                                       lake_var->snow.density;
            }

            /* Read soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                if (fscanf(init_state, " %lf",
                           &lake_var->energy.T[nidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            /* Read lake-specific variables */
            if (fscanf(init_state, " %hu", &lake_var->activenod) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->dz) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->surfdz) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->ldepth) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            for (node = 0; node <= lake_var->activenod; node++) {
                if (fscanf(init_state, " %lf",
                           &lake_var->surface[node]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            if (fscanf(init_state, " %lf", &lake_var->sarea) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->volume) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            for (node = 0; node < lake_var->activenod; node++) {
                if (fscanf(init_state, " %lf", &lake_var->temp[node]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            if (fscanf(init_state, " %lf", &lake_var->tempavg) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->areai) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->new_ice_area) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->ice_water_eq) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->hice) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->tempi) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->swe) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->surf_temp) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->pack_temp) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->coldcontent) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->surf_water) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->pack_water) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->SAlbedo) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
            if (fscanf(init_state, " %lf", &lake_var->sdepth) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
        }
    }

    // Check that soil moisture does not exceed maximum allowed
    for (veg = 0; veg <= Nveg; veg++) {
        for (band = 0; band < Nbands; band++) {
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                if (cell[veg][band].layer[lidx].moist >
                    soil_con->max_moist[lidx]) {
                    log_warn("Initial soil moisture (%f mm) exceeds "
                             "maximum (%f mm) in layer %zu for veg tile "
                             "%d and snow band %d.  Resetting to maximum.",
                             cell[veg][band].layer[lidx].moist,
                             soil_con->max_moist[lidx], lidx, veg, band);
                    for (frost_area = 0;
                         frost_area < options.Nfrost;
                         frost_area++) {
                        cell[veg][band].layer[lidx].ice[frost_area] *=
                            soil_con->max_moist[lidx] /
                            cell[veg][band].layer[lidx].moist;
                    }
                    cell[veg][band].layer[lidx].moist =
                        soil_con->max_moist[lidx];
                }

                for (frost_area = 0;
                     frost_area < options.Nfrost;
                     frost_area++) {
                    if (cell[veg][band].layer[lidx].ice[frost_area] >
                        cell[veg][band].layer[lidx].moist) {
                        cell[veg][band].layer[lidx].ice[frost_area] =
                            cell[veg][band].layer[lidx].moist;
                    }
                }
            }
        }

        // Override possible bad values of soil moisture under lake coming from state file
        // (ideally we wouldn't store these in the state file in the first place)
        if (options.LAKES && (int) veg == lake_con.lake_idx) {
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                lake_var->soil.layer[lidx].moist =
                    soil_con->max_moist[lidx];
                for (frost_area = 0;
                     frost_area < options.Nfrost;
                     frost_area++) {
                    if (lake_var->soil.layer[lidx].ice[frost_area] >
                        lake_var->soil.layer[lidx].moist) {
                        lake_var->soil.layer[lidx].ice[frost_area] =
                            lake_var->soil.layer[lidx].moist;
                    }
                }
            }
        }
    }
}
