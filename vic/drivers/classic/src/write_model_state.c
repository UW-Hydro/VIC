/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine saves the model state at second 0 of the date defined in the
 * global control file using STATEDAY, STATEMONTH, and STATEYEAR.  The saved
 * files can then be used to initialize the model to the same state as when the
 * files were created.
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
 * @brief    Saves the model state.
 *****************************************************************************/
void
write_model_state(all_vars_struct *all_vars,
                  int              Nveg,
                  int              cellnum,
                  filep_struct    *filep,
                  soil_con_struct *soil_con)
{
    extern option_struct options;

    double               tmpval;
    int                  veg;
    int                  band;
    size_t               lidx;
    size_t               nidx;
    int                  Nbands;
    int                  Nbytes;
    size_t               frost_area;

    cell_data_struct   **cell;
    snow_data_struct   **snow;
    energy_bal_struct  **energy;
    veg_var_struct     **veg_var;
    lake_var_struct      lake_var;
    int                  node;

    Nbands = options.SNOW_BAND;

    cell = all_vars->cell;
    veg_var = all_vars->veg_var;
    snow = all_vars->snow;
    energy = all_vars->energy;
    lake_var = all_vars->lake_var;

    /* write cell information */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&cellnum, sizeof(int), 1, filep->statefile);
        fwrite(&Nveg, sizeof(int), 1, filep->statefile);
        fwrite(&Nbands, sizeof(int), 1, filep->statefile);
    }
    else {
        fprintf(filep->statefile, "%i %i %i", cellnum, Nveg, Nbands);
    }
    // This stores the number of bytes from after this value to the end
    // of the line.  DO NOT CHANGE unless you have changed the values
    // written to the state file.
    // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
    if (options.STATE_FORMAT == BINARY) {
        Nbytes = options.Nnode * sizeof(double) + // dz_node
                 options.Nnode * sizeof(double) + // Zsum_node
                 (Nveg + 1) * Nbands * 2 * sizeof(int) + // veg & band
                 (Nveg + 1) * Nbands * options.Nlayer * sizeof(double) + // soil moisture
                 (Nveg +
                  1) * Nbands * options.Nlayer * options.Nfrost *
                 sizeof(double) +                                                     // soil ice
                 Nveg * Nbands * sizeof(double); // dew
        if (options.CARBON) {
            /* Carbon-specific state vars */
            Nbytes += Nveg * Nbands * 5 * sizeof(double); // AnnualNPP, AnnualNPPPrev, and 3 soil carbon storages
        }
        Nbytes += (Nveg + 1) * Nbands * sizeof(int) + // last_snow
                  (Nveg + 1) * Nbands * sizeof(char) + // MELTING
                  (Nveg + 1) * Nbands * sizeof(double) * 9 + // other snow parameters
                  (Nveg + 1) * Nbands * options.Nnode * sizeof(double) + // soil temperatures
                  (Nveg + 1) * Nbands * sizeof(double) + // Tfoliage
                  (Nveg + 1) * Nbands * sizeof(double) + // energy.LongUnderOut
                  (Nveg + 1) * Nbands * sizeof(double); // energy.snow_flux
        if (options.LAKES) {
            /* Lake/wetland tiles have lake-specific state vars */
            Nbytes += sizeof(int) + // activenod
                      sizeof(double) + // dz
                      sizeof(double) + // surfdz
                      sizeof(double) + // ldepth
                      (lake_var.activenod + 1) * sizeof(double) + // surface
                      sizeof(double) + // sarea
                      sizeof(double) + // volume
                      lake_var.activenod * sizeof(double) + // temp
                      sizeof(double) + // tempavg
                      sizeof(double) + // areai
                      sizeof(double) + // new_ice_area
                      sizeof(double) + // ice_water_eq
                      sizeof(double) + // hice
                      sizeof(double) + // tempi
                      sizeof(double) + // swe
                      sizeof(double) + // surf_temp
                      sizeof(double) + // pack_temp
                      sizeof(double) + // coldcontent
                      sizeof(double) + // surf_water
                      sizeof(double) + // pack_water
                      sizeof(double) + // SAlbedo
                      sizeof(double) + // sdepth
                      options.Nlayer * sizeof(double) + // soil moisture
                      options.Nlayer * options.Nfrost * sizeof(double) + // soil ice
                      sizeof(int) + // last_snow
                      sizeof(char) + // MELTING
                      sizeof(double) * 9 + // other snow parameters
                      options.Nnode * sizeof(double) // soil temperatures
            ;
            if (options.CARBON) {
                /* Carbon-specific state vars */
                Nbytes += 3 * sizeof(double); // 3 soil carbon storages
            }
        }
        fwrite(&Nbytes, sizeof(int), 1, filep->statefile);
    }

    /* Write soil thermal node deltas */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&soil_con->dz_node[nidx], sizeof(double), 1,
                   filep->statefile);
        }
        else {
            fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                    soil_con->dz_node[nidx]);
        }
    }
    /* Write soil thermal node depths */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&soil_con->Zsum_node[nidx], sizeof(double), 1,
                   filep->statefile);
        }
        else {
            fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                    soil_con->Zsum_node[nidx]);
        }
    }
    if (options.STATE_FORMAT == ASCII) {
        fprintf(filep->statefile, "\n");
    }

    /* Output for all vegetation types */
    for (veg = 0; veg <= Nveg; veg++) {
        /* Output for all snow bands */
        for (band = 0; band < Nbands; band++) {
            /* Write cell identification information */
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&veg, sizeof(int), 1, filep->statefile);
                fwrite(&band, sizeof(int), 1, filep->statefile);
            }
            else {
                fprintf(filep->statefile, "%i %i", veg, band);
            }

            /* Write total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                tmpval = cell[veg][band].layer[lidx].moist;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT, tmpval);
                }
            }

            /* Write average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    tmpval = cell[veg][band].layer[lidx].ice[frost_area];
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                }
            }

            if (veg < Nveg) {
                /* Write dew storage */
                tmpval = veg_var[veg][band].Wdew;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT, tmpval);
                }
                if (options.CARBON) {
                    /* Write cumulative NPP */
                    tmpval = veg_var[veg][band].AnnualNPP;
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                    tmpval = veg_var[veg][band].AnnualNPPPrev;
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                    /* Write soil carbon storages */
                    tmpval = cell[veg][band].CLitter;
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                    tmpval = cell[veg][band].CInter;
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                    tmpval = cell[veg][band].CSlow;
                    if (options.STATE_FORMAT == BINARY) {
                        fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                    }
                    else {
                        fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                                tmpval);
                    }
                }
            }

            /* Write snow data */
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg][band].last_snow, sizeof(int), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].MELTING, sizeof(char), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].coverage, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].swq, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].surf_temp, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].surf_water, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].pack_temp, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].pack_water, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].density, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].coldcontent, sizeof(double), 1,
                       filep->statefile);
                fwrite(&snow[veg][band].snow_canopy, sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " %i %i "
                        ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                        ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                        ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                        ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                        ASCII_STATE_FLOAT_FMT,
                        snow[veg][band].last_snow, (int)snow[veg][band].MELTING,
                        snow[veg][band].coverage, snow[veg][band].swq,
                        snow[veg][band].surf_temp, snow[veg][band].surf_water,
                        snow[veg][band].pack_temp, snow[veg][band].pack_water,
                        snow[veg][band].density, snow[veg][band].coldcontent,
                        snow[veg][band].snow_canopy);
            }

            /* Write soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&energy[veg][band].T[nidx], sizeof(double), 1,
                           filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            energy[veg][band].T[nidx]);
                }
            }

            /* Write foliage temperature */
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&energy[veg][band].Tfoliage, sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        energy[veg][band].Tfoliage);
            }

            /* Write outgoing longwave from understory */
            /* TO-DO: this is a flux. Saving it to the state file is a temporary solution! */
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&energy[veg][band].LongUnderOut, sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        energy[veg][band].LongUnderOut);
            }

            /* Write thermal flux through the snow pack */
            /* TO-DO: this is a flux. Saving it to the state file is a temporary solution! */
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&energy[veg][band].snow_flux, sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        energy[veg][band].snow_flux);
            }

            if (options.STATE_FORMAT == ASCII) {
                fprintf(filep->statefile, "\n");
            }
        }
    }

    if (options.LAKES) {
        if (options.STATE_FORMAT == BINARY) {
            /* Write total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                fwrite(&lake_var.soil.layer[lidx].moist, sizeof(double), 1,
                       filep->statefile);
            }

            /* Write average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    fwrite(&lake_var.soil.layer[lidx].ice[frost_area],
                           sizeof(double), 1, filep->statefile);
                }
            }
            if (options.CARBON) {
                /* Write soil carbon storages */
                tmpval = lake_var.soil.CLitter;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " %f", tmpval);
                }
                tmpval = lake_var.soil.CInter;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " %f", tmpval);
                }
                tmpval = lake_var.soil.CSlow;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " %f", tmpval);
                }
            }

            /* Write snow data */
            fwrite(&lake_var.snow.last_snow, sizeof(int), 1, filep->statefile);
            fwrite(&lake_var.snow.MELTING, sizeof(char), 1, filep->statefile);
            fwrite(&lake_var.snow.coverage, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.swq, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.snow.surf_temp, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.surf_water, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.pack_temp, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.pack_water, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.density, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.snow.coldcontent, sizeof(double), 1,
                   filep->statefile);
            fwrite(&lake_var.snow.snow_canopy, sizeof(double), 1,
                   filep->statefile);

            /* Write soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                fwrite(&lake_var.energy.T[nidx], sizeof(double), 1,
                       filep->statefile);
            }

            /* Write lake-specific variables */
            fwrite(&lake_var.activenod, sizeof(int), 1, filep->statefile);
            fwrite(&lake_var.dz, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.surfdz, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.ldepth, sizeof(double), 1, filep->statefile);
            for (node = 0; node <= lake_var.activenod; node++) {
                fwrite(&lake_var.surface[node], sizeof(double), 1,
                       filep->statefile);
            }
            fwrite(&lake_var.sarea, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.volume, sizeof(double), 1, filep->statefile);
            for (node = 0; node < lake_var.activenod; node++) {
                fwrite(&lake_var.temp[node], sizeof(double), 1,
                       filep->statefile);
            }
            fwrite(&lake_var.tempavg, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.areai, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.new_ice_area, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.ice_water_eq, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.hice, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.tempi, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.swe, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.surf_temp, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.pack_temp, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.coldcontent, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.surf_water, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.pack_water, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.SAlbedo, sizeof(double), 1, filep->statefile);
            fwrite(&lake_var.sdepth, sizeof(double), 1, filep->statefile);
        }
        else {
            /* Write total soil moisture */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                fprintf(filep->statefile, " %f",
                        lake_var.soil.layer[lidx].moist);
            }

            /* Write average ice content */
            for (lidx = 0; lidx < options.Nlayer; lidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    fprintf(filep->statefile, " %f",
                            lake_var.soil.layer[lidx].ice[frost_area]);
                }
            }

            /* Write snow data */
            fprintf(filep->statefile, " %i %i %f %f %f %f %f %f %f %f %f",
                    lake_var.snow.last_snow, (int)lake_var.snow.MELTING,
                    lake_var.snow.coverage, lake_var.snow.swq,
                    lake_var.snow.surf_temp, lake_var.snow.surf_water,
                    lake_var.snow.pack_temp, lake_var.snow.pack_water,
                    lake_var.snow.density, lake_var.snow.coldcontent,
                    lake_var.snow.snow_canopy);

            /* Write soil thermal node temperatures */
            for (nidx = 0; nidx < options.Nnode; nidx++) {
                fprintf(filep->statefile, " %f", lake_var.energy.T[nidx]);
            }

            /* Write lake-specific variables */
            fprintf(filep->statefile, " %d", lake_var.activenod);
            fprintf(filep->statefile, " %f", lake_var.dz);
            fprintf(filep->statefile, " %f", lake_var.surfdz);
            fprintf(filep->statefile, " %f", lake_var.ldepth);
            for (node = 0; node <= lake_var.activenod; node++) {
                fprintf(filep->statefile, " %f", lake_var.surface[node]);
            }
            fprintf(filep->statefile, " %f", lake_var.sarea);
            fprintf(filep->statefile, " %f", lake_var.volume);
            for (node = 0; node < lake_var.activenod; node++) {
                fprintf(filep->statefile, " %f", lake_var.temp[node]);
            }
            fprintf(filep->statefile, " %f", lake_var.tempavg);
            fprintf(filep->statefile, " %f", lake_var.areai);
            fprintf(filep->statefile, " %f", lake_var.new_ice_area);
            fprintf(filep->statefile, " %f", lake_var.ice_water_eq);
            fprintf(filep->statefile, " %f", lake_var.hice);
            fprintf(filep->statefile, " %f", lake_var.tempi);
            fprintf(filep->statefile, " %f", lake_var.swe);
            fprintf(filep->statefile, " %f", lake_var.surf_temp);
            fprintf(filep->statefile, " %f", lake_var.pack_temp);
            fprintf(filep->statefile, " %f", lake_var.coldcontent);
            fprintf(filep->statefile, " %f", lake_var.surf_water);
            fprintf(filep->statefile, " %f", lake_var.pack_water);
            fprintf(filep->statefile, " %f", lake_var.SAlbedo);
            fprintf(filep->statefile, " %f", lake_var.sdepth);

            fprintf(filep->statefile, "\n");
        }
    }
    /* Force file to be written */
    fflush(filep->statefile);
}
