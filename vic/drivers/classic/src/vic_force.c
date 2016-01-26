/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes atmospheric variables for both the model time step,
 * and the time step used by the snow algorithm (if different).
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
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize atmospheric variables for the model and snow time steps.
 *****************************************************************************/
void
vic_force(atmos_data_struct *atmos,
          dmy_struct        *dmy,
          FILE             **infile,
          veg_con_struct    *veg_con,
          veg_hist_struct  **veg_hist,
          soil_con_struct   *soil_con)
{
    extern option_struct       options;
    extern param_set_struct    param_set;
    extern global_param_struct global_param;
    extern parameters_struct   param;
    extern size_t              NR, NF;

    size_t                     i;
    size_t                     j;
    size_t                     v;
    size_t                     rec;
    size_t                     uidx;
    double                     t_offset;
    double                   **forcing_data;
    double                  ***veg_hist_data;
    double                     avgJulyAirTemp;
    double                    *Tfactor;
    bool                      *AboveTreeLine;

    /*******************************
       Check that required inputs were supplied
    *******************************/

    if (!param_set.TYPE[AIR_TEMP].SUPPLIED) {
        log_err("Air temperature must be supplied as a forcing");
    }
    if (!param_set.TYPE[PREC].SUPPLIED) {
        log_err("Precipitation must be supplied as a forcing");
    }
    if (!param_set.TYPE[SHORTWAVE].SUPPLIED) {
        log_err("Downward shortwave radiation must be supplied as a forcing");
    }
    if (!param_set.TYPE[LONGWAVE].SUPPLIED) {
        log_err("Downward longwave radiation must be supplied as a forcing");
    }
    if (!param_set.TYPE[PRESSURE].SUPPLIED) {
        log_err("Atmospheric pressure must be supplied as a forcing");
    }
    if (!param_set.TYPE[VP].SUPPLIED) {
        log_err("Vapor ressure must be supplied as a forcing");
    }
    if (!param_set.TYPE[WIND].SUPPLIED) {
        log_err("Wind speed must be supplied as a forcing");
    }

    /*******************************
       Miscellaneous initialization
    *******************************/

    /* Assign local copies of some variables */
    avgJulyAirTemp = soil_con->avgJulyAirTemp;
    Tfactor = soil_con->Tfactor;
    AboveTreeLine = soil_con->AboveTreeLine;

    /* Assign N_ELEM for veg-dependent forcings */
    if (param_set.TYPE[ALBEDO].SUPPLIED) {
        param_set.TYPE[ALBEDO].N_ELEM = veg_con[0].vegetat_type_num;
    }
    if (param_set.TYPE[LAI_IN].SUPPLIED) {
        param_set.TYPE[LAI_IN].N_ELEM = veg_con[0].vegetat_type_num;
    }
    if (param_set.TYPE[VEGCOVER].SUPPLIED) {
        param_set.TYPE[VEGCOVER].N_ELEM = veg_con[0].vegetat_type_num;
    }

    /*******************************
       read in meteorological data
    *******************************/

    forcing_data = read_forcing_data(infile, global_param, &veg_hist_data);

    log_info("Read meteorological forcing file");

    /****************************************************
       Variables in the atmos_data structure
    ****************************************************/

    t_offset = Tfactor[0];
    for (i = 1; i < options.SNOW_BAND; i++) {
        if (Tfactor[i] < t_offset) {
            t_offset = Tfactor[i];
        }
    }

    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (i = 0; i < NF; i++) {
            uidx = rec * NF + i;
            // temperature in Celsius
            atmos[rec].air_temp[i] = forcing_data[AIR_TEMP][uidx];
            // precipitation in mm/period
            atmos[rec].prec[i] = forcing_data[PREC][uidx];
            // downward shortwave in W/m2
            atmos[rec].shortwave[i] = forcing_data[SHORTWAVE][uidx];
            // downward longwave in W/m2
            atmos[rec].longwave[i] = forcing_data[LONGWAVE][uidx];
            // pressure in kPa
            atmos[rec].pressure[i] = forcing_data[PRESSURE][uidx] * PA_PER_KPA;
            // vapor pressure in kPa
            atmos[rec].vp[i] = forcing_data[VP][uidx] * PA_PER_KPA;
            // vapor pressure deficit in kPa
            atmos[rec].vpd[i] = svp(atmos[rec].air_temp[i]) - atmos[rec].vp[i];
            // air density
            atmos[rec].density[i] = air_density(atmos[rec].air_temp[i],
                                                atmos[rec].pressure[i]);
            // wind speed in m/s
            atmos[rec].wind[i] = forcing_data[WIND][uidx];
            // snow flag
            atmos[rec].snowflag[i] = will_it_snow(&(atmos[rec].air_temp[i]),
                                                  t_offset,
                                                  param.SNOW_MAX_SNOW_TEMP,
                                                  &(atmos[rec].prec[i]), 1);
            // Optional inputs
            if (options.LAKES) {
                // Channel inflow from upstream (into lake)
                if (param_set.TYPE[CHANNEL_IN].SUPPLIED) {
                    atmos[rec].channel_in[i] = forcing_data[CHANNEL_IN][uidx];
                }
                else {
                    atmos[rec].channel_in[i] = 0;
                }
            }
            if (options.CARBON) {
                // Atmospheric CO2 concentration
                atmos[rec].Catm[i] = forcing_data[CATM][uidx];
                // Fraction of shortwave that is direct
                atmos[rec].fdir[i] = forcing_data[FDIR][uidx];
                // photosynthetically active radiation
                atmos[rec].par[i] = forcing_data[PAR][uidx];
            }
        }
        if (NF > 1) {
            atmos[rec].air_temp[NR] = average(atmos[rec].air_temp, NF);
            // For precipitation put total
            atmos[rec].prec[NR] = average(atmos[rec].prec, NF) * NF;
            atmos[rec].shortwave[NR] = average(atmos[rec].shortwave, NF);
            atmos[rec].longwave[NR] = average(atmos[rec].longwave, NF);
            atmos[rec].pressure[NR] = average(atmos[rec].pressure, NF);
            atmos[rec].vp[NR] = average(atmos[rec].vp, NF);
            atmos[rec].vpd[NR] = average(atmos[rec].vpd, NF);
            atmos[rec].density[NR] = average(atmos[rec].density, NF);
            atmos[rec].wind[NR] = average(atmos[rec].wind, NF);
            atmos[rec].snowflag[NR] = false;
            for (i = 0; i < NF; i++) {
                if (atmos[rec].snowflag[i] == true) {
                    atmos[rec].snowflag[NR] = true;
                }
            }
            if (options.LAKES) {
                atmos[rec].channel_in[NR] =
                    average(atmos[rec].channel_in, NF) * NF;
            }
            if (options.CARBON) {
                atmos[rec].Catm[NR] = average(atmos[rec].Catm, NF);
                atmos[rec].fdir[NR] = average(atmos[rec].fdir, NF);
                atmos[rec].par[NR] = average(atmos[rec].par, NF);
            }
        }
    }

    /****************************************************
       Variables in the veg_hist structure
    ****************************************************/

    /* First, assign default climatology */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
            for (i = 0; i < NF; i++) {
                veg_hist[rec][v].albedo[i] =
                    veg_con[v].albedo[dmy[rec].month - 1];
                veg_hist[rec][v].LAI[i] =
                    veg_con[v].LAI[dmy[rec].month - 1];
                veg_hist[rec][v].vegcover[i] =
                    veg_con[v].vegcover[dmy[rec].month - 1];
            }
        }
    }

    /* Next, overwrite with veg_hist values, validate, and average */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
            for (i = 0; i < NF; i++) {
                uidx = rec * NF + i;
                if (param_set.TYPE[ALBEDO].SUPPLIED &&
                    options.ALB_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[ALBEDO][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].albedo[i] =
                            veg_hist_data[ALBEDO][v][uidx];
                    }
                }
                if (param_set.TYPE[LAI_IN].SUPPLIED &&
                    options.LAI_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[LAI_IN][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].LAI[i] =
                            veg_hist_data[LAI_IN][v][uidx];
                    }
                }
                if (param_set.TYPE[VEGCOVER].SUPPLIED &&
                    options.VEGCOVER_SRC == FROM_VEGHIST) {
                    if (veg_hist_data[VEGCOVER][v][uidx] != NODATA_VH) {
                        veg_hist[rec][v].vegcover[i] =
                            veg_hist_data[VEGCOVER][v][uidx];
                    }
                }
                // Check on vegcover
                if (veg_hist[rec][v].vegcover[i] < MIN_VEGCOVER) {
                    log_warn(
                        "rec %zu, veg %zu substep %zu vegcover %f < minimum of %f; setting = %f\n", rec, v, i,
                        veg_hist[rec][v].vegcover[i], MIN_VEGCOVER,
                        MIN_VEGCOVER);
                    veg_hist[rec][v].vegcover[i] = MIN_VEGCOVER;
                }
            }
            if (NF > 1) {
                veg_hist[rec][v].albedo[NR] = average(veg_hist[rec][v].albedo,
                                                      NF);
                veg_hist[rec][v].LAI[NR] = average(veg_hist[rec][v].LAI, NF);
                veg_hist[rec][v].vegcover[NR] = average(
                    veg_hist[rec][v].vegcover, NF);
            }
        }
    }

    /****************************************************
       Free forcing_data and veg_hist_data structures
    ****************************************************/

    for (i = 0; i < N_FORCING_TYPES; i++) {
        if (param_set.TYPE[i].SUPPLIED) {
            if (i != ALBEDO && i != LAI_IN && i != VEGCOVER) {
                free(forcing_data[i]);
            }
            else {
                for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                    free(veg_hist_data[i][j]);
                }
                free(veg_hist_data[i]);
            }
        }
    }
    free(forcing_data);
    free(veg_hist_data);

    /****************************************************
       Compute treeline based on July average temperature
    ****************************************************/

    if (options.COMPUTE_TREELINE) {
        if (!(options.JULY_TAVG_SUPPLIED && avgJulyAirTemp == -999)) {
            compute_treeline(atmos, dmy, avgJulyAirTemp, Tfactor,
                             AboveTreeLine);
        }
    }
}
