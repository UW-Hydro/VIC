/******************************************************************************
* @section DESCRIPTION
*
* This subroutine controls the model core, it solves both the energy and water
* balance models, as well as frozen soils.
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
******************************************************************************/

#include <vic_run.h>

veg_lib_struct *vic_run_veg_lib;

/******************************************************************************
* @brief        This subroutine controls the model core, it solves both the
*               energy and water balance models, as well as frozen soils.
******************************************************************************/
int
vic_run(atmos_data_struct   *atmos,
        all_vars_struct     *all_vars,
        dmy_struct          *dmy,
        global_param_struct *gp,
        lake_con_struct     *lake_con,
        soil_con_struct     *soil_con,
        veg_con_struct      *veg_con,
        veg_lib_struct      *veg_lib)
{
    extern option_struct     options;
    extern parameters_struct param;

    char                     overstory;
    int                      j;
    size_t                   lidx;
    unsigned short           iveg;
    size_t                   Nveg;
    unsigned short           veg_class;
    unsigned short           band;
    size_t                   Nbands;
    int                      ErrorFlag;
    double                   out_prec[2 * MAX_BANDS];
    double                   out_rain[2 * MAX_BANDS];
    double                   out_snow[2 * MAX_BANDS];
    double                   dp;
    double                   ice0[MAX_BANDS];
    double                   moist0[MAX_BANDS];
    double                   surf_atten;
    double                   wind_h;
    double                   height;
    double                   displacement[3];
    double                   roughness[3];
    double                   ref_height[3];
    double                  *aero_resist;
    double                   Cv;
    double                   Le;
    double                   Melt[2 * MAX_BANDS];
    double                   bare_albedo;
    double                   snow_inflow[MAX_BANDS];
    double                   rainonly;
    double                   sum_runoff;
    double                   sum_baseflow;
    double                   tmp_wind[3];
    double                   gauge_correction[2];
    double                   lag_one;
    double                   sigma_slope;
    double                   fetch;
    double                   lakefrac;
    double                   fraci;
    double                   wetland_runoff;
    double                   wetland_baseflow;
    double                   snowprec;
    double                   rainprec;
    size_t                   cidx;
    lake_var_struct         *lake_var;
    cell_data_struct       **cell;
    veg_var_struct         **veg_var;
    energy_bal_struct      **energy;
    snow_data_struct       **snow;

    // assign vic_run_veg_lib to veg_lib, so that the veg_lib for the correct
    // grid cell is used within vic_run. For simplicity sake, use vic_run_veg_lib
    // everywhere within vic_run
    vic_run_veg_lib = veg_lib;

    /* Allocate aero_resist array */
    aero_resist = calloc(3, sizeof(*aero_resist));

    /* set local pointers */
    cell = all_vars->cell;
    energy = all_vars->energy;
    lake_var = &all_vars->lake_var;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    Nbands = options.SNOW_BAND;

    /* Set number of vegetation types */
    Nveg = veg_con[0].vegetat_type_num;

    /** Set Damping Depth **/
    dp = soil_con->dp;

    /* Compute gauge undercatch correction factors
       - this assumes that the gauge is free of vegetation effects, so gauge
       correction is constant for the entire grid cell */
    if (options.CORRPREC && atmos->prec[NR] > 0) {
        correct_precip(gauge_correction, atmos->wind[NR], gp->wind_h,
                       soil_con->rough, soil_con->snow_rough);
    }
    else {
        gauge_correction[0] = 1;
        gauge_correction[1] = 1;
    }
    atmos->out_prec = 0;
    atmos->out_rain = 0;
    atmos->out_snow = 0;

    // Convert LAI from global to local
    for (iveg = 0; iveg < Nveg; iveg++) {
        veg_class = veg_con[iveg].veg_class;
        for (band = 0; band < Nbands; band++) {
            veg_var[iveg][band].LAI /= veg_var[iveg][band].fcanopy;
            veg_var[iveg][band].Wdew /= veg_var[iveg][band].fcanopy;
            veg_var[iveg][band].Wdmax = veg_var[iveg][band].LAI *
                                        param.VEG_LAI_WATER_FACTOR;
            snow[iveg][band].snow_canopy /= veg_var[iveg][band].fcanopy;
        }
    }

    /**************************************************
       Solve Energy and/or Water Balance for Each
       Vegetation Type
    **************************************************/
    for (iveg = 0; iveg <= Nveg; iveg++) {
        /** Solve Veg Type only if Coverage Greater than 0% **/
        if (veg_con[iveg].Cv > 0.0) {
            Cv = veg_con[iveg].Cv;
            Nbands = options.SNOW_BAND;

            /** Lake-specific processing **/
            if (veg_con[iveg].LAKE) {
                /* Update areai to equal new ice area from previous time step. */
                lake_var->areai = lake_var->new_ice_area;

                /* Compute lake fraction and ice-covered fraction */
                if (lake_var->areai < 0) {
                    lake_var->areai = 0;
                }
                if (lake_var->sarea > 0) {
                    fraci = lake_var->areai / lake_var->sarea;
                    if (fraci > 1.0) {
                        fraci = 1.0;
                    }
                }
                else {
                    fraci = 0.0;
                }
                lakefrac = lake_var->sarea / lake_con->basin[0];

                Nbands = 1;
                Cv *= (1 - lakefrac);

                if (Cv == 0) {
                    continue;
                }
            }

            /**************************************************
               Initialize Model Parameters
            **************************************************/

            for (band = 0; band < Nbands; band++) {
                if (soil_con->AreaFract[band] > 0) {
                    /* Initialize energy balance variables */
                    energy[iveg][band].shortwave = 0;
                    energy[iveg][band].longwave = 0.;

                    /* Initialize snow variables */
                    snow[iveg][band].vapor_flux = 0.;
                    snow[iveg][band].canopy_vapor_flux = 0.;
                    snow_inflow[band] = 0.;
                    Melt[band * 2] = 0.;
                }
            }

            /* Initialize precipitation storage */
            for (j = 0; j < 2 * MAX_BANDS; j++) {
                out_prec[j] = 0;
                out_rain[j] = 0;
                out_snow[j] = 0;
            }

            /** Define vegetation class number **/
            veg_class = veg_con[iveg].veg_class;

            /** Initialize other veg vars **/
            if (iveg < Nveg) {
                for (band = 0; band < Nbands; band++) {
                    veg_var[iveg][band].rc = param.HUGE_RESIST;
                }
            }

            /** Assign wind_h **/
            /** Note: this is ignored below **/
            wind_h = vic_run_veg_lib[veg_class].wind_h;

            /** Compute Surface Attenuation due to Vegetation Coverage **/
            surf_atten = (1 - veg_var[iveg][0].fcanopy) * 1.0 +
                         veg_var[iveg][0].fcanopy *
                         exp(-vic_run_veg_lib[veg_class].rad_atten *
                             veg_var[iveg][0].LAI);

            /* Initialize soil thermal properties for the top two layers */
            prepare_full_energy(iveg, all_vars, soil_con, moist0, ice0);

            /** Compute Bare (free of snow) Albedo **/
            if (iveg != Nveg) {
                bare_albedo = veg_var[iveg][0].albedo;
            }
            else {
                bare_albedo = param.ALBEDO_BARE_SOIL;
            }

            /*************************************
               Compute the aerodynamic resistance
            *************************************/

            /* Initialize wind speeds */
            tmp_wind[0] = atmos->wind[NR];
            tmp_wind[1] = MISSING;
            tmp_wind[2] = MISSING;

            /* Set surface descriptive variables */
            displacement[0] =
                vic_run_veg_lib[veg_class].displacement[dmy->month - 1];
            roughness[0] =
                vic_run_veg_lib[veg_class].roughness[dmy->month - 1];
            if (roughness[0] == 0) {
                roughness[0] = soil_con->rough;
            }
            overstory = vic_run_veg_lib[veg_class].overstory;

            /* Estimate vegetation height */
            height = calc_veg_height(displacement[0]);

            /* Estimate reference height */
            if (displacement[0] < wind_h) {
                ref_height[0] = wind_h;
            }
            else {
                ref_height[0] = displacement[0] + wind_h + roughness[0];
            }

            /* Compute aerodynamic resistance */
            ErrorFlag = CalcAerodynamic(overstory, height,
                                        vic_run_veg_lib[veg_class].trunk_ratio,
                                        soil_con->snow_rough, soil_con->rough,
                                        vic_run_veg_lib[veg_class].wind_atten,
                                        aero_resist, tmp_wind,
                                        displacement, ref_height,
                                        roughness);
            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            /* Initialize final aerodynamic resistance values */
            for (band = 0; band < Nbands; band++) {
                if (soil_con->AreaFract[band] > 0) {
                    cell[iveg][band].aero_resist[0] =
                        aero_resist[0];
                    cell[iveg][band].aero_resist[1] =
                        aero_resist[1];
                }
            }

            /******************************
               Compute nitrogen scaling factors and initialize other veg vars
            ******************************/
            if (options.CARBON && iveg < Nveg) {
                for (band = 0; band < Nbands; band++) {
                    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
                        veg_var[iveg][band].rsLayer[cidx] = param.HUGE_RESIST;
                    }
                    veg_var[iveg][band].aPAR = 0;
                    if (dmy->dayseconds == 0) {
                        calc_Nscale_factors(
                            vic_run_veg_lib[veg_class].NscaleFlag,
                            veg_con[iveg].CanopLayerBnd,
                            vic_run_veg_lib[veg_class].LAI[dmy->month - 1],
                            soil_con->lat,
                            soil_con->lng,
                            soil_con->time_zone_lng,
                            dmy->day_in_year,
                            veg_var[iveg][band].NscaleFactor);
                    }
                    if (dmy->month == 1 && dmy->day == 1) {
                        veg_var[iveg][band].AnnualNPPPrev =
                            veg_var[iveg][band].AnnualNPP;
                        veg_var[iveg][band].AnnualNPP = 0;
                    }
                }
            }

            /******************************
               Solve ground surface fluxes
            ******************************/

            for (band = 0; band < Nbands; band++) {
                if (soil_con->AreaFract[band] > 0) {
                    lag_one = veg_con[iveg].lag_one;
                    sigma_slope = veg_con[iveg].sigma_slope;
                    fetch = veg_con[iveg].fetch;

                    /* Initialize pot_evap */
                    cell[iveg][band].pot_evap = 0;

                    ErrorFlag = surface_fluxes(overstory, bare_albedo,
                                               ice0[band], moist0[band],
                                               surf_atten, &(Melt[band * 2]),
                                               &Le,
                                               aero_resist,
                                               displacement, gauge_correction,
                                               &out_prec[band * 2],
                                               &out_rain[band * 2],
                                               &out_snow[band * 2],
                                               ref_height, roughness,
                                               &snow_inflow[band],
                                               tmp_wind, veg_con[iveg].root,
                                               options.Nlayer, Nveg, band, dp,
                                               iveg, veg_class, atmos, dmy,
                                               &(energy[iveg][band]), gp,
                                               &(cell[iveg][band]),
                                               &(snow[iveg][band]),
                                               soil_con, &(veg_var[iveg][band]),
                                               lag_one, sigma_slope, fetch,
                                               veg_con[iveg].CanopLayerBnd);

                    if (ErrorFlag == ERROR) {
                        return (ERROR);
                    }

                    atmos->out_prec +=
                        out_prec[band * 2] * Cv * soil_con->AreaFract[band];
                    atmos->out_rain +=
                        out_rain[band * 2] * Cv * soil_con->AreaFract[band];
                    atmos->out_snow +=
                        out_snow[band * 2] * Cv * soil_con->AreaFract[band];

                    /********************************************************
                       Compute soil wetness and root zone soil moisture
                    ********************************************************/
                    cell[iveg][band].rootmoist = 0;
                    cell[iveg][band].wetness = 0;
                    for (lidx = 0; lidx < options.Nlayer; lidx++) {
                        if (veg_con->root[lidx] > 0) {
                            cell[iveg][band].rootmoist +=
                                cell[iveg][band].layer[lidx].moist;
                        }
                        cell[iveg][band].wetness +=
                            (cell[iveg][band].layer[lidx].moist -
                             soil_con->Wpwp[lidx]) /
                            (soil_con->porosity[lidx] * soil_con->depth[lidx] *
                             MM_PER_M - soil_con->Wpwp[lidx]);
                    }
                    cell[iveg][band].wetness /= options.Nlayer;
                } /** End non-zero area band **/
            } /** End Loop Through Elevation Bands **/
        } /** end non-zero area veg tile **/
    } /** end of vegetation loop **/

    /* Convert LAI back to global */
    for (iveg = 0; iveg < Nveg; iveg++) {
        for (band = 0; band < Nbands; band++) {
            veg_var[iveg][band].LAI *= veg_var[iveg][band].fcanopy;
            veg_var[iveg][band].Wdmax *= veg_var[iveg][band].fcanopy;
        }
    }

    free((char *) aero_resist);

    /****************************
       Run Lake Model
    ****************************/

    /** Compute total runoff and baseflow for all vegetation types
        within each snowband. **/
    if (options.LAKES && lake_con->lake_idx >= 0) {
        wetland_runoff = wetland_baseflow = 0;
        sum_runoff = sum_baseflow = 0;

        // Loop through all vegetation tiles
        for (iveg = 0; iveg <= Nveg; iveg++) {
            /** Solve Veg Tile only if Coverage Greater than 0% **/
            if (veg_con[iveg].Cv > 0.) {
                Cv = veg_con[iveg].Cv;
                Nbands = options.SNOW_BAND;
                if (veg_con[iveg].LAKE) {
                    Cv *= (1 - lakefrac);
                    Nbands = 1;
                }

                // Loop through snow elevation bands
                for (band = 0; band < Nbands; band++) {
                    if (soil_con->AreaFract[band] > 0) {
                        if (veg_con[iveg].LAKE) {
                            wetland_runoff += (cell[iveg][band].runoff *
                                               Cv * soil_con->AreaFract[band]);
                            wetland_baseflow += (cell[iveg][band].baseflow *
                                                 Cv *
                                                 soil_con->AreaFract[band]);
                            cell[iveg][band].runoff = 0;
                            cell[iveg][band].baseflow = 0;
                        }
                        else {
                            sum_runoff += (cell[iveg][band].runoff *
                                           Cv * soil_con->AreaFract[band]);
                            sum_baseflow += (cell[iveg][band].baseflow *
                                             Cv * soil_con->AreaFract[band]);
                            cell[iveg][band].runoff *= (1 - lake_con->rpercent);
                            cell[iveg][band].baseflow *=
                                (1 - lake_con->rpercent);
                        }
                    }
                }
            }
        }

        /** Run lake model **/
        iveg = lake_con->lake_idx;
        band = 0;
        lake_var->runoff_in =
            (sum_runoff * lake_con->rpercent +
             wetland_runoff) * soil_con->cell_area / MM_PER_M;                                               // m3
        lake_var->baseflow_in =
            (sum_baseflow * lake_con->rpercent +
             wetland_baseflow) * soil_con->cell_area / MM_PER_M;                                                 // m3
        lake_var->channel_in = atmos->channel_in[NR] * soil_con->cell_area /
                               MM_PER_M;                                        // m3
        lake_var->prec = atmos->prec[NR] * lake_var->sarea / MM_PER_M; // m3
        rainonly = calc_rainonly(atmos->air_temp[NR], atmos->prec[NR],
                                 param.SNOW_MAX_SNOW_TEMP,
                                 param.SNOW_MIN_RAIN_TEMP);
        if ((int) rainonly == ERROR) {
            return(ERROR);
        }

        /**********************************************************************
           Solve the energy budget for the lake.
        **********************************************************************/

        snowprec = gauge_correction[SNOW] * (atmos->prec[NR] - rainonly);
        rainprec = gauge_correction[SNOW] * rainonly;
        Cv = veg_con[iveg].Cv * lakefrac;
        atmos->out_prec += (snowprec + rainprec) * Cv;
        atmos->out_rain += rainprec * Cv;
        atmos->out_snow += snowprec * Cv;

        ErrorFlag = solve_lake(snowprec, rainprec, atmos->air_temp[NR],
                               atmos->wind[NR], atmos->vp[NR] / PA_PER_KPA,
                               atmos->shortwave[NR], atmos->longwave[NR],
                               atmos->vpd[NR] / PA_PER_KPA,
                               atmos->pressure[NR] / PA_PER_KPA,
                               atmos->density[NR], lake_var,
                               *soil_con, gp->dt, gp->wind_h, *dmy,
                               fraci);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        /**********************************************************************
           Solve the water budget for the lake.
        **********************************************************************/

        ErrorFlag = water_balance(lake_var, *lake_con, gp->dt, all_vars,
                                  iveg, band, lakefrac, *soil_con, *veg_con);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    } // end if (options.LAKES && lake_con->lake_idx >= 0)

    return (0);
}
