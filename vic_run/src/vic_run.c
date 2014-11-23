#include <vic_def.h>
#include <vic_run.h>

veg_lib_struct *vic_run_veg_lib;

int
vic_run(int                  gridcell,
        int                  rec,
        atmos_data_struct   *atmos,
        all_vars_struct     *all_vars,
        dmy_struct          *dmy,
        global_param_struct *gp,
        lake_con_struct     *lake_con,
        soil_con_struct     *soil_con,
        veg_con_struct      *veg_con,
        veg_lib_struct      *veg_lib,
        veg_hist_struct     *veg_hist)
/**********************************************************************
        vic_run

   This subroutine controls the model core, it solves both the energy
   and water balance models, as well as frozen soils.
**********************************************************************/
{
    extern option_struct options;
    char                 overstory;
    int                  j, p;
    int                  lidx;
    int                  iveg;
    int                  Nveg;
    int                  veg_class;
    int                  band;
    int                  Nbands;
    int                  ErrorFlag;
    double               out_prec[2 * MAX_BANDS];
    double               out_rain[2 * MAX_BANDS];
    double               out_snow[2 * MAX_BANDS];
    double               dp;
    double               ice0[MAX_BANDS];
    double               moist0[MAX_BANDS];
    double               surf_atten;
    double               wind_h;
    double               height;
    double               displacement[3];
    double               roughness[3];
    double               ref_height[3];
    double             **aero_resist;
    double               Cv;
    double               Le;
    double               Melt[2 * MAX_BANDS];
    double               bare_albedo;
    double               snow_inflow[MAX_BANDS];
    double               rainonly;
    double               sum_runoff;
    double               sum_baseflow;
    double               tmp_wind[3];
    double               gauge_correction[2];
    float                lag_one;
    float                sigma_slope;
    float                fetch;
    int                  pet_veg_class;
    double               lakefrac;
    double               fraci;
    double               wetland_runoff;
    double               wetland_baseflow;
    double               snowprec;
    double               rainprec;
    int                  cidx;
    lake_var_struct     *lake_var;
    cell_data_struct   **cell;
    veg_var_struct     **veg_var;
    energy_bal_struct  **energy;
    snow_data_struct   **snow;

    // assign vic_run_veg_lib to veg_lib, so that the veg_lib for the correct
    // grid cell is used within vic_run. For simplicity sake, use vic_run_veg_lib
    // everywhere within vic_run
    vic_run_veg_lib = veg_lib;

    /* Allocate aero_resist array */
    aero_resist = (double**)calloc(N_PET_TYPES + 1, sizeof(double*));
    for (p = 0; p < N_PET_TYPES + 1; p++) {
        aero_resist[p] = (double*)calloc(3, sizeof(double));
    }

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

    /* Assign current veg albedo and LAI */
    // Loop over vegetated tiles
    for (iveg = 0; iveg < Nveg; iveg++) {
        veg_class = veg_con[iveg].veg_class;
        if (veg_hist[iveg].vegcover[0] < MIN_VEGCOVER) {
            veg_hist[iveg].vegcover[0] = MIN_VEGCOVER;
        }
        for (band = 0; band < Nbands; band++) {
            veg_var[iveg][band].vegcover = veg_hist[iveg].vegcover[0];
            veg_var[iveg][band].albedo = veg_hist[iveg].albedo[0];
            veg_var[iveg][band].LAI = veg_hist[iveg].LAI[0];
            // Convert LAI from global to local
            veg_var[iveg][band].LAI /= veg_var[iveg][band].vegcover;
            veg_var[iveg][band].Wdew /= veg_var[iveg][band].vegcover;
            veg_var[iveg][band].Wdmax = veg_var[iveg][band].LAI *
                                        LAI_WATER_FACTOR;
            snow[iveg][band].snow_canopy /= veg_var[iveg][band].vegcover;
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
                    veg_var[iveg][band].rc = HUGE_RESIST;
                }
            }

            /** Assign wind_h **/
            /** Note: this is ignored below **/
            wind_h = vic_run_veg_lib[veg_class].wind_h;

            /** Compute Surface Attenuation due to Vegetation Coverage **/
            surf_atten = (1 - veg_var[iveg][0].vegcover) * 1.0 +
                         veg_var[iveg][0].vegcover *
                         exp(-vic_run_veg_lib[veg_class].rad_atten *
                             veg_var[iveg][0].LAI);

            /* Initialize soil thermal properties for the top two layers */
            prepare_full_energy(iveg, Nveg, options.Nnode, all_vars, soil_con,
                                moist0, ice0);

            /** Compute Bare (free of snow) Albedo **/
            bare_albedo = veg_var[iveg][0].albedo;

            /*************************************
               Compute the aerodynamic resistance
               for current veg cover and various
               types of potential evap
            *************************************/

            /* Loop over types of potential evap, plus current veg */
            /* Current veg will be last */
            for (p = 0; p < N_PET_TYPES + 1; p++) {
                /* Initialize wind speeds */
                tmp_wind[0] = atmos->wind[NR];
                tmp_wind[1] = -999.;
                tmp_wind[2] = -999.;

                /* Set surface descriptive variables */
                if (p < N_PET_TYPES_NON_NAT) {
                    pet_veg_class = vic_run_veg_lib[0].NVegLibTypes + p;
                }
                else {
                    pet_veg_class = veg_class;
                }
                displacement[0] =
                    vic_run_veg_lib[pet_veg_class].displacement[dmy[rec].month -
                                                                1];
                roughness[0] =
                    vic_run_veg_lib[pet_veg_class].roughness[dmy[rec].month -
                                                             1];
                overstory = vic_run_veg_lib[pet_veg_class].overstory;
                if (p >= N_PET_TYPES_NON_NAT) {
                    if (roughness[0] == 0) {
                        roughness[0] = soil_con->rough;
                    }
                }

                /* Estimate vegetation height */
                height = calc_veg_height(displacement[0]);

                /* Estimate reference height */
                if (displacement[0] < wind_h) {
                    ref_height[0] = wind_h;
                }
                else {
                    ref_height[0] = displacement[0] + wind_h + roughness[0];
                }

                /* Compute aerodynamic resistance over various surface types */
                /* Do this not only for current veg but also all types of PET */
                ErrorFlag = CalcAerodynamic(overstory, height,
                                            vic_run_veg_lib[pet_veg_class].trunk_ratio,
                                            soil_con->snow_rough, soil_con->rough,
                                            vic_run_veg_lib[pet_veg_class].wind_atten,
                                            aero_resist[p], tmp_wind,
                                            displacement, ref_height,
                                            roughness);
                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }

            /* Initialize final aerodynamic resistance values */
            for (band = 0; band < Nbands; band++) {
                if (soil_con->AreaFract[band] > 0) {
                    cell[iveg][band].aero_resist[0] =
                        aero_resist[N_PET_TYPES][0];
                    cell[iveg][band].aero_resist[1] =
                        aero_resist[N_PET_TYPES][1];
                }
            }

            /******************************
               Compute nitrogen scaling factors and initialize other veg vars
            ******************************/
            if (options.CARBON && iveg < Nveg) {
                for (band = 0; band < Nbands; band++) {
                    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
                        veg_var[iveg][band].rsLayer[cidx] = HUGE_RESIST;
                    }
                    veg_var[iveg][band].aPAR = 0;
                    if (dmy->hour == 0) {
                        calc_Nscale_factors(
                            vic_run_veg_lib[veg_class].NscaleFlag,
                            veg_con[iveg].CanopLayerBnd,
                            vic_run_veg_lib[veg_class].LAI[dmy[
                                                               rec].month - 1],
                            soil_con->lat,
                            soil_con->lng,
                            soil_con->time_zone_lng,
                            dmy[rec],
                            veg_var[iveg][band].NscaleFactor);
                    }
                    if (dmy[rec].month == 1 && dmy[rec].day == 1) {
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
                    for (p = 0; p < N_PET_TYPES; p++) {
                        cell[iveg][band].pot_evap[p] = 0;
                    }

                    ErrorFlag = surface_fluxes(overstory, bare_albedo, height,
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
                                               Nbands,
                                               options.Nlayer, Nveg, band, dp,
                                               iveg, rec, veg_class,
                                               atmos, dmy,
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
                             1000 - soil_con->Wpwp[lidx]);
                    }
                    cell[iveg][band].wetness /= options.Nlayer;
                } /** End non-zero area band **/
            } /** End Loop Through Elevation Bands **/
        } /** end non-zero area veg tile **/
    } /** end of vegetation loop **/

    /* Convert LAI back to global */
    if (rec >= 0) {
        for (iveg = 0; iveg < Nveg; iveg++) {
            for (band = 0; band < Nbands; band++) {
                veg_var[iveg][band].LAI *= veg_var[iveg][band].vegcover;
                veg_var[iveg][band].Wdmax *= veg_var[iveg][band].vegcover;
            }
        }
    }

    for (p = 0; p < N_PET_TYPES + 1; p++) {
        free((char *)aero_resist[p]);
    }
    free((char *)aero_resist);

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
             wetland_runoff) * soil_con->cell_area * 0.001;                                               // m3
        lake_var->baseflow_in =
            (sum_baseflow * lake_con->rpercent +
             wetland_baseflow) * soil_con->cell_area * 0.001;                                                 // m3
        lake_var->channel_in = atmos->channel_in[NR] * soil_con->cell_area *
                               0.001;                                        // m3
        lake_var->prec = atmos->prec[NR] * lake_var->sarea * 0.001; // m3
        rainonly = calc_rainonly(atmos->air_temp[NR], atmos->prec[NR],
                                 gp->MAX_SNOW_TEMP, gp->MIN_RAIN_TEMP);
        if ((int)rainonly == ERROR) {
            return(ERROR);
        }

        /**********************************************************************
           Solve the energy budget for the lake.
        **********************************************************************/

        snowprec = gauge_correction[SNOW] * (atmos->prec[NR] - rainonly);
        rainprec = gauge_correction[SNOW] * rainonly;
        atmos->out_prec += (snowprec + rainprec) * lake_con->Cl[0] * lakefrac;
        atmos->out_rain += rainprec * lake_con->Cl[0] * lakefrac;
        atmos->out_snow += snowprec * lake_con->Cl[0] * lakefrac;

        ErrorFlag = solve_lake(snowprec, rainprec, atmos->air_temp[NR],
                               atmos->wind[NR], atmos->vp[NR] / 1000.,
                               atmos->shortwave[NR], atmos->longwave[NR],
                               atmos->vpd[NR] / 1000.,
                               atmos->pressure[NR] / 1000.,
                               atmos->density[NR], lake_var, *lake_con,
                               *soil_con, gp->dt, rec, gp->wind_h, dmy[rec],
                               fraci);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }

        /**********************************************************************
           Solve the water budget for the lake.
        **********************************************************************/

        ErrorFlag = water_balance(lake_var, *lake_con, gp->dt, all_vars, rec,
                                  iveg, band, lakefrac, *soil_con, *veg_con);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    } // end if (options.LAKES && lake_con->lake_idx >= 0)

    return (0);
}
