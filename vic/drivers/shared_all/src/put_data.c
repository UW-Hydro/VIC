/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine converts data units, and stores finalized values in an array
 * for later output to the output files.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine converts data units, and stores finalized values in
 *           an array for later output to the output files.
 *****************************************************************************/
void
put_data(all_vars_struct   *all_vars,
         force_data_struct *force,
         soil_con_struct   *soil_con,
         veg_con_struct    *veg_con,
         veg_lib_struct    *veg_lib,
         lake_con_struct   *lake_con,
         double           **out_data,
         save_data_struct  *save_data,
         timer_struct      *timer)
{
    extern global_param_struct global_param;
    extern option_struct       options;
    extern parameters_struct   param;

    size_t                     veg;
    size_t                     index;
    size_t                     band;
    size_t                     Nbands;
    bool                       overstory;
    bool                       HasVeg;
    bool                       IsWet;
    bool                      *AboveTreeLine;
    double                    *AreaFract;
    double                    *depth;
    double                    *frost_fract;
    double                     frost_slope;
    double                     Cv;
    double                     Clake;
    double                     cv_baresoil;
    double                     cv_veg;
    double                     cv_overstory;
    double                     cv_snow;
    double                     inflow;
    double                     outflow;
    double                     storage;
    double                    *TreeAdjustFactor;
    double                     ThisAreaFract;
    double                     ThisTreeAdjust;
    size_t                     i;
    double                     dt_sec;

    cell_data_struct         **cell;
    energy_bal_struct        **energy;
    lake_var_struct            lake_var;
    snow_data_struct         **snow;
    veg_var_struct           **veg_var;

    cell = all_vars->cell;
    energy = all_vars->energy;
    lake_var = all_vars->lake_var;
    snow = all_vars->snow;
    veg_var = all_vars->veg_var;

    AboveTreeLine = soil_con->AboveTreeLine;
    AreaFract = soil_con->AreaFract;
    depth = soil_con->depth;
    frost_fract = soil_con->frost_fract;
    frost_slope = soil_con->frost_slope;
    dt_sec = global_param.dt;

    // Compute treeline adjustment factors
    TreeAdjustFactor = calloc(options.SNOW_BAND, sizeof(*TreeAdjustFactor));
    check_alloc_status(TreeAdjustFactor, "Memory allocation error.");
    for (band = 0; band < options.SNOW_BAND; band++) {
        if (AboveTreeLine[band]) {
            Cv = 0;
            for (veg = 0; veg < veg_con[0].vegetat_type_num; veg++) {
                if (veg_lib[veg_con[veg].veg_class].overstory) {
                    if (options.LAKES && veg_con[veg].LAKE) {
                        if (band == 0) {
                            // Fraction of tile that is flooded
                            Clake = lake_var.sarea / lake_con->basin[0];
                            Cv += veg_con[veg].Cv * (1 - Clake);
                        }
                    }
                    else {
                        Cv += veg_con[veg].Cv;
                    }
                }
            }
            TreeAdjustFactor[band] = 1. / (1. - Cv);
        }
        else {
            TreeAdjustFactor[band] = 1.;
        }
        if (TreeAdjustFactor[band] != 1) {
            log_warn("Tree adjust factor for band %zu is equal to %f.",
                     band, TreeAdjustFactor[band]);
        }
    }

    cv_baresoil = 0;
    cv_veg = 0;
    cv_overstory = 0;
    cv_snow = 0;

    // Initialize output data to zero
    zero_output_list(out_data);

    // Set output versions of input forcings
    out_data[OUT_AIR_TEMP][0] = force->air_temp[NR];
    out_data[OUT_DENSITY][0] = force->density[NR];
    out_data[OUT_LWDOWN][0] = force->longwave[NR];
    out_data[OUT_PREC][0] = force->out_prec;  // mm over grid cell
    out_data[OUT_PRESSURE][0] = force->pressure[NR] / PA_PER_KPA;
    out_data[OUT_QAIR][0] = CONST_EPS * force->vp[NR] /
                            force->pressure[NR];
    out_data[OUT_RAINF][0] = force->out_rain;   // mm over grid cell
    out_data[OUT_REL_HUMID][0] = FRACT_TO_PERCENT * force->vp[NR] /
                                 (force->vp[NR] + force->vpd[NR]);
    if (options.LAKES && lake_con->Cl[0] > 0) {
        out_data[OUT_LAKE_CHAN_IN][0] = force->channel_in[NR];  // mm over grid cell
    }
    else {
        out_data[OUT_LAKE_CHAN_IN][0] = 0;
    }
    out_data[OUT_SWDOWN][0] = force->shortwave[NR];
    out_data[OUT_SNOWF][0] = force->out_snow;   // mm over grid cell
    out_data[OUT_VP][0] = force->vp[NR] / PA_PER_KPA;
    out_data[OUT_VPD][0] = force->vpd[NR] / PA_PER_KPA;
    out_data[OUT_WIND][0] = force->wind[NR];
    if (options.CARBON) {
        out_data[OUT_CATM][0] = force->Catm[NR] / PPM_to_MIXRATIO;
        out_data[OUT_FDIR][0] = force->fdir[NR];
        out_data[OUT_PAR][0] = force->par[NR];
    }
    else {
        out_data[OUT_CATM][0] = MISSING;
        out_data[OUT_FDIR][0] = MISSING;
        out_data[OUT_PAR][0] = MISSING;
    }

    /****************************************
       Store Output for all Vegetation Types (except lakes)
    ****************************************/
    for (veg = 0; veg <= veg_con[0].vegetat_type_num; veg++) {
        Cv = veg_con[veg].Cv;
        Clake = 0;
        Nbands = options.SNOW_BAND;
        IsWet = false;

        if (veg < veg_con[0].vegetat_type_num) {
            HasVeg = true;
        }
        else {
            HasVeg = false;
        }

        if (Cv > 0) {
            // Check if this is lake/wetland tile
            if (options.LAKES && veg_con[veg].LAKE) {
                Clake = lake_var.sarea / lake_con->basin[0];
                Nbands = 1;
                IsWet = true;
            }

            overstory = veg_lib[veg_con[veg].veg_class].overstory;

            /*********************************
               Store Output for all Bands
            *********************************/
            for (band = 0; band < Nbands; band++) {
                ThisAreaFract = AreaFract[band];
                ThisTreeAdjust = TreeAdjustFactor[band];
                if (IsWet) {
                    ThisAreaFract = 1;
                    ThisTreeAdjust = 1;
                }

                if (ThisAreaFract > 0. && (veg == veg_con[0].vegetat_type_num ||
                                           (!AboveTreeLine[band] ||
                                            (AboveTreeLine[band] &&
                                             !overstory)))) {
                    /** compute running totals of various landcovers **/
                    if (HasVeg) {
                        cv_veg += Cv * ThisAreaFract * ThisTreeAdjust;
                    }
                    else {
                        cv_baresoil += Cv * ThisAreaFract * ThisTreeAdjust;
                    }
                    if (overstory) {
                        cv_overstory += Cv * ThisAreaFract * ThisTreeAdjust;
                    }
                    if (snow[veg][band].swq > 0.0) {
                        cv_snow += Cv * ThisAreaFract * ThisTreeAdjust;
                    }

                    /*********************************
                       Record Water Balance Terms
                    *********************************/
                    collect_wb_terms(cell[veg][band],
                                     veg_var[veg][band],
                                     snow[veg][band],
                                     Cv,
                                     ThisAreaFract,
                                     ThisTreeAdjust,
                                     HasVeg,
                                     (1 - Clake),
                                     overstory,
                                     frost_fract,
                                     out_data);

                    /**********************************
                       Record Energy Balance Terms
                    **********************************/
                    collect_eb_terms(energy[veg][band],
                                     snow[veg][band],
                                     cell[veg][band],
                                     Cv,
                                     ThisAreaFract,
                                     ThisTreeAdjust,
                                     HasVeg,
                                     0,
                                     (1 - Clake),
                                     overstory,
                                     band,
                                     frost_fract,
                                     frost_slope,
                                     out_data);

                    // Store Wetland-Specific Variables
                    if (IsWet) {
                        // Wetland soil temperatures
                        for (i = 0; i < options.Nnode; i++) {
                            out_data[OUT_SOIL_TNODE_WL][i] =
                                energy[veg][band].T[i];
                        }
                    }

                    /**********************************
                       Record Lake Variables
                    **********************************/
                    if (IsWet) {
                        // Override some variables of soil under lake with those of wetland
                        // This is for those variables whose lake values shouldn't be included
                        // in grid cell average
                        // Note: doing this for eb terms will lead to reporting of eb errors
                        // this should be fixed when we implement full thermal solution beneath lake
                        for (i = 0; i < MAX_FRONTS; i++) {
                            lake_var.energy.fdepth[i] =
                                energy[veg][band].fdepth[i];
                            lake_var.energy.tdepth[i] =
                                energy[veg][band].fdepth[i];
                        }
                        for (i = 0; i < options.Nnode; i++) {
                            lake_var.energy.ice[i] = energy[veg][band].ice[i];
                            lake_var.energy.T[i] = energy[veg][band].T[i];
                        }
                        lake_var.soil.pot_evap =
                            cell[veg][band].pot_evap;
                        lake_var.soil.rootmoist = cell[veg][band].rootmoist;
                        lake_var.energy.deltaH = energy[veg][band].deltaH;
                        lake_var.energy.fusion = energy[veg][band].fusion;
                        lake_var.energy.grnd_flux = energy[veg][band].grnd_flux;


                        /*********************************
                           Record Water Balance Terms
                        *********************************/
                        collect_wb_terms(lake_var.soil,
                                         veg_var[0][0],
                                         lake_var.snow,
                                         Cv,
                                         ThisAreaFract,
                                         ThisTreeAdjust,
                                         0,
                                         Clake,
                                         overstory,
                                         frost_fract,
                                         out_data);

                        /**********************************
                           Record Energy Balance Terms
                        **********************************/
                        collect_eb_terms(lake_var.energy,
                                         lake_var.snow,
                                         lake_var.soil,
                                         Cv,
                                         ThisAreaFract,
                                         ThisTreeAdjust,
                                         0,
                                         1,
                                         Clake,
                                         overstory,
                                         band,
                                         frost_fract,
                                         frost_slope,
                                         out_data);

                        // Store Lake-Specific Variables

                        // Lake ice
                        if (lake_var.new_ice_area > 0.0) {
                            out_data[OUT_LAKE_ICE][0] =
                                (lake_var.ice_water_eq /
                                 lake_var.new_ice_area) * CONST_RHOICE /
                                CONST_RHOFW;
                            out_data[OUT_LAKE_ICE_TEMP][0] =
                                lake_var.tempi;
                            out_data[OUT_LAKE_ICE_HEIGHT][0] =
                                lake_var.hice;
                            out_data[OUT_LAKE_SWE][0] = lake_var.swe /
                                                        lake_var.areai;       // m over lake ice
                            out_data[OUT_LAKE_SWE_V][0] = lake_var.swe;  // m3
                        }
                        else {
                            out_data[OUT_LAKE_ICE][0] = 0.0;
                            out_data[OUT_LAKE_ICE_TEMP][0] = 0.0;
                            out_data[OUT_LAKE_ICE_HEIGHT][0] = 0.0;
                            out_data[OUT_LAKE_SWE][0] = 0.0;
                            out_data[OUT_LAKE_SWE_V][0] = 0.0;
                        }
                        out_data[OUT_LAKE_DSWE_V][0] = lake_var.swe -
                                                       lake_var.swe_save;       // m3
                        // same as OUT_LAKE_MOIST
                        out_data[OUT_LAKE_DSWE][0] =
                            (lake_var.swe - lake_var.swe_save) * MM_PER_M /
                            soil_con->cell_area;

                        // Lake dimensions
                        out_data[OUT_LAKE_AREA_FRAC][0] = Cv * Clake;
                        out_data[OUT_LAKE_DEPTH][0] = lake_var.ldepth;
                        out_data[OUT_LAKE_SURF_AREA][0] = lake_var.sarea;
                        if (out_data[OUT_LAKE_SURF_AREA][0] > 0) {
                            out_data[OUT_LAKE_ICE_FRACT][0] =
                                lake_var.new_ice_area /
                                out_data[OUT_LAKE_SURF_AREA][0];
                        }
                        else {
                            out_data[OUT_LAKE_ICE_FRACT][0] = 0.;
                        }
                        out_data[OUT_LAKE_VOLUME][0] = lake_var.volume;
                        out_data[OUT_LAKE_DSTOR_V][0] = lake_var.volume -
                                                        lake_var.
                                                        volume_save;
                        // mm over gridcell
                        out_data[OUT_LAKE_DSTOR][0] =
                            (lake_var.volume - lake_var.volume_save) *
                            MM_PER_M /
                            soil_con->cell_area;

                        // Other lake characteristics
                        out_data[OUT_LAKE_SURF_TEMP][0] = lake_var.temp[0];
                        if (out_data[OUT_LAKE_SURF_AREA][0] > 0) {
                            // mm over gridcell
                            out_data[OUT_LAKE_MOIST][0] =
                                (lake_var.volume / soil_con->cell_area) *
                                MM_PER_M;
                            // same as OUT_LAKE_MOIST
                            out_data[OUT_SURFSTOR][0] =
                                (lake_var.volume / soil_con->cell_area) *
                                MM_PER_M;
                        }
                        else {
                            out_data[OUT_LAKE_MOIST][0] = 0;
                            out_data[OUT_SURFSTOR][0] = 0;
                        }

                        // Lake moisture fluxes
                        out_data[OUT_LAKE_BF_IN_V][0] =
                            lake_var.baseflow_in;  // m3
                        out_data[OUT_LAKE_BF_OUT_V][0] =
                            lake_var.baseflow_out;  // m3
                        out_data[OUT_LAKE_CHAN_IN_V][0] =
                            lake_var.channel_in;  // m3
                        out_data[OUT_LAKE_CHAN_OUT_V][0] =
                            lake_var.runoff_out;  // m3
                        out_data[OUT_LAKE_EVAP_V][0] = lake_var.evapw;  // m3
                        out_data[OUT_LAKE_PREC_V][0] = lake_var.prec;  // m3
                        out_data[OUT_LAKE_RCHRG_V][0] = lake_var.recharge;  // m3
                        out_data[OUT_LAKE_RO_IN_V][0] = lake_var.runoff_in;  // m3
                        out_data[OUT_LAKE_VAPFLX_V][0] =
                            lake_var.vapor_flux;  // m3
                        out_data[OUT_LAKE_BF_IN][0] =
                            lake_var.baseflow_in * MM_PER_M /
                            soil_con->cell_area;  // mm over gridcell
                        out_data[OUT_LAKE_BF_OUT][0] =
                            lake_var.baseflow_out * MM_PER_M /
                            soil_con->cell_area;  // mm over gridcell
                        out_data[OUT_LAKE_CHAN_OUT][0] =
                            lake_var.runoff_out * MM_PER_M /
                            soil_con->cell_area;  // mm over gridcell
                        // mm over gridcell
                        out_data[OUT_LAKE_EVAP][0] = lake_var.evapw * MM_PER_M /
                                                     soil_con->cell_area;
                        // mm over gridcell
                        out_data[OUT_LAKE_RCHRG][0] = lake_var.recharge *
                                                      MM_PER_M /
                                                      soil_con->cell_area;
                        // mm over gridcell
                        out_data[OUT_LAKE_RO_IN][0] = lake_var.runoff_in *
                                                      MM_PER_M /
                                                      soil_con->cell_area;
                        out_data[OUT_LAKE_VAPFLX][0] =
                            lake_var.vapor_flux * MM_PER_M /
                            soil_con->cell_area;  // mm over gridcell
                    } // End if options.LAKES etc.
                } // End if ThisAreaFract etc.
            } // End loop over bands
        } // End if Cv > 0
    } // End loop over veg


    /*****************************************
       Finish aggregation of special-case variables
    *****************************************/
    // Normalize quantities that aren't present over entire grid cell
    if (cv_baresoil > 0) {
        out_data[OUT_BARESOILT][0] /= cv_baresoil;
    }
    if (cv_veg > 0) {
        out_data[OUT_VEGT][0] /= cv_veg;
    }
    if (cv_overstory > 0) {
        out_data[OUT_AERO_COND2][0] /= cv_overstory;
    }
    if (cv_snow > 0) {
        out_data[OUT_SALBEDO][0] /= cv_snow;
        out_data[OUT_SNOW_SURF_TEMP][0] /= cv_snow;
        out_data[OUT_SNOW_PACK_TEMP][0] /= cv_snow;
    }

    // Radiative temperature
    out_data[OUT_RAD_TEMP][0] = pow(out_data[OUT_RAD_TEMP][0], 0.25);

    // Aerodynamic conductance and resistance
    if (out_data[OUT_AERO_COND1][0] > DBL_EPSILON) {
        out_data[OUT_AERO_RESIST1][0] = 1 / out_data[OUT_AERO_COND1][0];
    }
    else {
        out_data[OUT_AERO_RESIST1][0] = param.HUGE_RESIST;
    }
    if (out_data[OUT_AERO_COND2][0] > DBL_EPSILON) {
        out_data[OUT_AERO_RESIST2][0] = 1 / out_data[OUT_AERO_COND2][0];
    }
    else {
        out_data[OUT_AERO_RESIST2][0] = param.HUGE_RESIST;
    }
    if (out_data[OUT_AERO_COND][0] > DBL_EPSILON) {
        out_data[OUT_AERO_RESIST][0] = 1 / out_data[OUT_AERO_COND][0];
    }
    else {
        out_data[OUT_AERO_RESIST][0] = param.HUGE_RESIST;
    }

    /*****************************************
       Compute derived variables
    *****************************************/
    // Water balance terms
    out_data[OUT_DELSOILMOIST][0] = 0;
    for (index = 0; index < options.Nlayer; index++) {
        out_data[OUT_SOIL_MOIST][index] =
            out_data[OUT_SOIL_LIQ][index] +
            out_data[OUT_SOIL_ICE][index];
        out_data[OUT_DELSOILMOIST][0] +=
            out_data[OUT_SOIL_MOIST][index];

        out_data[OUT_SMLIQFRAC][index] = out_data[OUT_SOIL_LIQ][index] /
                                         out_data[OUT_SOIL_MOIST][index];
        out_data[OUT_SMFROZFRAC][index] = 1 - out_data[OUT_SMLIQFRAC][index];
        out_data[OUT_SOIL_LIQ_FRAC][index] = out_data[OUT_SOIL_LIQ][index] /
                                             (depth[index] *
                                              MM_PER_M);
        out_data[OUT_SOIL_ICE_FRAC][index] = out_data[OUT_SOIL_ICE][index] /
                                             (depth[index] *
                                              MM_PER_M);
    }
    out_data[OUT_DELSOILMOIST][0] -= save_data->total_soil_moist;
    out_data[OUT_DELSWE][0] = out_data[OUT_SWE][0] +
                              out_data[OUT_SNOW_CANOPY][0] -
                              save_data->swe;
    out_data[OUT_DELINTERCEPT][0] = out_data[OUT_WDEW][0] - save_data->wdew;
    out_data[OUT_DELSURFSTOR][0] = out_data[OUT_SURFSTOR][0] -
                                   save_data->surfstor;

    // Energy terms
    out_data[OUT_REFREEZE][0] =
        (out_data[OUT_RFRZ_ENERGY][0] / CONST_LATICE) * dt_sec;
    out_data[OUT_R_NET][0] = out_data[OUT_SWNET][0] + out_data[OUT_LWNET][0];

    // Save current moisture state for use in next time step
    save_data->total_soil_moist = 0;
    for (index = 0; index < options.Nlayer; index++) {
        save_data->total_soil_moist += out_data[OUT_SOIL_MOIST][index];
    }
    save_data->surfstor = out_data[OUT_SURFSTOR][0];
    save_data->swe = out_data[OUT_SWE][0] + out_data[OUT_SNOW_CANOPY][0];
    save_data->wdew = out_data[OUT_WDEW][0];

    // Carbon Terms
    if (options.CARBON) {
        out_data[OUT_RHET][0] *= dt_sec / SEC_PER_DAY;  // convert to gC/m2d
        out_data[OUT_NEE][0] = out_data[OUT_NPP][0] - out_data[OUT_RHET][0];
    }

    /********************
       Check Water Balance
    ********************/
    inflow = out_data[OUT_PREC][0] + out_data[OUT_LAKE_CHAN_IN][0];  // mm over grid cell
    outflow = out_data[OUT_EVAP][0] + out_data[OUT_RUNOFF][0] +
              out_data[OUT_BASEFLOW][0];  // mm over grid cell
    storage = 0.;
    for (index = 0; index < options.Nlayer; index++) {
        storage += out_data[OUT_SOIL_LIQ][index] +
                   out_data[OUT_SOIL_ICE][index];
    }
    storage += out_data[OUT_SWE][0] + out_data[OUT_SNOW_CANOPY][0] +
               out_data[OUT_WDEW][0] + out_data[OUT_SURFSTOR][0];
    out_data[OUT_WATER_ERROR][0] = calc_water_balance_error(inflow,
                                                            outflow,
                                                            storage,
                                                            save_data->total_moist_storage);

    // Store total storage for next timestep
    save_data->total_moist_storage = storage;

    /********************
       Check Energy Balance
    ********************/
    if (options.FULL_ENERGY) {
        out_data[OUT_ENERGY_ERROR][0] = \
            calc_energy_balance_error(out_data[OUT_SWNET][0] +
                                      out_data[OUT_LWNET][0],
                                      out_data[OUT_LATENT][0] +
                                      out_data[OUT_LATENT_SUB][0],
                                      out_data[OUT_SENSIBLE][0] +
                                      out_data[OUT_ADV_SENS][0],
                                      out_data[OUT_GRND_FLUX][0] +
                                      out_data[OUT_DELTAH][0] +
                                      out_data[OUT_FUSION][0],
                                      out_data[OUT_ADVECTION][0] -
                                      out_data[OUT_DELTACC][0] +
                                      out_data[OUT_SNOW_FLUX][0] +
                                      out_data[OUT_RFRZ_ENERGY][0]);
    }
    else {
        out_data[OUT_ENERGY_ERROR][0] = MISSING;
    }

    free((char *) (TreeAdjustFactor));

    // vic_run run time
    out_data[OUT_TIME_VICRUN_WALL][0] = timer->delta_wall;
    out_data[OUT_TIME_VICRUN_CPU][0] = timer->delta_cpu;
}

/******************************************************************************
 * @brief    This routine collects water balance terms.
 *****************************************************************************/
void
collect_wb_terms(cell_data_struct cell,
                 veg_var_struct   veg_var,
                 snow_data_struct snow,
                 double           Cv,
                 double           AreaFract,
                 double           TreeAdjustFactor,
                 bool             HasVeg,
                 double           lakefactor,
                 bool             overstory,
                 double          *frost_fract,
                 double         **out_data)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   AreaFactor;
    double                   tmp_evap;
    double                   tmp_cond1;
    double                   tmp_cond2;
    double                   tmp_moist;
    double                   tmp_ice;
    size_t                   index;
    size_t                   frost_area;

    AreaFactor = Cv * AreaFract * TreeAdjustFactor * lakefactor;

    /** record evaporation components **/
    tmp_evap = 0.0;
    for (index = 0; index < options.Nlayer; index++) {
        tmp_evap += cell.layer[index].evap;
        if (HasVeg) {
            out_data[OUT_EVAP_BARE][0] += cell.layer[index].evap *
                                          cell.layer[index].bare_evap_frac
                                          *
                                          AreaFactor;
            out_data[OUT_TRANSP_VEG][0] += cell.layer[index].evap *
                                           (1 -
                                            cell.layer[index].
                                            bare_evap_frac) * AreaFactor;
        }
        else {
            out_data[OUT_EVAP_BARE][0] += cell.layer[index].evap *
                                          AreaFactor;
        }
    }
    tmp_evap += snow.vapor_flux * MM_PER_M;
    out_data[OUT_SUB_SNOW][0] += snow.vapor_flux * MM_PER_M * AreaFactor;
    out_data[OUT_SUB_SURFACE][0] += snow.surface_flux * MM_PER_M *
                                    AreaFactor;
    out_data[OUT_SUB_BLOWING][0] += snow.blowing_flux * MM_PER_M *
                                    AreaFactor;
    if (HasVeg) {
        tmp_evap += snow.canopy_vapor_flux * MM_PER_M;
        out_data[OUT_SUB_CANOP][0] += snow.canopy_vapor_flux * MM_PER_M *
                                      AreaFactor;
    }
    if (HasVeg) {
        tmp_evap += veg_var.canopyevap;
        out_data[OUT_EVAP_CANOP][0] += veg_var.canopyevap * AreaFactor;
    }
    out_data[OUT_EVAP][0] += tmp_evap * AreaFactor;  // mm over gridcell

    /** record potential evap **/
    out_data[OUT_PET][0] += cell.pot_evap * AreaFactor;

    /** record saturated area fraction **/
    out_data[OUT_ASAT][0] += cell.asat * AreaFactor;

    /** record runoff **/
    out_data[OUT_RUNOFF][0] += cell.runoff * AreaFactor;

    /** record baseflow **/
    out_data[OUT_BASEFLOW][0] += cell.baseflow * AreaFactor;

    /** record inflow **/
    out_data[OUT_INFLOW][0] += (cell.inflow) * AreaFactor;

    /** record canopy interception **/
    if (HasVeg) {
        out_data[OUT_WDEW][0] += veg_var.Wdew * AreaFactor;
    }

    /** record LAI **/
    out_data[OUT_LAI][0] += veg_var.LAI * AreaFactor;

    /** record fcanopy **/
    out_data[OUT_FCANOPY][0] += veg_var.fcanopy * AreaFactor;

    /** record aerodynamic conductance and resistance **/
    if (cell.aero_resist[0] > DBL_EPSILON) {
        tmp_cond1 = (1 / cell.aero_resist[0]) * AreaFactor;
    }
    else {
        tmp_cond1 = param.HUGE_RESIST;
    }
    out_data[OUT_AERO_COND1][0] += tmp_cond1;
    if (overstory) {
        if (cell.aero_resist[1] > DBL_EPSILON) {
            tmp_cond2 = (1 / cell.aero_resist[1]) * AreaFactor;
        }
        else {
            tmp_cond2 = param.HUGE_RESIST;
        }
        out_data[OUT_AERO_COND2][0] += tmp_cond2;
        out_data[OUT_AERO_COND][0] += tmp_cond2;
    }
    else {
        out_data[OUT_AERO_COND][0] += tmp_cond1;
    }

    /** record layer moistures **/
    for (index = 0; index < options.Nlayer; index++) {
        tmp_moist = cell.layer[index].moist;
        tmp_ice = 0;
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            tmp_ice +=
                (cell.layer[index].ice[frost_area] * frost_fract[frost_area]);
        }
        tmp_moist -= tmp_ice;

        out_data[OUT_SOIL_LIQ][index] += tmp_moist * AreaFactor;
        out_data[OUT_SOIL_ICE][index] += tmp_ice * AreaFactor;
    }
    out_data[OUT_SOIL_WET][0] += cell.wetness * AreaFactor;
    out_data[OUT_ROOTMOIST][0] += cell.rootmoist * AreaFactor;

    /** record water table position **/
    out_data[OUT_ZWT][0] += cell.zwt * AreaFactor;
    out_data[OUT_ZWT_LUMPED][0] += cell.zwt_lumped * AreaFactor;

    /** record layer temperatures **/
    for (index = 0; index < options.Nlayer; index++) {
        out_data[OUT_SOIL_TEMP][index] += cell.layer[index].T * AreaFactor;
    }

    /*****************************
       Record Snow Pack Variables
    *****************************/

    /** record snow water equivalence **/
    out_data[OUT_SWE][0] += snow.swq * AreaFactor * MM_PER_M;

    /** record snowpack depth **/
    out_data[OUT_SNOW_DEPTH][0] += snow.depth * AreaFactor * CM_PER_M;

    /** record snowpack albedo, temperature **/
    if (snow.swq > 0.0) {
        out_data[OUT_SALBEDO][0] += snow.albedo * AreaFactor;
        out_data[OUT_SNOW_SURF_TEMP][0] += snow.surf_temp * AreaFactor;
        out_data[OUT_SNOW_PACK_TEMP][0] += snow.pack_temp * AreaFactor;
    }

    /** record canopy intercepted snow **/
    if (HasVeg) {
        out_data[OUT_SNOW_CANOPY][0] += (snow.snow_canopy) * AreaFactor *
                                        MM_PER_M;
    }

    /** record snowpack melt **/
    out_data[OUT_SNOW_MELT][0] += snow.melt * AreaFactor;

    /** record snow cover fraction **/
    out_data[OUT_SNOW_COVER][0] += snow.coverage * AreaFactor;

    /*****************************
       Record Carbon Cycling Variables
    *****************************/
    if (options.CARBON) {
        out_data[OUT_APAR][0] += veg_var.aPAR * AreaFactor;
        out_data[OUT_GPP][0] += veg_var.GPP * CONST_MWC / MOLE_PER_KMOLE *
                                CONST_CDAY *
                                AreaFactor;
        out_data[OUT_RAUT][0] += veg_var.Raut * CONST_MWC /
                                 MOLE_PER_KMOLE * CONST_CDAY *
                                 AreaFactor;
        out_data[OUT_NPP][0] += veg_var.NPP * CONST_MWC / MOLE_PER_KMOLE *
                                CONST_CDAY *
                                AreaFactor;
        out_data[OUT_LITTERFALL][0] += veg_var.Litterfall * AreaFactor;
        out_data[OUT_RHET][0] += cell.RhTot * AreaFactor;
        out_data[OUT_CLITTER][0] += cell.CLitter * AreaFactor;
        out_data[OUT_CINTER][0] += cell.CInter * AreaFactor;
        out_data[OUT_CSLOW][0] += cell.CSlow * AreaFactor;
    }
}

/******************************************************************************
 * @brief    This routine collects energy balance terms.
 *****************************************************************************/
void
collect_eb_terms(energy_bal_struct energy,
                 snow_data_struct  snow,
                 cell_data_struct  cell_wet,
                 double            Cv,
                 double            AreaFract,
                 double            TreeAdjustFactor,
                 bool              HasVeg,
                 bool              IsWet,
                 double            lakefactor,
                 bool              overstory,
                 int               band,
                 double           *frost_fract,
                 double            frost_slope,
                 double          **out_data)
{
    extern option_struct options;
    double               AreaFactor;
    double               tmp_fract;
    double               rad_temp;
    double               surf_temp;
    size_t               index;
    size_t               frost_area;

    AreaFactor = Cv * AreaFract * TreeAdjustFactor * lakefactor;

    /**********************************
       Record Frozen Soil Variables
    **********************************/

    /** record freezing and thawing front depths **/
    if (options.FROZEN_SOIL) {
        for (index = 0; index < MAX_FRONTS; index++) {
            if (energy.fdepth[index] != MISSING) {
                out_data[OUT_FDEPTH][index] += energy.fdepth[index] *
                                               AreaFactor * CM_PER_M;
            }
            if (energy.tdepth[index] != MISSING) {
                out_data[OUT_TDEPTH][index] += energy.tdepth[index] *
                                               AreaFactor * CM_PER_M;
            }
        }
    }

    tmp_fract = 0;
    for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
        if (cell_wet.layer[0].ice[frost_area]) {
            tmp_fract += frost_fract[frost_area];
        }
    }
    out_data[OUT_SURF_FROST_FRAC][0] += tmp_fract * AreaFactor;

    tmp_fract = 0;
    if ((energy.T[0] + frost_slope / 2.) > 0) {
        if ((energy.T[0] - frost_slope / 2.) <= 0) {
            tmp_fract +=
                linear_interp(0, (energy.T[0] + frost_slope / 2.),
                              (energy.T[0] - frost_slope / 2.), 1,
                              0) * AreaFactor;
        }
    }
    else {
        tmp_fract += 1 * AreaFactor;
    }

    /**********************************
       Record Energy Balance Variables
    **********************************/

    /** record surface radiative temperature **/
    if (overstory && snow.snow && !(options.LAKES && IsWet)) {
        rad_temp = energy.Tfoliage + CONST_TKFRZ;
    }
    else {
        rad_temp = energy.Tsurf + CONST_TKFRZ;
    }

    /** record surface skin temperature **/
    surf_temp = energy.Tsurf;

    /** record landcover temperature **/
    if (!HasVeg) {
        // landcover is bare soil
        out_data[OUT_BARESOILT][0] +=
            (rad_temp - CONST_TKFRZ) * AreaFactor;
    }
    else {
        // landcover is vegetation
        if (overstory && !snow.snow) {
            // here, rad_temp will be wrong since it will pick the understory temperature
            out_data[OUT_VEGT][0] += energy.Tfoliage * AreaFactor;
        }
        else {
            out_data[OUT_VEGT][0] += (rad_temp - CONST_TKFRZ) * AreaFactor;
        }
    }

    /** record mean surface temperature [C]  **/
    out_data[OUT_SURF_TEMP][0] += surf_temp * AreaFactor;

    /** record thermal node temperatures **/
    for (index = 0; index < options.Nnode; index++) {
        out_data[OUT_SOIL_TNODE][index] += energy.T[index] * AreaFactor;
    }
    if (IsWet) {
        for (index = 0; index < options.Nnode; index++) {
            out_data[OUT_SOIL_TNODE_WL][index] = energy.T[index];
        }
    }

    /** record temperature flags  **/
    out_data[OUT_SURFT_FBFLAG][0] += energy.Tsurf_fbflag * AreaFactor;
    for (index = 0; index < options.Nnode; index++) {
        out_data[OUT_SOILT_FBFLAG][index] += energy.T_fbflag[index] *
                                             AreaFactor;
    }
    out_data[OUT_SNOWT_FBFLAG][0] += snow.surf_temp_fbflag * AreaFactor;
    out_data[OUT_TFOL_FBFLAG][0] += energy.Tfoliage_fbflag * AreaFactor;
    out_data[OUT_TCAN_FBFLAG][0] += energy.Tcanopy_fbflag * AreaFactor;

    /** record net shortwave radiation **/
    out_data[OUT_SWNET][0] += energy.NetShortAtmos * AreaFactor;

    /** record net longwave radiation **/
    out_data[OUT_LWNET][0] += energy.NetLongAtmos * AreaFactor;

    /** record incoming longwave radiation at ground surface (under veg) **/
    if (snow.snow && overstory) {
        out_data[OUT_IN_LONG][0] += energy.LongOverIn * AreaFactor;
    }
    else {
        out_data[OUT_IN_LONG][0] += energy.LongUnderIn * AreaFactor;
    }

    /** record albedo **/
    if (snow.snow && overstory) {
        out_data[OUT_ALBEDO][0] += energy.AlbedoOver * AreaFactor;
    }
    else {
        out_data[OUT_ALBEDO][0] += energy.AlbedoUnder * AreaFactor;
    }

    /** record latent heat flux **/
    out_data[OUT_LATENT][0] -= energy.AtmosLatent * AreaFactor;

    /** record latent heat flux from sublimation **/
    out_data[OUT_LATENT_SUB][0] -= energy.AtmosLatentSub * AreaFactor;

    /** record sensible heat flux **/
    out_data[OUT_SENSIBLE][0] -= energy.AtmosSensible * AreaFactor;

    /** record ground heat flux (+ heat storage) **/
    out_data[OUT_GRND_FLUX][0] -= energy.grnd_flux * AreaFactor;

    /** record heat storage **/
    out_data[OUT_DELTAH][0] -= energy.deltaH * AreaFactor;

    /** record heat of fusion **/
    out_data[OUT_FUSION][0] -= energy.fusion * AreaFactor;

    /** record radiative effective temperature [K],
        emissivities set = 1.0  **/
    out_data[OUT_RAD_TEMP][0] +=
        ((rad_temp) * (rad_temp) * (rad_temp) * (rad_temp)) * AreaFactor;

    /** record snowpack cold content **/
    out_data[OUT_DELTACC][0] += energy.deltaCC * AreaFactor;

    /** record snowpack advection **/
    if (snow.snow && overstory) {
        out_data[OUT_ADVECTION][0] += energy.canopy_advection * AreaFactor;
    }
    out_data[OUT_ADVECTION][0] += energy.advection * AreaFactor;

    /** record snow energy flux **/
    out_data[OUT_SNOW_FLUX][0] += energy.snow_flux * AreaFactor;

    /** record refreeze energy **/
    if (snow.snow && overstory) {
        out_data[OUT_RFRZ_ENERGY][0] += energy.canopy_refreeze *
                                        AreaFactor;
    }
    out_data[OUT_RFRZ_ENERGY][0] += energy.refreeze_energy * AreaFactor;

    /** record melt energy **/
    out_data[OUT_MELT_ENERGY][0] += energy.melt_energy * AreaFactor;

    /** record advected sensible heat energy **/
    if (!overstory) {
        out_data[OUT_ADV_SENS][0] -= energy.advected_sensible * AreaFactor;
    }

    /**********************************
       Record Band-Specific Variables
    **********************************/

    /** record band snow water equivalent **/
    out_data[OUT_SWE_BAND][band] += snow.swq * Cv * lakefactor * MM_PER_M;

    /** record band snowpack depth **/
    out_data[OUT_SNOW_DEPTH_BAND][band] += snow.depth * Cv * lakefactor *
                                           CM_PER_M;

    /** record band canopy intercepted snow **/
    if (HasVeg) {
        out_data[OUT_SNOW_CANOPY_BAND][band] += (snow.snow_canopy) * Cv *
                                                lakefactor * MM_PER_M;
    }

    /** record band snow melt **/
    out_data[OUT_SNOW_MELT_BAND][band] += snow.melt * Cv * lakefactor;

    /** record band snow coverage **/
    out_data[OUT_SNOW_COVER_BAND][band] += snow.coverage * Cv * lakefactor;

    /** record band cold content **/
    out_data[OUT_DELTACC_BAND][band] += energy.deltaCC * Cv * lakefactor;

    /** record band advection **/
    out_data[OUT_ADVECTION_BAND][band] += energy.advection * Cv *
                                          lakefactor;

    /** record band snow flux **/
    out_data[OUT_SNOW_FLUX_BAND][band] += energy.snow_flux * Cv *
                                          lakefactor;

    /** record band refreeze energy **/
    out_data[OUT_RFRZ_ENERGY_BAND][band] += energy.refreeze_energy * Cv *
                                            lakefactor;

    /** record band melt energy **/
    out_data[OUT_MELT_ENERGY_BAND][band] += energy.melt_energy * Cv *
                                            lakefactor;

    /** record band advected sensble heat **/
    out_data[OUT_ADV_SENS_BAND][band] -= energy.advected_sensible * Cv *
                                         lakefactor;

    /** record surface layer temperature **/
    out_data[OUT_SNOW_SURFT_BAND][band] += snow.surf_temp * Cv *
                                           lakefactor;

    /** record pack layer temperature **/
    out_data[OUT_SNOW_PACKT_BAND][band] += snow.pack_temp * Cv *
                                           lakefactor;

    /** record latent heat of sublimation **/
    out_data[OUT_LATENT_SUB_BAND][band] += energy.latent_sub * Cv *
                                           lakefactor;

    /** record band net downwards shortwave radiation **/
    out_data[OUT_SWNET_BAND][band] += energy.NetShortAtmos * Cv *
                                      lakefactor;

    /** record band net downwards longwave radiation **/
    out_data[OUT_LWNET_BAND][band] += energy.NetLongAtmos * Cv *
                                      lakefactor;

    /** record band albedo **/
    if (snow.snow && overstory) {
        out_data[OUT_ALBEDO_BAND][band] += energy.AlbedoOver * Cv *
                                           lakefactor;
    }
    else {
        out_data[OUT_ALBEDO_BAND][band] += energy.AlbedoUnder * Cv *
                                           lakefactor;
    }

    /** record band net latent heat flux **/
    out_data[OUT_LATENT_BAND][band] -= energy.latent * Cv * lakefactor;

    /** record band net sensible heat flux **/
    out_data[OUT_SENSIBLE_BAND][band] -= energy.sensible * Cv * lakefactor;

    /** record band net ground heat flux **/
    out_data[OUT_GRND_FLUX_BAND][band] -= energy.grnd_flux * Cv *
                                          lakefactor;
}

/******************************************************************************
 * @brief    Initialize the save data structure.
 *****************************************************************************/
void
initialize_save_data(all_vars_struct   *all_vars,
                     force_data_struct *atmos,
                     soil_con_struct   *soil_con,
                     veg_con_struct    *veg_con,
                     veg_lib_struct    *veg_lib,
                     lake_con_struct   *lake_con,
                     double           **out_data,
                     save_data_struct  *save_data,
                     timer_struct      *timer)
{
    // Calling put data will populate the save data storage terms
    put_data(all_vars, atmos, soil_con, veg_con, veg_lib, lake_con,
             out_data, save_data, timer);

    zero_output_list(out_data);
}

/******************************************************************************
 * @brief    This subroutine computes the overall model water balance, and
 *           warns the model user if large errors are found.
 *****************************************************************************/
double
calc_water_balance_error(double inflow,
                         double outflow,
                         double storage,
                         double last_storage)
{
    double error;

    error = inflow - outflow - (storage - last_storage);

    return(error);
}

/******************************************************************************
 * @brief    This subroutine computes the overall model energy balance.
 *****************************************************************************/
double
calc_energy_balance_error(double net_rad,
                          double latent,
                          double sensible,
                          double grnd_flux,
                          double snow_fluxes)
{
    double error;

    error = net_rad - latent - sensible - grnd_flux + snow_fluxes;

    return(error);
}
