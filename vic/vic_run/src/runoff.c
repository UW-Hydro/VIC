/******************************************************************************
* @section DESCRIPTION
*
* Calculate infiltration and runoff from the surface, gravity driven drainage
* between all soil layers, and generates baseflow from the bottom layer.
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
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Calculate infiltration and runoff from the surface, gravity driven
*           drainage between all soil layers, and generates baseflow from the
*           bottom layer.
******************************************************************************/
int
runoff(cell_data_struct  *cell,
       energy_bal_struct *energy,
       soil_con_struct   *soil_con,
       double             ppt,
       double            *frost_fract,
       int                Nnodes)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    size_t                     lindex;
    size_t                     time_step;
    int                        last_index;
    int                        tmplayer;
    int                        fidx;
    int                        ErrorFlag;
    double                     A, frac;
    double                     tmp_runoff;
    double                     inflow;
    double                     resid_moist[MAX_LAYERS]; // residual moisture (mm)
    double                     org_moist[MAX_LAYERS]; // total soil moisture (liquid and frozen) at beginning of this function (mm)
    double                     avail_liq[MAX_LAYERS][MAX_FROST_AREAS]; // liquid soil moisture available for evap/drainage (mm)
    double                     liq[MAX_LAYERS]; // current liquid soil moisture (mm)
    double                     ice[MAX_LAYERS]; // current frozen soil moisture (mm)
    double                     moist[MAX_LAYERS]; // current total soil moisture (liquid and frozen) (mm)
    double                     max_moist[MAX_LAYERS]; // maximum storable moisture (liquid and frozen) (mm)
    double                     Ksat[MAX_LAYERS];
    double                     Q12[MAX_LAYERS - 1];
    double                     Dsmax;
    double                     tmp_inflow;
    double                     tmp_moist;
    double                     tmp_moist_for_runoff[MAX_LAYERS];
    double                     tmp_liq;
    double                     dt_inflow;
    double                     dt_runoff;
    double                     runoff[MAX_FROST_AREAS];
    double                     tmp_dt_runoff[MAX_FROST_AREAS];
    double                     baseflow[MAX_FROST_AREAS];
    double                     dt_baseflow;
    double                     rel_moist;
    double                     evap[MAX_LAYERS][MAX_FROST_AREAS];
    double                     sum_liq;
    double                     evap_fraction;
    double                     evap_sum;
    layer_data_struct         *layer;
    layer_data_struct          tmp_layer;
    unsigned short             runoff_steps_per_dt;

    /** Set Residual Moisture **/
    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        resid_moist[lindex] = soil_con->resid_moist[lindex] *
                              soil_con->depth[lindex] * MM_PER_M;
    }

    /** Allocate and Set Values for Soil Sublayers **/
    layer = cell->layer;

    cell->runoff = 0;
    cell->baseflow = 0;
    cell->asat = 0;

    runoff_steps_per_dt = global_param.runoff_steps_per_day /
                          global_param.model_steps_per_day;

    for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
        baseflow[fidx] = 0;
    }

    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        evap[lindex][0] = layer[lindex].evap / (double) runoff_steps_per_dt;
        org_moist[lindex] = layer[lindex].moist;
        layer[lindex].moist = 0;
        if (evap[lindex][0] > 0) { // if there is positive evaporation
            sum_liq = 0;
            // compute available soil moisture for each frost sub area.
            for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
                avail_liq[lindex][fidx] =
                    (org_moist[lindex] - layer[lindex].ice[fidx] -
                     resid_moist[lindex]);
                if (avail_liq[lindex][fidx] < 0) {
                    avail_liq[lindex][fidx] = 0;
                }
                sum_liq += avail_liq[lindex][fidx] *
                           frost_fract[fidx];
            }
            // compute fraction of available soil moisture that is evaporated
            if (sum_liq > 0) {
                evap_fraction = evap[lindex][0] / sum_liq;
            }
            else {
                evap_fraction = 1.0;
            }
            // distribute evaporation between frost sub areas by percentage
            evap_sum = evap[lindex][0];
            for (fidx = (int)options.Nfrost - 1; fidx >= 0; fidx--) {
                evap[lindex][fidx] = avail_liq[lindex][fidx] * evap_fraction;
                avail_liq[lindex][fidx] -= evap[lindex][fidx];
                evap_sum -= evap[lindex][fidx] * frost_fract[fidx];
            }
        }
        else {
            for (fidx = (int)options.Nfrost - 1; fidx > 0; fidx--) {
                evap[lindex][fidx] = evap[lindex][0];
            }
        }
    }

    for (fidx = 0; fidx < (int)options.Nfrost; fidx++) {
        /** ppt = amount of liquid water coming to the surface **/
        inflow = ppt;

        /**************************************************
           Initialize Variables
        **************************************************/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            Ksat[lindex] = soil_con->Ksat[lindex] /
                           global_param.runoff_steps_per_day;

            /** Set Layer Liquid Moisture Content **/
            liq[lindex] = org_moist[lindex] - layer[lindex].ice[fidx];

            /** Set Layer Frozen Moisture Content **/
            ice[lindex] = layer[lindex].ice[fidx];

            /** Set Layer Maximum Moisture Content **/
            max_moist[lindex] = soil_con->max_moist[lindex];
        }

        /******************************************************
           Runoff Based on Soil Moisture Level of Upper Layers
        ******************************************************/

        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_moist_for_runoff[lindex] = (liq[lindex] + ice[lindex]);
        }
        compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, inflow, &A,
                                &(runoff[fidx]));

        // save dt_runoff based on initial runoff estimate,
        // since we will modify total runoff below for the case of completely saturated soil
        tmp_dt_runoff[fidx] = runoff[fidx] /
                              (double) runoff_steps_per_dt;

        /**************************************************
           Compute Flow Between Soil Layers ()
        **************************************************/

        dt_inflow = inflow / (double) runoff_steps_per_dt;

        Dsmax = soil_con->Dsmax / global_param.runoff_steps_per_day;

        for (time_step = 0; time_step < runoff_steps_per_dt; time_step++) {
            inflow = dt_inflow;

            /*************************************
               Compute Drainage between Sublayers
            *************************************/

            for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
                /** Brooks & Corey relation for hydraulic conductivity **/

                if ((tmp_liq = liq[lindex] - evap[lindex][fidx]) <
                    resid_moist[lindex]) {
                    tmp_liq = resid_moist[lindex];
                }

                if (tmp_liq > resid_moist[lindex]) {
                    Q12[lindex] = calc_Q12(Ksat[lindex], tmp_liq,
                                           resid_moist[lindex],
                                           soil_con->max_moist[lindex],
                                           soil_con->expt[lindex]);
                }
                else {
                    Q12[lindex] = 0.;
                }
            }

            /**************************************************
               Solve for Current Soil Layer Moisture, and
               Check Versus Maximum and Minimum Moisture Contents.
            **************************************************/

            last_index = 0;
            for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
                if (lindex == 0) {
                    dt_runoff = tmp_dt_runoff[fidx];
                }
                else {
                    dt_runoff = 0;
                }

                /* transport moisture for all sublayers **/

                tmp_inflow = 0.;

                /** Update soil layer moisture content **/
                liq[lindex] = liq[lindex] +
                              (inflow - dt_runoff) -
                              (Q12[lindex] + evap[lindex][fidx]);

                /** Verify that soil layer moisture is less than maximum **/
                if ((liq[lindex] + ice[lindex]) > max_moist[lindex]) {
                    tmp_inflow = (liq[lindex] + ice[lindex]) -
                                 max_moist[lindex];
                    liq[lindex] = max_moist[lindex] - ice[lindex];

                    if (lindex == 0) {
                        Q12[lindex] += tmp_inflow;
                        tmp_inflow = 0;
                    }
                    else {
                        tmplayer = lindex;
                        while (tmp_inflow > 0) {
                            tmplayer--;
                            if (tmplayer < 0) {
                                /** If top layer saturated, add to runoff **/
                                runoff[fidx] += tmp_inflow;
                                tmp_inflow = 0;
                            }
                            else {
                                /** else add excess soil moisture to next higher layer **/
                                liq[tmplayer] += tmp_inflow;
                                if ((liq[tmplayer] + ice[tmplayer]) >
                                    max_moist[tmplayer]) {
                                    tmp_inflow =
                                        ((liq[tmplayer] +
                                          ice[tmplayer]) - max_moist[tmplayer]);
                                    liq[tmplayer] = max_moist[tmplayer] -
                                                    ice[tmplayer];
                                }
                                else {
                                    tmp_inflow = 0;
                                }
                            }
                        }
                    } /** end trapped excess moisture **/
                } /** end check if excess moisture in top layer **/

                /** verify that current layer moisture is greater than minimum **/
                if (liq[lindex] < 0) {
                    /** liquid cannot fall below 0 **/
                    Q12[lindex] += liq[lindex];
                    liq[lindex] = 0;
                }
                if ((liq[lindex] + ice[lindex]) < resid_moist[lindex]) {
                    /** moisture cannot fall below minimum **/
                    Q12[lindex] +=
                        (liq[lindex] + ice[lindex]) - resid_moist[lindex];
                    liq[lindex] = resid_moist[lindex] - ice[lindex];
                }

                inflow = (Q12[lindex] + tmp_inflow);
                Q12[lindex] += tmp_inflow;

                last_index++;
            } /* end loop through soil layers */

            /**************************************************
               Compute Baseflow
            **************************************************/

            /** ARNO model for the bottom soil layer (based on bottom
                soil layer moisture from previous time step) **/

            lindex = options.Nlayer - 1;

            /** Compute relative moisture **/
            rel_moist =
                (liq[lindex] -
                 resid_moist[lindex]) /
                (soil_con->max_moist[lindex] - resid_moist[lindex]);

            /** Compute baseflow as function of relative moisture **/
            frac = Dsmax * soil_con->Ds / soil_con->Ws;
            dt_baseflow = frac * rel_moist;
            if (rel_moist > soil_con->Ws) {
                frac = (rel_moist - soil_con->Ws) / (1 - soil_con->Ws);
                dt_baseflow += Dsmax * (1 - soil_con->Ds / soil_con->Ws) * pow(
                    frac, soil_con->c);
            }

            /** Make sure baseflow isn't negative **/
            if (dt_baseflow < 0) {
                dt_baseflow = 0;
            }

            /** Extract baseflow from the bottom soil layer **/

            liq[lindex] +=
                Q12[lindex - 1] - (evap[lindex][fidx] + dt_baseflow);

            /** Check Lower Sub-Layer Moistures **/
            tmp_moist = 0;

            /* If soil moisture has gone below minimum, take water out
             * of baseflow and add back to soil to make up the difference
             * Note: this may lead to negative baseflow, in which case we will
             * reduce evap to make up for it */
            if ((liq[lindex] + ice[lindex]) < resid_moist[lindex]) {
                dt_baseflow +=
                    (liq[lindex] + ice[lindex]) - resid_moist[lindex];
                liq[lindex] = resid_moist[lindex] - ice[lindex];
            }

            if ((liq[lindex] + ice[lindex]) > max_moist[lindex]) {
                /* soil moisture above maximum */
                tmp_moist = ((liq[lindex] + ice[lindex]) - max_moist[lindex]);
                liq[lindex] = max_moist[lindex] - ice[lindex];
                tmplayer = lindex;
                while (tmp_moist > 0) {
                    tmplayer--;
                    if (tmplayer < 0) {
                        /** If top layer saturated, add to runoff **/
                        runoff[fidx] += tmp_moist;
                        tmp_moist = 0;
                    }
                    else {
                        /** else if sublayer exists, add excess soil moisture **/
                        liq[tmplayer] += tmp_moist;
                        if ((liq[tmplayer] + ice[tmplayer]) >
                            max_moist[tmplayer]) {
                            tmp_moist =
                                ((liq[tmplayer] +
                                  ice[tmplayer]) - max_moist[tmplayer]);
                            liq[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
                        }
                        else {
                            tmp_moist = 0;
                        }
                    }
                }
            }

            baseflow[fidx] += dt_baseflow;
        } /* end of sub-dt time step loop */

        /** If negative baseflow, reduce evap accordingly **/
        if (baseflow[fidx] < 0) {
            layer[lindex].evap += baseflow[fidx];
            baseflow[fidx] = 0;
        }

        /** Recompute Asat based on final moisture level of upper layers **/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_moist_for_runoff[lindex] = (liq[lindex] + ice[lindex]);
        }
        compute_runoff_and_asat(soil_con, tmp_moist_for_runoff, 0, &A,
                                &tmp_runoff);

        /** Store tile-wide values **/
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            layer[lindex].moist +=
                ((liq[lindex] + ice[lindex]) * frost_fract[fidx]);
        }
        cell->asat += A * frost_fract[fidx];
        cell->runoff += runoff[fidx] * frost_fract[fidx];
        cell->baseflow += baseflow[fidx] * frost_fract[fidx];
    }

    /** Compute water table depth **/
    wrap_compute_zwt(soil_con, cell);

    /** Recompute Thermal Parameters Based on New Moisture Distribution **/
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            tmp_layer = cell->layer[lindex];
            moist[lindex] = tmp_layer.moist;
        }

        ErrorFlag = distribute_node_moisture_properties(energy->moist,
                                                        energy->ice,
                                                        energy->kappa_node,
                                                        energy->Cs_node,
                                                        soil_con->Zsum_node,
                                                        energy->T,
                                                        soil_con->max_moist_node,
                                                        soil_con->expt_node,
                                                        soil_con->bubble_node,
                                                        moist, soil_con->depth,
                                                        soil_con->soil_dens_min,
                                                        soil_con->bulk_dens_min,
                                                        soil_con->quartz,
                                                        soil_con->soil_density,
                                                        soil_con->bulk_density,
                                                        soil_con->organic, Nnodes,
                                                        options.Nlayer,
                                                        soil_con->FS_ACTIVE);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }
    return (0);
}

/******************************************************************************
* @brief    Calculate the saturated area and runoff
******************************************************************************/
void
compute_runoff_and_asat(soil_con_struct *soil_con,
                        double          *moist,
                        double           inflow,
                        double          *A,
                        double          *runoff)
{
    extern option_struct options;
    double               top_moist; // total moisture (liquid and frozen) in topmost soil layers (mm)
    double               top_max_moist; // maximum storable moisture (liquid and frozen) in topmost soil layers (mm)
    size_t               lindex;
    double               ex;
    double               max_infil;
    double               i_0;
    double               basis;

    top_moist = 0.;
    top_max_moist = 0.;
    for (lindex = 0; lindex < options.Nlayer - 1; lindex++) {
        top_moist += moist[lindex];
        top_max_moist += soil_con->max_moist[lindex];
    }
    if (top_moist > top_max_moist) {
        top_moist = top_max_moist;
    }

    /** A as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
    ex = soil_con->b_infilt / (1.0 + soil_con->b_infilt);
    *A = 1.0 - pow((1.0 - top_moist / top_max_moist), ex);

    max_infil = (1.0 + soil_con->b_infilt) * top_max_moist;
    i_0 = max_infil * (1.0 - pow((1.0 - *A), (1.0 / soil_con->b_infilt)));

    /** equation (3a) Wood et al. **/

    if (inflow == 0.0) {
        *runoff = 0.0;
    }
    else if (max_infil == 0.0) {
        *runoff = inflow;
    }
    else if ((i_0 + inflow) > max_infil) {
        *runoff = inflow - top_max_moist + top_moist;
    }
    /** equation (3b) Wood et al. (wrong in paper) **/
    else {
        basis = 1.0 - (i_0 + inflow) / max_infil;
        *runoff = (inflow - top_max_moist + top_moist +
                   top_max_moist *
                   pow(basis, 1.0 * (1.0 + soil_con->b_infilt)));
    }
    if (*runoff < 0.) {
        *runoff = 0.;
    }
}

/******************************************************************************
* @brief    Calculate drainage between two layers
******************************************************************************/
double
calc_Q12(double Ksat,
         double init_moist,
         double resid_moist,
         double max_moist,
         double expt)
{
    double Q12;

    Q12 = init_moist - pow(pow(init_moist - resid_moist, 1.0 - expt) -
                           Ksat /
                           pow(max_moist - resid_moist, expt) * (1.0 - expt),
                           1.0 / (1.0 - expt)) - resid_moist;

    return Q12;
}
