/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the model state (energy balance, water balance, and
 * snow components).
 *
 * If a state file is provided to the model than its
 * contents are checked to see if it agrees with the current simulation set-up,
 * if so it is used to initialize the model state.  If no state file is
 * provided the model initializes all variables with defaults and the user
 * should expect to throw out the beginning of the simulation period as model
 * start-up.
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
 * @brief    Initialize the model state (energy balance, water balance, and
 * snow components).
 *****************************************************************************/
int
initialize_model_state(all_vars_struct *all_vars,
                       size_t           Nveg,
                       size_t           Nnodes,
                       double           surf_temp,
                       soil_con_struct *soil_con,
                       veg_con_struct  *veg_con)
{
    extern global_param_struct global_param;
    extern option_struct       options;

    char                       FIRST_VEG;
    size_t                     veg;
    size_t                     index;
    size_t                     lidx;
    size_t                     band;
    size_t                     frost_area;
    int                        ErrorFlag;
    double                     Cv;
    double                     Zsum, dp;
    double                     tmpdp, tmpadj, Bexp;
    double                     moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
    double                     ice[MAX_VEG][MAX_BANDS][MAX_LAYERS][
        MAX_FROST_AREAS];
    double                     dt_thresh;

    cell_data_struct         **cell;
    energy_bal_struct        **energy;

    cell = all_vars->cell;
    energy = all_vars->energy;

    // Initialize soil depths
    dp = soil_con->dp;

    FIRST_VEG = true;

    // increase initial soil surface temperature if air is very cold
    if (surf_temp < -1.) {
        surf_temp = -1.;
    }

    /************************************************************************
       CASE 1: initial conditions files provided -- handled in vic_restore()
    ************************************************************************/

    /************************************************************************
       CASE 2: Initialize soil if using quick heat flux, and no initial
       soil properties file given
    ************************************************************************/

    if (!options.INIT_STATE && options.QUICK_FLUX) {
        Nnodes = options.Nnode;

        /* Initialize soil node thicknesses */
        soil_con->dz_node[0] = soil_con->depth[0];
        soil_con->dz_node[1] = soil_con->depth[0];
        soil_con->dz_node[2] = 2. * (dp - 1.5 * soil_con->depth[0]);
        soil_con->Zsum_node[0] = 0;
        soil_con->Zsum_node[1] = soil_con->depth[0];
        soil_con->Zsum_node[2] = dp;

        for (veg = 0; veg <= Nveg; veg++) {
            // Initialize soil for existing vegetation types
            Cv = veg_con[veg].Cv;

            if (Cv > 0) {
                for (band = 0; band < options.SNOW_BAND; band++) {
                    /* Initialize soil node temperatures */
                    energy[veg][band].T[0] = surf_temp;
                    energy[veg][band].T[1] = surf_temp;
                    energy[veg][band].T[2] = soil_con->avg_temp;

                    /* Initialize soil layer moisture and ice contents */
                    for (lidx = 0; lidx < options.Nlayer; lidx++) {
                        moist[veg][band][lidx] =
                            cell[veg][band].layer[lidx].moist;
                        for (frost_area = 0;
                             frost_area < options.Nfrost;
                             frost_area++) {
                            ice[veg][band][lidx][frost_area] = 0.;
                        }
                    }
                }
            }
        }
    }

    /*****************************************************************
       CASE 3: Initialize Energy Balance Variables if not using quick
       ground heat flux, and no Initial Condition File Given
    *****************************************************************/
    else if (!options.INIT_STATE && !options.QUICK_FLUX) {
        for (veg = 0; veg <= Nveg; veg++) {
            // Initialize soil for existing vegetation types
            Cv = veg_con[veg].Cv;

            if (Cv > 0) {
                for (band = 0; band < options.SNOW_BAND; band++) {
                    if (!options.EXP_TRANS) {
                        /* Initialize soil node temperatures and thicknesses
                           Nodes set at surface, the depth of the first layer,
                           twice the depth of the first layer, and at the
                           damping depth.  Extra nodes are placed equal distance
                           between the damping depth and twice the depth of the
                           first layer. */

                        soil_con->dz_node[0] = soil_con->depth[0];
                        soil_con->dz_node[1] = soil_con->depth[0];
                        soil_con->dz_node[2] = soil_con->depth[0];

                        soil_con->Zsum_node[0] = 0;
                        soil_con->Zsum_node[1] = soil_con[0].depth[0];
                        Zsum = 2. * soil_con[0].depth[0];
                        soil_con->Zsum_node[2] = Zsum;
                        tmpdp = dp - soil_con[0].depth[0] * 2.5;
                        tmpadj = 3.5;
                        for (index = 3; index < Nnodes - 1; index++) {
                            if (FIRST_VEG) {
                                soil_con->dz_node[index] = tmpdp /
                                                           (((double) Nnodes -
                                                             tmpadj));
                            }
                            Zsum += (soil_con->dz_node[index] +
                                     soil_con->dz_node[index - 1]) / 2.;
                            soil_con->Zsum_node[index] = Zsum;
                        }
                        energy[veg][band].T[0] = surf_temp;
                        for (index = 1; index < Nnodes; index++) {
                            energy[veg][band].T[index] = soil_con->avg_temp;
                        }
                        if (FIRST_VEG) {
                            FIRST_VEG = false;
                            soil_con->dz_node[Nnodes - 1] = (dp - Zsum -
                                                             soil_con->dz_node[
                                                                 Nnodes - 2] /
                                                             2.) * 2.;
                            Zsum += (soil_con->dz_node[Nnodes - 2] +
                                     soil_con->dz_node[Nnodes - 1]) / 2.;
                            soil_con->Zsum_node[Nnodes - 1] = Zsum;
                            if ((int) (Zsum * MM_PER_M + 0.5) !=
                                (int) (dp * MM_PER_M + 0.5)) {
                                log_err("Sum of thermal node thicknesses (%f) "
                                        "in initialize_model_state do not "
                                        "equal dp (%f), check initialization "
                                        "procedure", Zsum, dp);
                            }
                        }
                    }
                    else { // exponential grid transformation, EXP_TRANS = TRUE
                           // calculate exponential function parameter */
                        if (FIRST_VEG) {
                            // to force Zsum=dp at bottom node
                            Bexp = logf(dp + 1.) / (double) (Nnodes - 1);
                            // validate Nnodes by requiring that there be at
                            // least 3 nodes in the top 50cm
                            if (Nnodes < 5 * logf(dp + 1.) + 1) {
                                log_err("The number of soil thermal nodes (%zu) "
                                        "is too small for the supplied damping "
                                        "depth (%f) with EXP_TRANS set to "
                                        "TRUE, leading to fewer than 3 nodes "
                                        "in the top 50 cm of the soil column.  "
                                        "For EXP_TRANS=TRUE, Nnodes and dp "
                                        "must follow the relationship:\n"
                                        "5*ln(dp+1)<Nnodes-1\n"
                                        "Either set Nnodes to at least %d in "
                                        "the global param file or reduce "
                                        "damping depth to %f in the soil "
                                        "parameter file.  Or set EXP_TRANS to "
                                        "FALSE in the global parameter file.",
                                        Nnodes, dp,
                                        (int) (5 * logf(
                                                   dp + 1.)) + 2,
                                        exp(0.2 * (Nnodes - 1)) + 1);
                            }
                            for (index = 0; index <= Nnodes - 1; index++) {
                                soil_con->Zsum_node[index] =
                                    expf(Bexp * index) - 1.;
                            }
                            if (soil_con->Zsum_node[0] > soil_con->depth[0]) {
                                log_err("Depth of first thermal node (%f) in "
                                        "initialize_model_state is greater "
                                        "than depth of first soil layer (%f); "
                                        "increase the number of nodes or "
                                        "decrease the thermal damping depth "
                                        "dp (%f)",
                                        soil_con->Zsum_node[0],
                                        soil_con->depth[0], dp);
                            }
                        }

                        // top node
                        index = 0;
                        if (FIRST_VEG) {
                            soil_con->dz_node[index] =
                                soil_con->Zsum_node[index +
                                                    1] -
                                soil_con->Zsum_node[index];
                        }
                        energy[veg][band].T[index] = surf_temp;
                        // middle nodes
                        for (index = 1; index < Nnodes - 1; index++) {
                            if (FIRST_VEG) {
                                soil_con->dz_node[index] =
                                    (soil_con->Zsum_node[index +
                                                         1] -
                                     soil_con->Zsum_node[index]) / 2. +
                                    (soil_con->Zsum_node[index] -
                                     soil_con->Zsum_node[index - 1]) / 2.;
                            }
                            energy[veg][band].T[index] = soil_con->avg_temp;
                        }
                        // bottom node
                        index = Nnodes - 1;
                        if (FIRST_VEG) {
                            soil_con->dz_node[index] =
                                soil_con->Zsum_node[index] -
                                soil_con->Zsum_node[index - 1];
                        }
                        energy[veg][band].T[index] = soil_con->avg_temp;
                    } // end if !EXP_TRANS

                    // initialize moisture and ice for each soil layer
                    for (lidx = 0; lidx < options.Nlayer; lidx++) {
                        moist[veg][band][lidx] =
                            cell[veg][band].layer[lidx].moist;
                        for (frost_area = 0;
                             frost_area < options.Nfrost;
                             frost_area++) {
                            ice[veg][band][lidx][frost_area] = 0.;
                        }
                    }
                }
            }
        }
    }

    /*********************************
       CASE 4: Unknown option
    *********************************/
    else if (!options.INIT_STATE) {
        for (veg = 0; veg <= Nveg; veg++) {
            // Initialize soil for existing vegetation types
            Cv = veg_con[veg].Cv;

            if (Cv > 0) {
                for (band = 0; band < options.SNOW_BAND; band++) {
                    // Initialize soil for existing snow elevation bands
                    if (soil_con->AreaFract[band] > 0.) {
                        for (index = 0; index < options.Nlayer; index++) {
                            soil_con->dz_node[index] = 1.;
                        }
                    }
                }
            }
        }
    }

    /******************************************
       Initialize soil thermal node properties
    ******************************************/

    FIRST_VEG = true;
    for (veg = 0; veg <= Nveg; veg++) {
        // Initialize soil for existing vegetation types
        Cv = veg_con[veg].Cv;

        if (Cv > 0) {
            for (band = 0; band < options.SNOW_BAND; band++) {
                // Initialize soil for existing snow elevation bands
                if (soil_con->AreaFract[band] > 0.) {
                    /** Set soil properties for all soil nodes **/
                    if (FIRST_VEG) {
                        FIRST_VEG = false;
                        set_node_parameters(soil_con->Zsum_node,
                                            soil_con->max_moist_node,
                                            soil_con->expt_node,
                                            soil_con->bubble_node,
                                            soil_con->alpha, soil_con->beta,
                                            soil_con->gamma, soil_con->depth,
                                            soil_con->max_moist, soil_con->expt,
                                            soil_con->bubble,
                                            Nnodes, options.Nlayer);
                    }

                    // set soil moisture properties for all soil thermal nodes
                    ErrorFlag =
                        distribute_node_moisture_properties(
                            energy[veg][band].moist,
                            energy[veg][band].ice,
                            energy[veg][band].kappa_node,
                            energy[veg][band].Cs_node,
                            soil_con->Zsum_node,
                            energy[veg][band].T,
                            soil_con->max_moist_node,
                            soil_con->expt_node,
                            soil_con->bubble_node,
                            moist[veg][band],
                            soil_con->depth,
                            soil_con->soil_dens_min,
                            soil_con->bulk_dens_min,
                            soil_con->quartz,
                            soil_con->soil_density,
                            soil_con->bulk_density,
                            soil_con->organic,
                            Nnodes, options.Nlayer,
                            soil_con->FS_ACTIVE);
                    if (ErrorFlag == ERROR) {
                        return (ErrorFlag);
                    }

                    // Check node spacing v time step
                    // (note this is only approximate since heat capacity and
                    // conductivity can change considerably during the
                    // simulation depending on soil moisture and ice content)
                    if ((options.FROZEN_SOIL &&
                         !options.QUICK_FLUX) && !options.IMPLICIT) {
                        // in seconds
                        dt_thresh = 0.5 * energy[veg][band].Cs_node[1] /
                                    energy[veg][band].kappa_node[1] *
                                    pow((soil_con->dz_node[1]),
                                        2);
                        if (global_param.dt > dt_thresh) {
                            log_err("You are currently running FROZEN SOIL "
                                    "with an explicit method (IMPLICIT is "
                                    "set to FALSE).  For the explicit method "
                                    "to be stable, time step %f seconds is too "
                                    "large for the given thermal node spacing "
                                    "%f m, soil heat capacity %f J/m3/K, and "
                                    "soil thermal conductivity %f J/m/s/K.  "
                                    "Either set IMPLICIT to TRUE in your "
                                    "global parameter file (this is the "
                                    "recommended action), or decrease time "
                                    "step length to <= %f seconds, or decrease "
                                    "the number of soil thermal nodes.",
                                    global_param.dt,
                                    soil_con->dz_node[1],
                                    energy[veg][band].Cs_node[1],
                                    energy[veg][band].kappa_node[1], dt_thresh);
                        }
                    }

                    /* initialize layer moistures and ice contents */
                    for (lidx = 0; lidx < options.Nlayer; lidx++) {
                        cell[veg][band].layer[lidx].moist =
                            moist[veg][band][lidx];
                        for (frost_area = 0;
                             frost_area < options.Nfrost;
                             frost_area++) {
                            cell[veg][band].layer[lidx].ice[frost_area] =
                                ice[veg][band][lidx][frost_area];
                        }
                    }
                    if (options.QUICK_FLUX) {
                        ErrorFlag =
                            estimate_layer_ice_content_quick_flux(
                                cell[veg][band].layer,
                                soil_con->depth, soil_con->dp,
                                energy[
                                    veg][band].T[0], energy[veg][band].T[1],
                                soil_con->avg_temp, soil_con->max_moist,
                                soil_con->expt, soil_con->bubble,
                                soil_con->frost_fract, soil_con->frost_slope,
                                soil_con->FS_ACTIVE);
                    }
                    else {
                        ErrorFlag = estimate_layer_ice_content(
                            cell[veg][band].layer,
                            soil_con->Zsum_node,
                            energy[veg][band].T,
                            soil_con->depth,
                            soil_con->max_moist,
                            soil_con->expt,
                            soil_con->bubble,
                            soil_con->frost_fract,
                            soil_con->frost_slope,
                            Nnodes, options.Nlayer,
                            soil_con->FS_ACTIVE);
                    }

                    /* Find freezing and thawing front depths */
                    if (!options.QUICK_FLUX && soil_con->FS_ACTIVE) {
                        find_0_degree_fronts(&energy[veg][band],
                                             soil_con->Zsum_node,
                                             energy[veg][band].T, Nnodes);
                    }
                }
            }
        }
    }

    return(0);
}
