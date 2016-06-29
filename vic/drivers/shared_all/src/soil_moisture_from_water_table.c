/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute soil moistures for various values of water table depth.
 *
 * Here we use the relationship (e.g., Letts et al., 2000)
 * w(z) = { ((zwt-z)/bubble)**(-1/b), z <  zwt-bubble
 * { 1.0,                      z >= zwt-bubble
 * where
 * z      = depth below surface [cm]
 * w(z)   = relative moisture at depth z given by
 * (moist(z) - resid_moist) / (max_moist - resid_moist)
 * zwt    = depth of water table below surface [cm]
 * bubble = bubbling pressure [cm]
 * b      = 0.5*(expt-3)
 * Note that zwt-bubble = depth of the free water surface, i.e.
 * position below which soil is completely saturated.
 *
 * This assumes water in unsaturated zone above water table
 * is always in equilibrium between gravitational and matric
 * tension (e.g., Frolking et al, 2002).
 *
 * So, to find the soil moisture value in a layer corresponding
 * to a given water table depth zwt, we integrate w(z) over the
 * whole layer:
 *
 * w_avg = average w over whole layer = (integral of w*dz) / layer depth
 *
 * Then,
 * layer moisture = w_avg * (max_moist - resid_moist) + resid_moist
 *
 * Instead of the zwt defined above, will actually report free
 * water surface elevation zwt' = -(zwt-bubble).  I.e. zwt' < 0
 * below the soil surface, and marks the point of saturation
 * rather than pressure = 1 atm.
 *
 * Do this for each layer individually and also for a) the top N-1 layers
 * lumped together, and b) the entire soil column lumped together.
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
 * @brief    Compute soil moistures for various values of water table depth
 *****************************************************************************/
void
soil_moisture_from_water_table(soil_con_struct *soil_con,
                               size_t           nlayers)
{
    size_t i;
    size_t j;
    double b;
    double b_save;
    double bubble;
    double bub_save;
    double tmp_depth;
    double tmp_depth2;
    double tmp_depth2_save;
    double tmp_max_moist;
    double tmp_moist;
    double tmp_resid_moist;
    double w_avg;
    double zwt_prime;
    double zwt_prime_eff;

    // Individual layers
    tmp_depth = 0;
    for (j = 0; j < nlayers; j++) {
        b = 0.5 * (soil_con->expt[j] - 3);
        bubble = soil_con->bubble[j];
        tmp_resid_moist = soil_con->resid_moist[j] *
                          soil_con->depth[j] * MM_PER_M; // mm
        // depth of free water surface below top of layer (not yet elevation)
        zwt_prime = 0;
        for (i = 0; i < MAX_ZWTVMOIST; i++) {
            // elevation (cm) relative to soil surface
            soil_con->zwtvmoist_zwt[j][i] = -tmp_depth * CM_PER_M - zwt_prime;
            w_avg =
                (soil_con->depth[j] * CM_PER_M - zwt_prime -
                 (b / (b - 1)) * bubble *
                 (1 - pow((zwt_prime + bubble) / bubble, (b - 1) / b))) /
                (soil_con->depth[j] * CM_PER_M); // cm
            if (w_avg < 0) {
                w_avg = 0;
            }
            if (w_avg > 1) {
                w_avg = 1;
            }
            soil_con->zwtvmoist_moist[j][i] = w_avg *
                                              (soil_con->max_moist[j] -
                                               tmp_resid_moist) +
                                              tmp_resid_moist;
            zwt_prime += soil_con->depth[j] * CM_PER_M / (MAX_ZWTVMOIST - 1); // cm
        }
        tmp_depth += soil_con->depth[j];
    }

    /* Top N-1 layers lumped together (with average soil properties) */
    tmp_depth = 0;
    b = 0;
    bubble = 0;
    tmp_max_moist = 0;
    tmp_resid_moist = 0;
    for (j = 0; j < nlayers - 1; j++) {
        b += 0.5 * (soil_con->expt[j] - 3) * soil_con->depth[j];
        bubble += soil_con->bubble[j] * soil_con->depth[j];
        tmp_max_moist += soil_con->max_moist[j];   // total max_moist
        // total resid_moist in mm
        tmp_resid_moist += soil_con->resid_moist[j] * soil_con->depth[j] *
                           MM_PER_M;
        tmp_depth += soil_con->depth[j];
    }
    b /= tmp_depth;     // average b
    bubble /= tmp_depth;     // average bubble
    // depth of free water surface below top of layer (not yet elevation)
    zwt_prime = 0;
    for (i = 0; i < MAX_ZWTVMOIST; i++) {
        // elevation (cm) relative to soil surface
        soil_con->zwtvmoist_zwt[nlayers][i] = -zwt_prime;
        w_avg =
            (tmp_depth * CM_PER_M - zwt_prime - (b / (b - 1)) * bubble *
             (1 - pow((zwt_prime + bubble) / bubble, (b - 1) / b))) /
            (tmp_depth * CM_PER_M); // cm
        if (w_avg < 0) {
            w_avg = 0;
        }
        if (w_avg > 1) {
            w_avg = 1;
        }
        soil_con->zwtvmoist_moist[nlayers][i] =
            w_avg * (tmp_max_moist - tmp_resid_moist) + tmp_resid_moist;
        zwt_prime += tmp_depth * CM_PER_M / (MAX_ZWTVMOIST - 1); // in cm
    }

    // Compute zwt by taking total column soil moisture and filling column
    // from bottom up
    tmp_depth = 0;
    for (j = 0; j < nlayers; j++) {
        tmp_depth += soil_con->depth[j];
    }
    // depth of free water surface below soil surface (not yet elevation)
    zwt_prime = 0;
    for (i = 0; i < MAX_ZWTVMOIST; i++) {
        // elevation (cm) relative to soil surface
        soil_con->zwtvmoist_zwt[nlayers + 1][i] = -zwt_prime;
        // Integrate w_avg in pieces
        if (zwt_prime == 0) {
            tmp_moist = 0;
            for (j = 0; j < nlayers; j++) {
                tmp_moist += soil_con->max_moist[j];
            }
            soil_con->zwtvmoist_moist[nlayers + 1][i] = tmp_moist;
        }
        else {
            tmp_moist = 0;
            j = nlayers - 1;
            tmp_depth2 = tmp_depth - soil_con->depth[j];
            while (j > 0 && zwt_prime <= tmp_depth2 * CM_PER_M) {
                tmp_moist += soil_con->max_moist[j];
                j--;
                tmp_depth2 -= soil_con->depth[j];
            }
            w_avg =
                (tmp_depth2 * CM_PER_M + soil_con->depth[j] * CM_PER_M -
                 zwt_prime) / (soil_con->depth[j] * CM_PER_M);
            b = 0.5 * (soil_con->expt[j] - 3);
            bubble = soil_con->bubble[j];
            tmp_resid_moist = soil_con->resid_moist[j] * soil_con->depth[j] *
                              MM_PER_M;
            w_avg += -(b / (b - 1)) * bubble *
                     (1 -
                      pow((zwt_prime + bubble - tmp_depth2 * CM_PER_M) / bubble,
                          (b - 1) / b)) / (soil_con->depth[j] * CM_PER_M);
            tmp_moist += w_avg * (soil_con->max_moist[j] - tmp_resid_moist) +
                         tmp_resid_moist;
            b_save = b;
            bub_save = bubble;
            tmp_depth2_save = tmp_depth2;
            while (j > 0) {
                j--;
                tmp_depth2 -= soil_con->depth[j];
                b = 0.5 * (soil_con->expt[j] - 3);
                bubble = soil_con->bubble[j];
                tmp_resid_moist =
                    soil_con->resid_moist[j] * soil_con->depth[j] * MM_PER_M;
                zwt_prime_eff = tmp_depth2_save * CM_PER_M - bubble + bubble *
                                pow((zwt_prime + bub_save - tmp_depth2_save *
                                     CM_PER_M) / bub_save,
                                    b / b_save);
                w_avg = -(b / (b - 1)) * bubble *
                        (1 -
                         pow((zwt_prime_eff + bubble - tmp_depth2 *
                              CM_PER_M) / bubble,
                             (b - 1) / b)) / (soil_con->depth[j] * CM_PER_M);
                tmp_moist += w_avg *
                             (soil_con->max_moist[j] -
                              tmp_resid_moist) + tmp_resid_moist;
                b_save = b;
                bub_save = bubble;
                tmp_depth2_save = tmp_depth2;
            }
            soil_con->zwtvmoist_moist[nlayers + 1][i] = tmp_moist;
        }
        zwt_prime += tmp_depth * CM_PER_M / (MAX_ZWTVMOIST - 1); // in cm
    }
}
