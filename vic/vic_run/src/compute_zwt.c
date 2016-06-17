/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute spatial average water table position (zwt).
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

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute spatial average water table position (zwt).  Water table
 *           position is measured in cm and is negative below the soil surface.
 *****************************************************************************/
double
compute_zwt(soil_con_struct *soil_con,
            int              lindex,
            double           moist)
{
    int    i;
    double zwt;

    zwt = MISSING;

    /** Compute zwt using soil moisture v zwt curve **/
    i = MAX_ZWTVMOIST - 1;
    while (i >= 1 && moist > soil_con->zwtvmoist_moist[lindex][i]) {
        i--;
    }
    if (i == MAX_ZWTVMOIST - 1) {
        if (moist < soil_con->zwtvmoist_moist[lindex][i]) {
            zwt = 999; // 999 indicates water table not present in this layer
        }
        else if (moist == soil_con->zwtvmoist_moist[lindex][i]) {
            zwt = soil_con->zwtvmoist_zwt[lindex][i]; // Just barely enough water for a water table
        }
    }
    else {
        zwt =
            soil_con->zwtvmoist_zwt[lindex][i +
                                            1] +
            (soil_con->zwtvmoist_zwt[lindex][i] -
             soil_con->zwtvmoist_zwt[lindex][i +
                                             1]) *
            (moist -
             soil_con->zwtvmoist_moist[lindex][i +
                                               1]) /
            (soil_con->zwtvmoist_moist[lindex][i] -
             soil_con->zwtvmoist_moist[lindex][i + 1]);                                                                                                                                                                                                    // interpolate to find water table level
    }

    return(zwt);
}

/******************************************************************************
 * @brief    Function to compute spatial average water table position (zwt) for
 *           individual layers as well as various total-column versions of zwt.
 *           Water table position is measured in cm and is negative below the
 *           soil surface.
 *****************************************************************************/
void
wrap_compute_zwt(soil_con_struct  *soil_con,
                 cell_data_struct *cell)
{
    extern option_struct options;

    size_t               lindex;
    short                idx;
    double               total_depth;
    double               tmp_depth;
    double               tmp_moist;

    /** Compute total soil column depth **/
    total_depth = 0;
    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        total_depth += soil_con->depth[lindex];
    }

    /** Compute each layer's zwt using soil moisture v zwt curve **/
    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        cell->layer[lindex].zwt =
            compute_zwt(soil_con, lindex, cell->layer[lindex].moist);
    }
    if (cell->layer[options.Nlayer - 1].zwt == 999) {
        cell->layer[options.Nlayer - 1].zwt = -total_depth * CM_PER_M;                                     // in cm
    }
    /** Compute total soil column's zwt; this will be the zwt of the lowest layer that isn't completely saturated **/
    idx = options.Nlayer - 1;
    tmp_depth = total_depth;
    while (idx >= 0 && soil_con->max_moist[idx] -
           cell->layer[idx].moist <= DBL_EPSILON) {
        tmp_depth -= soil_con->depth[idx];
        idx--;
    }
    if (idx < 0) {
        cell->zwt = 0;
    }
    else if (idx < (short) (options.Nlayer - 1)) {
        if (cell->layer[idx].zwt != 999) {
            cell->zwt = cell->layer[idx].zwt;
        }
        else {
            cell->zwt = -tmp_depth * CM_PER_M;
        }
    }
    else {
        cell->zwt = cell->layer[idx].zwt;
    }

    /** Compute total soil column's zwt_lumped; this will be the zwt of all N layers lumped together. **/
    tmp_moist = 0;
    for (lindex = 0; lindex < options.Nlayer; lindex++) {
        tmp_moist += cell->layer[lindex].moist;
    }
    cell->zwt_lumped = compute_zwt(soil_con, options.Nlayer + 1, tmp_moist);
    if (cell->zwt_lumped == 999) {
        cell->zwt_lumped = -total_depth * CM_PER_M;                      // in cm;
    }
}
