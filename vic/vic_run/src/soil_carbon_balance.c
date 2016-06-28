/******************************************************************************
* @section DESCRIPTION
*
* Calculate soil carbon balance
*
* @section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014  The Land Surface Hydrology Group, Department of Civil
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
* @brief    Calculate soil carbon balance.
*
* @Note     We assume that respiration in the litter pool results in a mix of
*           end products: CO2, which is lost to the atmosphere, and larger
*           molecules, which become part of the soil carbon pools (intermediate
*           and slow).  The fraction of respired carbon (by mass) that ends up
*           as CO2 is given by fAir.  All respiration from the intermediate and
*           slow pools is assumed to go to the atmosphere.
******************************************************************************/
void
soil_carbon_balance(soil_con_struct   *soil_con,
                    energy_bal_struct *energy,
                    cell_data_struct  *cell,
                    veg_var_struct    *veg_var)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern parameters_struct   param;

    size_t                     i;
    size_t                     Nnodes;
    double                    *dZ = NULL;
    double                    *dZCum = NULL;
    double                     dZTot;
    double                    *T = NULL;
    double                    *w = NULL;
    double                     tmp_double;
    double                     b;
    double                     wtd;
    double                     w0;
    double                     w1;

    // Find subset of thermal nodes that span soil hydrologic layers
    dZTot = 0;
    for (i = 0; i < options.Nlayer; i++) {
        dZTot += soil_con->depth[i];
    }
    i = 0;
    while (i < options.Nnode - 1 && soil_con->Zsum_node[i] < dZTot) {
        i++;
    }
    Nnodes = i;
    if (soil_con->Zsum_node[i] > dZTot) {
        Nnodes--;
    }
    dZ = calloc(Nnodes, sizeof(*dZ));
    check_alloc_status(dZ, "Memory allocation error");
    dZCum = calloc(Nnodes, sizeof(*dZCum));
    check_alloc_status(dZCum, "Memory allocation error");
    T = calloc(Nnodes, sizeof(*T));
    check_alloc_status(T, "Memory allocation error");
    w = calloc(Nnodes, sizeof(*w));
    check_alloc_status(w, "Memory allocation error");

    // Assign node thicknesses and temperatures for subset
    dZTot = 0;
    for (i = 0; i < Nnodes; i++) {
        dZ[i] = soil_con->dz_node[i] * MM_PER_M; // mm
        dZTot += dZ[i];
        dZCum[i] = dZTot;
        T[i] = energy->T[i];
    }

    // Compute node relative moistures based on lumped water table depth
    for (i = 0; i < Nnodes; i++) {
        b = 0.5 * (soil_con->expt_node[i] - 3);
        wtd = -(cell->zwt_lumped) * 10; // mm, positive downwards
        if (wtd > dZCum[i]) {
            if (i > 0) {
                w0 = pow(
                    (wtd + soil_con->bubble_node[i] -
                     dZCum[i - 1]) / soil_con->bubble_node[i], -1 / b);
            }
            else {
                w0 = pow(
                    (wtd + soil_con->bubble_node[i]) / soil_con->bubble_node[i],
                    -1 / b);
            }
            w1 = pow(
                (wtd + soil_con->bubble_node[i] -
                 dZCum[i]) / soil_con->bubble_node[i], -1 / b);
            w[i] = 0.5 * (w0 + w1);
        }
        else if ((i == 0 && wtd > 0) || (i > 0 && wtd > dZCum[i - 1])) {
            if (i > 0) {
                w0 = pow(
                    (wtd + soil_con->bubble_node[i] -
                     dZCum[i - 1]) / soil_con->bubble_node[i], -1 / b);
                tmp_double = 0.5 * (dZCum[i - 1] + wtd);
                w1 = pow(
                    (wtd + soil_con->bubble_node[i] -
                     tmp_double) / soil_con->bubble_node[i], -1 / b);
                w[i] =
                    (0.5 *
                     (w0 +
                      w1) *
                     (tmp_double -
                      dZCum[i -
                            1]) + 0.5 *
                     (w1 +
                      1.0) *
                     (wtd -
                      tmp_double) + 1.0 *
                     (dZCum[i] - wtd)) / (dZCum[i] - dZCum[i - 1]);
            }
            else {
                w0 = pow(
                    (wtd + soil_con->bubble_node[i]) / soil_con->bubble_node[i],
                    -1 / b);
                tmp_double = 0.5 * (0 + wtd);
                w1 = pow(
                    (wtd + soil_con->bubble_node[i] -
                     tmp_double) / soil_con->bubble_node[i], -1 / b);
                w[i] =
                    (0.5 *
                     (w0 +
                      w1) *
                     (tmp_double -
                      0) + 0.5 *
                     (w1 +
                      1.0) *
                     (wtd -
                      tmp_double) + 1.0 * (dZCum[i] - wtd)) / (dZCum[i] - 0);
            }
        }
        else {
            w[i] = 1.0;
        }
    }

    // Compute carbon fluxes out of soil (evaluate fluxes at subset of thermal nodes and sum them)
    compute_soil_resp(Nnodes, dZ, dZTot, global_param.dt, T, w, cell->CLitter,
                      cell->CInter, cell->CSlow, &(cell->RhLitter),
                      &(cell->RhInter), &(cell->RhSlow));
    cell->RhLitter2Atm = param.SRESP_FAIR * cell->RhLitter;
    cell->RhTot = cell->RhLitter2Atm + cell->RhInter + cell->RhSlow;

    // Compute balances of soil carbon pools
    // Assume previous year's NPP enters soil evenly throughout current year
    veg_var->Litterfall = veg_var->AnnualNPPPrev /
                          (CONST_DDAYS_PER_YEAR * SEC_PER_DAY /
                           global_param.dt);
    cell->CLitter += veg_var->Litterfall - cell->RhLitter;
    cell->CInter +=
        (1 -
         param.SRESP_FAIR) * cell->RhLitter * param.SRESP_FINTER -
        cell->RhInter;
    cell->CSlow +=
        (1 -
         param.SRESP_FAIR) * cell->RhLitter *
        (1 - param.SRESP_FINTER) - cell->RhSlow;

    // Free temporary dynamic arrays
    free((char *) dZ);
    free((char *) dZCum);
    free((char *) T);
    free((char *) w);
}
