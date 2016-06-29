/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate soil respiration (heterotrophic respiration, or Rh).
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
 * @brief    Calculate soil respiration (heterotrophic respiration, or Rh).
 *****************************************************************************/
void
compute_soil_resp(int     Nnodes,
                  double *dZ,
                  double  dZTot,
                  double  dt,
                  double *T,
                  double *w,
                  double  CLitter,
                  double  CInter,
                  double  CSlow,
                  double *RhLitter,
                  double *RhInterTot,
                  double *RhSlowTot)
{
    extern parameters_struct param;

    int                      i;
    double                   Tref;
    double                  *TK = NULL;
    double                   fTLitter;
    double                  *fTSoil = NULL;
    double                   fMLitter;
    double                  *fMSoil = NULL;
    double                  *CInterNode = NULL;
    double                  *CSlowNode = NULL;
    double                  *RhInter = NULL;
    double                  *RhSlow = NULL;

    /* Allocate temp arrays */
    TK = calloc(Nnodes, sizeof(*TK));
    check_alloc_status(TK, "Memory allocation error.");
    fTSoil = calloc(Nnodes, sizeof(*fTSoil));
    check_alloc_status(fTSoil, "Memory allocation error.");
    fMSoil = calloc(Nnodes, sizeof(*fMSoil));
    check_alloc_status(fMSoil, "Memory allocation error.");
    CInterNode = calloc(Nnodes, sizeof(*CInterNode));
    check_alloc_status(CInterNode, "Memory allocation error.");
    CSlowNode = calloc(Nnodes, sizeof(*CSlowNode));
    check_alloc_status(CSlowNode, "Memory allocation error.");
    RhInter = calloc(Nnodes, sizeof(*RhInter));
    check_alloc_status(RhInter, "Memory allocation error.");
    RhSlow = calloc(Nnodes, sizeof(*RhSlow));
    check_alloc_status(RhSlow, "Memory allocation error.");

    /* Compute Lloyd-Taylor temperature dependence */
    Tref = 10. + CONST_TKFRZ; /* reference temperature of 10 C */
    for (i = 0; i < Nnodes; i++) {
        TK[i] = T[i] + CONST_TKFRZ;
        if (TK[i] < param.SRESP_T0_LT) {
            TK[i] = param.SRESP_T0_LT;
        }
    }
    fTLitter =
        exp(param.SRESP_E0_LT *
            (1 / (Tref - param.SRESP_T0_LT) - 1 / (TK[0] - param.SRESP_T0_LT)));
    for (i = 0; i < Nnodes; i++) {
        fTSoil[i] =
            exp(param.SRESP_E0_LT *
                (1 /
                 (Tref - param.SRESP_T0_LT) - 1 / (TK[i] - param.SRESP_T0_LT)));
    }

    /* Compute moisture dependence */
    for (i = 0; i < Nnodes; i++) {
        if (w[i] < param.SRESP_WMINFM) {
            w[i] = param.SRESP_WMINFM;
        }
        if (w[i] > param.SRESP_WMAXFM) {
            w[i] = param.SRESP_WMAXFM;
        }
    }
    if (w[0] <= param.SRESP_WOPTFM) {
        fMLitter =
            (w[0] -
             param.SRESP_WMINFM) *
            (w[0] -
             param.SRESP_WMAXFM) /
            ((w[0] -
              param.SRESP_WMINFM) *
             (w[0] -
              param.SRESP_WMAXFM) -
             (w[0] - param.SRESP_WOPTFM) * (w[0] - param.SRESP_WOPTFM));
    }
    else {
        fMLitter = param.SRESP_RHSAT +
                   (1 -
                    param.SRESP_RHSAT) *
                   (w[0] -
                    param.SRESP_WMINFM) *
                   (w[0] -
                    param.SRESP_WMAXFM) /
                   ((w[0] -
                     param.SRESP_WMINFM) *
                    (w[0] -
                     param.SRESP_WMAXFM) -
                    (w[0] - param.SRESP_WOPTFM) * (w[0] - param.SRESP_WOPTFM));
    }
    if (fMLitter > 1.0) {
        fMLitter = 1.0;
    }
    if (fMLitter < 0.0) {
        fMLitter = 0.0;
    }
    for (i = 0; i < Nnodes; i++) {
        if (w[i] <= param.SRESP_WOPTFM) {
            fMSoil[i] =
                (w[i] -
                 param.SRESP_WMINFM) *
                (w[i] -
                 param.SRESP_WMAXFM) /
                ((w[i] -
                  param.SRESP_WMINFM) *
                 (w[i] -
                  param.SRESP_WMAXFM) -
                 (w[i] - param.SRESP_WOPTFM) * (w[i] - param.SRESP_WOPTFM));
        }
        else {
            fMSoil[i] = param.SRESP_RHSAT +
                        (1 -
                         param.SRESP_RHSAT) *
                        (w[i] -
                         param.SRESP_WMINFM) *
                        (w[i] -
                         param.SRESP_WMAXFM) /
                        ((w[i] -
                          param.SRESP_WMINFM) *
                         (w[i] -
                          param.SRESP_WMAXFM) -
                         (w[i] -
                          param.SRESP_WOPTFM) * (w[i] - param.SRESP_WOPTFM));
        }
        if (fMSoil[i] > 1.0) {
            fMSoil[i] = 1.0;
        }
        if (fMSoil[i] < 0.0) {
            fMSoil[i] = 0.0;
        }
    }

    /* Compute C per node */
    for (i = 0; i < Nnodes; i++) {
        CInterNode[i] = CInter * dZ[i] / dZTot;
        CSlowNode[i] = CSlow * dZ[i] / dZTot;
    }

    /* Compute Rh for various pools, nodes; C fluxes in [gC/m2d] */
    *RhLitter = param.SRESP_RFACTOR *
                (fTLitter * fMLitter /
                 (param.SRESP_TAULITTER * CONST_DDAYS_PER_YEAR * SEC_PER_DAY /
                  dt)) * CLitter;
    *RhInterTot = 0;
    *RhSlowTot = 0;
    for (i = 0; i < Nnodes; i++) {
        RhInter[i] = param.SRESP_RFACTOR *
                     (fTSoil[i] * fMSoil[i] /
                      (param.SRESP_TAUINTER * CONST_DDAYS_PER_YEAR *
                       HOURS_PER_DAY / dt)) * CInterNode[i];
        RhSlow[i] = param.SRESP_RFACTOR *
                    (fTSoil[i] * fMSoil[i] /
                     (param.SRESP_TAUSLOW * CONST_DDAYS_PER_YEAR *
                      HOURS_PER_DAY / dt)) * CSlowNode[i];
        *RhInterTot += RhInter[i];
        *RhSlowTot += RhSlow[i];
    }

    free((char *) TK);
    free((char *) fTSoil);
    free((char *) fMSoil);
    free((char *) CInterNode);
    free((char *) CSlowNode);
    free((char *) RhInter);
    free((char *) RhSlow);
}
