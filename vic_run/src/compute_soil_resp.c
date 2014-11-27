/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate soil respiration (heterotrophic respiration, or Rh).
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
    int                  i;
    double               Tref;
    double              *TK;
    double               fTLitter;
    double              *fTSoil;
    double               fMLitter;
    double              *fMSoil;
    double              *CInterNode;
    double              *CSlowNode;
    double              *RhInter;
    double              *RhSlow;

    /* Allocate temp arrays */
    TK = (double*)calloc(Nnodes, sizeof(double));
    fTSoil = (double*)calloc(Nnodes, sizeof(double));
    fMSoil = (double*)calloc(Nnodes, sizeof(double));
    CInterNode = (double*)calloc(Nnodes, sizeof(double));
    CSlowNode = (double*)calloc(Nnodes, sizeof(double));
    RhInter = (double*)calloc(Nnodes, sizeof(double));
    RhSlow = (double*)calloc(Nnodes, sizeof(double));

    /* Compute Lloyd-Taylor temperature dependence */
    Tref = 10 + KELVIN; /* reference temperature of 10 C */
    for (i = 0; i < Nnodes; i++) {
        TK[i] = T[i] + KELVIN;
        if (TK[i] < T0_LT) {
            TK[i] = T0_LT;
        }
    }
    fTLitter = exp(E0_LT * (1 / (Tref - T0_LT) - 1 / (TK[0] - T0_LT)));
    for (i = 0; i < Nnodes; i++) {
        fTSoil[i] = exp(E0_LT * (1 / (Tref - T0_LT) - 1 / (TK[i] - T0_LT)));
    }

    /* Compute moisture dependence */
    for (i = 0; i < Nnodes; i++) {
        if (w[i] < wminFM) {
            w[i] = wminFM;
        }
        if (w[i] > wmaxFM) {
            w[i] = wmaxFM;
        }
    }
    if (w[0] <= woptFM) {
        fMLitter =
            (w[0] -
             wminFM) *
            (w[0] -
             wmaxFM) /
            ((w[0] -
              wminFM) * (w[0] - wmaxFM) - (w[0] - woptFM) * (w[0] - woptFM));
    }
    else {
        fMLitter = Rhsat +
                   (1 -
                    Rhsat) *
                   (w[0] -
                    wminFM) *
                   (w[0] -
                    wmaxFM) /
                   ((w[0] -
                     wminFM) *
                    (w[0] - wmaxFM) - (w[0] - woptFM) * (w[0] - woptFM));
    }
    if (fMLitter > 1.0) {
        fMLitter = 1.0;
    }
    if (fMLitter < 0.0) {
        fMLitter = 0.0;
    }
    for (i = 0; i < Nnodes; i++) {
        if (w[i] <= woptFM) {
            fMSoil[i] =
                (w[i] -
                 wminFM) *
                (w[i] -
                 wmaxFM) /
                ((w[i] -
                  wminFM) *
                 (w[i] - wmaxFM) - (w[i] - woptFM) * (w[i] - woptFM));
        }
        else {
            fMSoil[i] = Rhsat +
                        (1 -
                         Rhsat) *
                        (w[i] -
                         wminFM) *
                        (w[i] -
                         wmaxFM) /
                        ((w[i] -
                          wminFM) *
                         (w[i] - wmaxFM) - (w[i] - woptFM) * (w[i] - woptFM));
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
    *RhLitter = Rfactor *
                (fTLitter * fMLitter /
                 (tauLitter * 365.25 * 24 / dt)) * CLitter;
    *RhInterTot = 0;
    *RhSlowTot = 0;
    for (i = 0; i < Nnodes; i++) {
        RhInter[i] = Rfactor *
                     (fTSoil[i] * fMSoil[i] /
                      (tauInter * 365.25 * 24 / dt)) * CInterNode[i];
        RhSlow[i] = Rfactor *
                    (fTSoil[i] * fMSoil[i] /
                     (tauSlow * 365.25 * 24 / dt)) * CSlowNode[i];
        *RhInterTot += RhInter[i];
        *RhSlowTot += RhSlow[i];
    }

    free((char *)TK);
    free((char *)fTSoil);
    free((char *)fMSoil);
    free((char *)CInterNode);
    free((char *)CSlowNode);
    free((char *)RhInter);
    free((char *)RhSlow);
}
