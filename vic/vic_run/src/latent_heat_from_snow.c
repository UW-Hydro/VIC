/******************************************************************************
 * @section DESCRIPTION
 *
 * Computes the latent heat from the snowpack.
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
 * @brief    Computes the latent heat from the snowpack.
 *****************************************************************************/
void
latent_heat_from_snow(double  AirDens,
                      double  EactAir,
                      double  Lv,
                      double  Press,
                      double  Ra,
                      double  TMean,
                      double  Vpd,
                      double *LatentHeat,
                      double *LatentHeatSublimation,
                      double *VaporMassFlux,
                      double *BlowingMassFlux,
                      double *SurfaceMassFlux)
{
    double EsSnow;
    double Ls;

    EsSnow = svp(TMean);

    // SurfaceMassFlux and BlowingMassFlux in kg/m2s

    *SurfaceMassFlux = AirDens * (CONST_EPS / Press) * (EactAir - EsSnow) / Ra;

    if (Vpd == 0.0 && *SurfaceMassFlux < 0.0) {
        *SurfaceMassFlux = 0.0;
    }

    /* Calculate total latent heat flux */

    *VaporMassFlux = *SurfaceMassFlux + *BlowingMassFlux;

    if (TMean >= 0.0) {
        /* Melt conditions: use latent heat of vaporization */
        *LatentHeat = Lv * (*VaporMassFlux);
        *LatentHeatSublimation = 0;
    }
    else {
        /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
        Ls = calc_latent_heat_of_sublimation(TMean);
        *LatentHeatSublimation = Ls * (*VaporMassFlux);
        *LatentHeat = 0;
    }
}
