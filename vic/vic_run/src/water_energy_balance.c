/******************************************************************************
* @section DESCRIPTION
*
* Calculate snow accumulation and melt for the lake model.
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
* @brief        Calculate snow accumulation and melt for the lake model.
******************************************************************************/
int
water_energy_balance(int     numnod,
                     double *surface,
                     double *evapw,
                     double  dt,
                     double  dz,
                     double  surfdz,
                     double  lat,
                     double  Tcutoff,
                     double  Tair,
                     double  wind,
                     double  pressure,
                     double  vp,
                     double  air_density,
                     double  longwave,
                     double  shortwave,
                     double  wind_h,
                     double *Qh,
                     double *Qle,
                     double *LWnet,
                     double *T,
                     double *water_density,
                     double *deltaH,
                     double *energy_ice_formation,
                     double  fracprv,
                     double *new_ice_area,
                     double *cp,
                     double *new_ice_height,
                     double *energy_out_bottom,
                     double *new_ice_water_eq,
                     double  lvolume)
{
    extern parameters_struct param;

    double                   Ts;
    double                   Tcutk;
    double                   Tskin;
    double                   Tnew[MAX_LAKE_NODES];
    int                      k;
    double                   sumjouli;

    double                   Tmean;
    int                      iterations;
    double                   Le;
    double                   jouleold;
    double                   joulenew;


    double                   de[MAX_LAKE_NODES];
    double                   epsilon = 0.0001;

    /* Calculate the surface energy balance for water surface temp = 0.0 */

    Tmean = MISSING;
    Ts = T[0];
    iterations = 0;
    for (k = 0; k < numnod; k++) {
        Tnew[k] = T[k];
    }

    energycalc(T, &jouleold, numnod, dz, surfdz, surface, cp, water_density);

    while ((fabs(Tmean - Ts) > epsilon) && iterations < param.LAKE_MAX_ITER) {
        if (iterations == 0) {
            Ts = T[0];
        }
        else {
            Ts = Tmean;
        }

        /* ....................................................................
         * Pass the skin temperature of the lake in Kelvin since the
         * Stefan-Boltzmann formula uses K.
         * ....................................................................*/

        Tskin = Ts + CONST_TKFRZ;
        Tcutk = Tcutoff + CONST_TKFRZ;

        /* ....................................................................
         * Send an ice height of 0 meters to latsens for the calculation
         * of latent and sensible heat over liquid water.
         * ....................................................................*/

        latsens(Tskin, Tcutk, 0.0, Tair, wind, pressure, vp, air_density,
                evapw, Qh, wind_h);

        /**********************************************************************
           Compute the Latent Heat Flux
        **********************************************************************/

        Le = calc_latent_heat_of_vaporization(Tair); /*J/kg */

        *Qle = -1. * (*evapw) * Le;                  /* W/m^2 */

        /* --------------------------------------------------------------------
         * Calculate the outgoing long wave fluxes, positive downwards.
         * -------------------------------------------------------------------- */

        *LWnet = longwave - calc_outgoing_longwave(Tskin, param.EMISS_H2O);

        /*************************************************************
           Use a Triadiagonal Matric to Explicitly Solve for
           Temperatures at Water Thermal Nodes
        *************************************************************/

        /* --------------------------------------------------------------------
         * Calculate the eddy diffusivity.
         * -------------------------------------------------------------------- */

        eddy(1, wind, water_density, de, lat, numnod, dz, surfdz);

        /* --------------------------------------------------------------------
         * Calculate the lake temperatures at different levels for the
         * new timestep.
         * -------------------------------------------------------------------- */

        temp_area(shortwave * param.LAKE_A1, shortwave * param.LAKE_A2,
                  *Qle + *Qh + *LWnet, T, Tnew,
                  water_density, de, dt, surface, numnod,
                  dz, surfdz, &joulenew, cp, energy_out_bottom);

        /* Surface temperature < 0.0, then ice will form. */
        if (Tnew[0] < Tcutoff) {
            iceform(energy_ice_formation, Tnew, Tcutoff, fracprv, new_ice_area,
                    numnod, dt, dz, surfdz, cp, surface, new_ice_height,
                    new_ice_water_eq, lvolume);

            energycalc(Tnew, &sumjouli, numnod, dz, surfdz, surface, cp,
                       water_density);
            *deltaH = (sumjouli - jouleold) /
                      (surface[0] * dt);
        }
        else {
            *deltaH = (joulenew - jouleold) /
                      (surface[0] * dt);
            *energy_ice_formation = 0.0;
        }

        Tmean = (Tnew[0] + T[0]) / 2.;

        iterations += 1;
    }

    if (fabs(Tmean - Ts) <= epsilon) {
        // Temperature reached convergence
        for (k = 0; k < numnod; k++) {
            T[k] = Tnew[k];
        }
        return(0);
    }
    else {
        Tskin = T[0] + CONST_TKFRZ;
        Tcutk = Tcutoff + CONST_TKFRZ;
        latsens(Tskin, Tcutk, 0.0, Tair, wind, pressure, vp, air_density,
                evapw, Qh, wind_h);

        Le = calc_latent_heat_of_vaporization(Tair); /*J/kg */
        *Qle = -1. * (*evapw) * Le;                  /* W/m^2 */

        *LWnet = longwave - calc_outgoing_longwave(Tskin, param.EMISS_H2O);

        if (T[0] < Tcutoff) {
            iceform(energy_ice_formation, T, Tcutoff, fracprv, new_ice_area,
                    numnod, dt, dz, surfdz, cp, surface, new_ice_height,
                    new_ice_water_eq, lvolume);
        }
        else {
            *energy_ice_formation = 0.0;
        }

        *deltaH = 0.0;
        return(0);
    }
}
