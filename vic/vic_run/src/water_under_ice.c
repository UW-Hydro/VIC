/******************************************************************************
* \section DESCRIPTION
*
* Calculates temperatures of water column under ice.
*
* \section LICENSE
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
* \brief        Calculates temperatures of water column under ice.
******************************************************************************/
int
water_under_ice(int     freezeflag,
                double  sw_ice,
                double  wind,
                double *Ti,
                double *water_density,
                double  lat,
                int     numnod,
                double  dz,
                double  surfdz,
                double  Tcutoff,
                double *qw,
                double *surface,
                double *deltaH,
                double *water_cp,
                int     mixdepth,
                double  hice,
                double  sdepth,
                double  dt,
                double *energy_out_bottom)
{
    extern parameters_struct param;

    double                   Tnew[MAX_LAKE_NODES];
    int                      k;
    int                      iterations;
    double                   jouleold;
    double                   joulenew;
    double                   de[MAX_LAKE_NODES];
    double                   epsilon = 0.0001;
    double                   qw_init, qw_mean, qw_final;
    double                   sw_underice_visible, sw_underice_nir;

    iterations = 0;

    for (k = 0; k < numnod; k++) {
        Tnew[k] = Ti[k];
    }

    // compute the eddy diffusivity
    eddy(freezeflag, wind, water_density, de, lat, numnod, dz, surfdz);

    // estimate the flux out of the water
    qw_init = 0.57 * (Ti[0] - Tcutoff) / (surfdz / 2.);
    *qw = qw_init;
    qw_mean = MISSING;

    energycalc(Ti, &jouleold, numnod, dz, surfdz, surface, water_cp,
               water_density);

    while ((fabs(qw_mean - *qw) >
            epsilon) && iterations < param.LAKE_MAX_ITER) {
        if (iterations == 0) {
            *qw = qw_init;
        }
        else {
            *qw = qw_mean;
        }

        // compute shortwave that transmitted through the lake ice
        sw_underice_visible = param.LAKE_A1 * sw_ice *
                              exp(-1. * (param.LAKE_LAMISW * hice +
                                         param.LAKE_LAMSSW * sdepth));
        sw_underice_nir = param.LAKE_A2 * sw_ice *
                          exp(-1. * (param.LAKE_LAMILW * hice +
                                     param.LAKE_LAMSLW * sdepth));

        /* --------------------------------------------------------------------
         * Calculate the lake temperatures at different levels for the
         * new timestep.
         * -------------------------------------------------------------------- */

        temp_area(sw_underice_visible, sw_underice_nir, -1. * (*qw), Ti, Tnew,
                  water_density, de, dt, surface, numnod,
                  dz, surfdz, &joulenew, water_cp, energy_out_bottom);

        // recompute storage of heat in the lake
        *deltaH = (joulenew - jouleold) / (surface[0] * dt);

        /* --------------------------------------------------------------------
         * Do the convective mixing of the lake water.
         * -------------------------------------------------------------------- */

        tracer_mixer(Tnew, &mixdepth, surface, numnod, dz, surfdz, water_cp);

        qw_final = 0.57 * (Tnew[0] - Tcutoff) / (surfdz / 2.);

        qw_mean = (qw_final + *qw) / 2.;

        iterations += 1;
    }

    if (fabs(qw_mean - *qw) <= epsilon) {
        // Temperature reached convergence
        for (k = 0; k < numnod; k++) {
            Ti[k] = Tnew[k];
        }
        *qw = qw_mean;
        return(0);
    }
    else {
        *qw = 0.0;
        for (k = 0; k < numnod; k++) {
            Ti[k] = Tcutoff;
        }
        energycalc(Ti, &joulenew, numnod, dz, surfdz, surface, water_cp,
                   water_density);
        *deltaH = (joulenew - jouleold) / (surface[0] * dt);
        return(0);
    }
}
