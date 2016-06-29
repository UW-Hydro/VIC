/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine uses Xu Liangs 3-layer energy balance formulation to estimate
 * the temperature between the first and second layers.  Formerly calculated
 * independently in each of the surface energy balance equation routines.
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
 * @brief    3-layer energy balance formulation to estimate the temperature
 *           between the first and second layers.
 *****************************************************************************/
double
estimate_T1(double Ts,
            double T1_old,
            double T2,
            double D1,
            double D2,
            double kappa1,
            double kappa2,
            double Cs2,
            double dp,
            double delta_t)
{
    double C1;
    double C2;
    double C3;
    double T1;

    C1 = Cs2 * dp / D2 * (1. - exp(-D2 / dp));
    C2 = -(1. - exp(D1 / dp)) * exp(-D2 / dp);
    C3 = kappa1 / D1 - kappa2 / D1 + kappa2 / D1 * exp(-D1 / dp);

    T1 = (kappa1 / 2. / D1 / D2 * (Ts) + C1 / delta_t * T1_old +
          (2. * C2 - 1. + exp(-D1 / dp)) * kappa2 / 2. / D1 / D2 * T2) /
         (C1 / delta_t + kappa2 / D1 / D2 * C2 + C3 / 2. / D2);

    return(T1);
}
