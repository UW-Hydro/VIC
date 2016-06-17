/******************************************************************************
* @section DESCRIPTION
*
* Calculate values related to the saturated vapor pressure.
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
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine computes the saturated vapor pressure
*
* @note         Handbook of Hydrology eqn 4.2.2.
******************************************************************************/
double
svp(double temp)
{
    extern parameters_struct param;

    double                   SVP;

    SVP = param.SVP_A * exp((param.SVP_B * temp) / (param.SVP_C + temp));

    if (temp < 0) {
        SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;
    }

    return (SVP * PA_PER_KPA);
}

/******************************************************************************
* @brief        This routine computes the gradient of d(svp)/dT
*
* @note         Handbook of Hydrology eqn 4.2.3
******************************************************************************/
double
svp_slope(double temp)
{
    extern parameters_struct param;

    return (param.SVP_B * param.SVP_C) / ((param.SVP_C + temp) *
                                          (param.SVP_C + temp)) * svp(temp);
}
