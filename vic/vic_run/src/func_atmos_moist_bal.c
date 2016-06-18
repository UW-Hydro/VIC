/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine solves the atmospheric exchange moisture balance.
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
 * @brief    This routine solves the atmospheric exchange moisture balance.
 *****************************************************************************/
double
func_atmos_moist_bal(double  VPcanopy,
                     va_list ap)
{
    double  InLatentHeat;
    double  Lv;
    double  Ra;
    double  atmos_density;
    double  gamma;
    double  vp; // atmospheric vapor pressure

    double *LatentHeat;

    // internal routine variables
    double  Error;

    // extract variables from va_arg
    InLatentHeat = (double)  va_arg(ap, double);
    Lv = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    gamma = (double)  va_arg(ap, double);
    vp = (double)  va_arg(ap, double);

    LatentHeat = (double *) va_arg(ap, double *);

    // compute sensible heat flux between canopy and atmosphere
    (*LatentHeat) = Lv * calc_sensible_heat(atmos_density, vp, VPcanopy,
                                            gamma * Ra);

    // compute energy balance error
    Error = InLatentHeat - (*LatentHeat);

    return (Error);
}
