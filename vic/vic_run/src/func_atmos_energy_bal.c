/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine solves the atmospheric exchange energy balance.
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
 * @brief    This routine solves the atmospheric exchange energy balance.
 *****************************************************************************/
double
func_atmos_energy_bal(double  Tcanopy,
                      va_list ap)
{
    double  Ra;
    double  Tair;
    double  atmos_density;
    double  InSensible;

    double *SensibleHeat;

    // internal routine variables
    double  Error;

    // extract variables from va_arg
    Ra = (double)  va_arg(ap, double);
    Tair = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    InSensible = (double)  va_arg(ap, double);
    SensibleHeat = (double *)va_arg(ap, double *);

    // compute sensible heat flux between canopy and atmosphere
    (*SensibleHeat) = calc_sensible_heat(atmos_density, Tair, Tcanopy, Ra);

    // compute energy balance error
    Error = InSensible - (*SensibleHeat);

    return (Error);
}
