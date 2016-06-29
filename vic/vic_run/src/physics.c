/******************************************************************************
 * @section DESCRIPTION
 *
 * This set of functions calculate basic physical quantities that are
 * frequently calculated throughout the model
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
 * @brief    Compute the latent heat of sublimation.
 * @note     Eq. 3.19, Bras 1990
 *****************************************************************************/
double
calc_latent_heat_of_sublimation(double temp)
{
    double ls;

    ls = (677. - 0.07 * temp) * JOULES_PER_CAL * GRAMS_PER_KG;

    return(ls);
}

/******************************************************************************
 * @brief    This subroutine computes the latent heat of vaporization.
 * @note     Eq. 4.2.1 in Handbook of Hydrology.
 *****************************************************************************/
double
calc_latent_heat_of_vaporization(double temp)
{
    double lv;

    lv = CONST_LATVAP - 2361. * temp;

    return(lv);
}

/******************************************************************************
 * @brief    Compute outgoing longwave using the Stefan-Boltzman Law.
 *****************************************************************************/
double
calc_outgoing_longwave(double temp,
                       double emiss)
{
    double lwout;

    lwout = emiss * CONST_STEBOL * temp * temp * temp * temp;

    return(lwout);
}

/******************************************************************************
 * @brief    Calculate scale height based on average temperature in the column
 *****************************************************************************/
double
calc_scale_height(double tair,
                  double elevation)
{
    extern parameters_struct param;

    double                   h;

    h = CONST_RDAIR / CONST_G *
        ((tair + CONST_TKFRZ) + 0.5 * elevation * param.LAPSE_RATE);

    return(h);
}

/******************************************************************************
 * @brief    Compute the sensible heat flux.
 *****************************************************************************/
double
calc_sensible_heat(double atmos_density,
                   double t1,
                   double t0,
                   double Ra)
{
    double sensible;

    sensible = CONST_CPMAIR * atmos_density * (t1 - t0) / Ra;

    return(sensible);
}
