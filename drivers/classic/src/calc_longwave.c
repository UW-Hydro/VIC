/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine estimates long wave radiation.
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
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine estimates long wave radiation based on the fractional
 *           cloud cover (tskc), the current air temperature (C), and the
 *           atmospheric vapor pressure (Pa).
 *****************************************************************************/
void
calc_longwave(double *longwave,
              double  tskc,
              double  air_temp,
              double  vp)
{
    extern option_struct options;
    double               emissivity;
    double               emissivity_clear;
    double               cloudfactor;
    double               cloudfrac;
    double               x;

    air_temp += CONST_TKFRZ; // convert to Kelvin
    vp /= BAR_PER_KPA; // convert to mbar

    /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
       Addison-Wesley, 1990, p. 42-45 */

    if (options.LW_TYPE == LW_TVA) {
        emissivity_clear = 0.740 + 0.0049 * vp; // TVA (1972) - see Bras 2.35
    }
    else if (options.LW_TYPE == LW_ANDERSON) {
        emissivity_clear = 0.68 + 0.036 * pow(vp, 0.5); // Anderson (1964)
    }
    else if (options.LW_TYPE == LW_BRUTSAERT) {
        x = vp / air_temp;
        emissivity_clear = 1.24 * pow(x, 0.14285714); // Brutsaert (1975)
    }
    else if (options.LW_TYPE == LW_SATTERLUND) {
        emissivity_clear = 1.08 * (1 - exp(-1 * pow(vp, (air_temp / 2016)))); // Satterlund (1979)
    }
    else if (options.LW_TYPE == LW_IDSO) {
        emissivity_clear = 0.7 + 5.95e-5*vp*exp(1500 / air_temp); // Idso (1981)
    }
    else if (options.LW_TYPE == LW_PRATA) {
        x = 46.5 * vp / air_temp;
        emissivity_clear = 1 - (1 + x) * exp(-1 * pow((1.2 + 3 * x), 0.5)); // Prata (1996)
    }

    if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
        /* Assume emissivity of clouds is 1.0, and that total emissivity is weighted
           average of cloud emission plus clear-sky emission, weighted by fraction of
           sky occupied by each (method of Deardorff, 1978) */
        cloudfrac = tskc;
        emissivity = cloudfrac * 1.0 + (1 - cloudfrac) * emissivity_clear; // Deardorff (1978)
    }
    else {
        cloudfactor = 1.0 + 0.17 * tskc * tskc; // see Bras 2.43
        emissivity = cloudfactor * emissivity_clear;
    }

    *longwave = calc_outgoing_longwave(air_temp, emissivity);
}
