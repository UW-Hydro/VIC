/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute potential evaporation.
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
 * @brief    Compute potential evaporation.
 *****************************************************************************/
void
compute_pot_evap(size_t  model_steps_per_day,
                 double  rsmin,
                 double  albedo,
                 double  shortwave,
                 double  net_longwave,
                 double  RGL,
                 double  tair,
                 double  vpd,
                 double  lai,
                 double  elevation,
                 double *aero_resist,
                 char    overstory,
                 double  rarc,
                 double  fcanopy,
                 double  ra_soil,
                 double *pot_evap)
{
    extern parameters_struct param;

    double                   net_short;
    double                   gsm_inv;
    bool                     ref_crop;
    double                   rc;
    double                   net_rad;
    double                   ra_veg;
    double                   Epot_veg;
    double                   Epot_soil;

    /************************************************
       Estimate potential evap using penman equation
    ************************************************/

    net_short = (1.0 - albedo) * shortwave;
    gsm_inv = 1.0;
    ref_crop = false;
    rc = calc_rc(rsmin, net_short, RGL, tair, vpd, lai, gsm_inv, ref_crop);
    net_rad = net_short + net_longwave;
    if (!overstory) {
        ra_veg = aero_resist[0];
    }
    else {
        ra_veg = aero_resist[1];
    }
    Epot_veg =
        penman(tair, elevation, net_rad, vpd, ra_veg, rc, rarc) /
        model_steps_per_day;

    Epot_soil =
        penman(tair, elevation, net_rad, vpd, ra_soil, 0.0,
               param.SOIL_RARC) / model_steps_per_day;

    *pot_evap = fcanopy * Epot_veg + (1 - fcanopy) * Epot_soil;
}
