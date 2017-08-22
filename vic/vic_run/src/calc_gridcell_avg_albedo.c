/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes gridcell-averaged albedo.
 *
 * Note: this routine is specifically designed for the CESM driver, for WRF,
 * but has been implemented in other drivers as well for the sake of consistency.
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
 * @brief    Compute gridcell-averaged albedo.
 *****************************************************************************/
void
calc_gridcell_avg_albedo(double             *albedo,
                         double              shortwave,
                         size_t              Nveg,
                         bool                overstory,
                         energy_bal_struct **energy,
                         snow_data_struct  **snow,
                         veg_con_struct     *veg_con,
                         soil_con_struct    *soil_con)
{
    extern option_struct options;
    size_t               veg;
    size_t               band;
    double               Cv;
    double               AreaFactor;
    double               TreeAdjustFactor = 1.;
    double               lakefactor = 1;
    double               swnet;

    swnet = 0;
    *albedo = 0;

    for (veg = 0; veg <= Nveg; veg++) {
        Cv = veg_con[veg].Cv;
        if (Cv > 0) {
            for (band = 0; band < options.SNOW_BAND; band++) {
                if (soil_con->AreaFract[band] > 0.) {
                    // TO-DO: account for treeline and lake factors
                    AreaFactor = (Cv * soil_con->AreaFract[band] *
                                  TreeAdjustFactor * lakefactor);
                    swnet += AreaFactor * energy[veg][band].NetShortAtmos;
                }
            }
        }
    }

    // compute gridcell-averaged albedo using average shortwave
    if (shortwave > 0) {
        // use average shortwave for albedo calculation
        *albedo = 1. - (swnet / shortwave);
    }
    else {
        // use vegetation, snow or bare soil albedo
        for (veg = 0; veg <= Nveg; veg++) {
            Cv = veg_con[veg].Cv;
            if (Cv > 0) {
                for (band = 0; band < options.SNOW_BAND; band++) {
                    if (soil_con->AreaFract[band] > 0.) {
                        // TO-DO: account for treeline and lake factors
                        AreaFactor = (Cv * soil_con->AreaFract[band] *
                                      TreeAdjustFactor * lakefactor);
                        if (snow[veg][band].snow && overstory) {
                            // use snow canopy albedo
                            *albedo += AreaFactor *
                                       energy[veg][band].AlbedoOver;
                        }
                        else {
                            // use surface albedo
                            *albedo += AreaFactor *
                                       energy[veg][band].AlbedoUnder;
                        }
                    }
                }
            }
        }
    }
}
