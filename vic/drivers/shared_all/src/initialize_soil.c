/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the soil variable arrays for each new grid cell.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine initializes the soil variable arrays for each new
 *           grid cell.
 *****************************************************************************/
void
initialize_soil(cell_data_struct **cell,
                size_t             veg_num)
{
    extern option_struct options;

    size_t               veg, band, lindex, frost_area;

    for (veg = 0; veg <= veg_num; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            // Prognostic states
            cell[veg][band].aero_resist[0] = 0.0;
            cell[veg][band].aero_resist[1] = 0.0;
            cell[veg][band].CLitter = 0.0;
            cell[veg][band].CInter = 0.0;
            cell[veg][band].CSlow = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].Cs = 0.0;
                cell[veg][band].layer[lindex].T = 0.0;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg][band].layer[lindex].ice[frost_area] = 0.0;
                }
                cell[veg][band].layer[lindex].kappa = 0.0;
                cell[veg][band].layer[lindex].moist = 0.0;
                cell[veg][band].layer[lindex].phi = 0.0;
            }
            cell[veg][band].rootmoist = 0.0;
            cell[veg][band].wetness = 0.0;
            // Fluxes
            cell[veg][band].pot_evap = 0.0;
            cell[veg][band].baseflow = 0.0;
            cell[veg][band].runoff = 0.0;
            cell[veg][band].inflow = 0.0;
            cell[veg][band].RhLitter = 0.0;
            cell[veg][band].RhLitter2Atm = 0.0;
            cell[veg][band].RhInter = 0.0;
            cell[veg][band].RhSlow = 0.0;
            cell[veg][band].RhTot = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].esoil = 0.0;
                cell[veg][band].layer[lindex].transp = 0.0;
                cell[veg][band].layer[lindex].evap = 0.0;
            }
        }
    }
}
