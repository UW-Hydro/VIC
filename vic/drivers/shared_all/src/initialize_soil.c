/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the soil variable arrays for each new grid cell.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine initializes the soil variable arrays for each new
 *           grid cell.
 *****************************************************************************/
void
initialize_soil(cell_data_struct **cell,
                soil_con_struct   *soil_con,
                size_t             veg_num)
{
    extern option_struct options;

    size_t               veg, band, lindex, frost_area;
    double               dummy_double;

    for (veg = 0; veg <= veg_num; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            // Prognostic states
            cell[veg][band].aero_resist[0] = dummy_double;
            cell[veg][band].aero_resist[1] = dummy_double;
            cell[veg][band].CLitter = dummy_double;
            cell[veg][band].CInter = dummy_double;
            cell[veg][band].CSlow = dummy_double;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].Cs = dummy_double;
                cell[veg][band].layer[lindex].T = dummy_double;
                cell[veg][band].layer[lindex].bare_evap_frac = dummy_double;
                cell[veg][band].layer[lindex].evap = dummy_double;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg][band].layer[lindex].ice[frost_area] = dummy_double;
                }
                cell[veg][band].layer[lindex].kappa = dummy_double;
                cell[veg][band].layer[lindex].moist = dummy_double;
                cell[veg][band].layer[lindex].phi = dummy_double;
            }
            cell[veg][band].rootmoist = dummy_double;
            cell[veg][band].wetness = dummy_double;
            // Fluxes
            cell[veg][band].pot_evap = dummy_double;
            cell[veg][band].baseflow = dummy_double;
            cell[veg][band].runoff = dummy_double;
            cell[veg][band].inflow = dummy_double;
            cell[veg][band].RhLitter = dummy_double;
            cell[veg][band].RhLitter2Atm = dummy_double;
            cell[veg][band].RhInter = dummy_double;
            cell[veg][band].RhSlow = dummy_double;
            cell[veg][band].RhTot = dummy_double;
        }
    }
}
