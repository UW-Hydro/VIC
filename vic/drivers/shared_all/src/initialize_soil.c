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
    double               tmp_moist[MAX_LAYERS];
    double               tmp_runoff;

    for (veg = 0; veg <= veg_num; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            cell[veg][band].baseflow = 0;
            cell[veg][band].runoff = 0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].evap = 0;
                cell[veg][band].layer[lindex].moist =
                    soil_con->init_moist[lindex];
                if (cell[veg][band].layer[lindex].moist >
                    soil_con->max_moist[lindex]) {
                    cell[veg][band].layer[lindex].moist =
                        soil_con->max_moist[lindex];
                }
                tmp_moist[lindex] = cell[veg][band].layer[lindex].moist;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg][band].layer[lindex].ice[frost_area] = 0;
                }
            }
            compute_runoff_and_asat(soil_con, tmp_moist, 0,
                                    &(cell[veg][band].asat), &tmp_runoff);
            wrap_compute_zwt(soil_con, &(cell[veg][band]));
            cell[veg][band].CLitter = 0;
            cell[veg][band].CInter = 0;
            cell[veg][band].CSlow = 0;
        }
    }
    for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
        if (options.Nfrost == 1) {
            soil_con->frost_fract[frost_area] = 1.;
        }
        else if (options.Nfrost == 2) {
            soil_con->frost_fract[frost_area] = 0.5;
        }
        else {
            soil_con->frost_fract[frost_area] = 1. / (options.Nfrost - 1);
            if (frost_area == 0 || frost_area == options.Nfrost - 1) {
                soil_con->frost_fract[frost_area] /= 2.;
            }
        }
    }
}
