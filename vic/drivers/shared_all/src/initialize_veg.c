/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initailizes the vegetation variable array.
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
 * @brief    This routine initailizes the vegetation variable array.
 *****************************************************************************/
void
initialize_veg(veg_var_struct **veg_var,
               size_t           Nveg)
{
    extern option_struct options;

    size_t               i, j, k;

    for (i = 0; i < Nveg; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Prognostic states
            veg_var[i][j].albedo = 0.0;
            veg_var[i][j].displacement = 0.0;
            veg_var[i][j].fcanopy = 0.0;
            veg_var[i][j].LAI = 0.0;
            veg_var[i][j].roughness = 0.0;
            veg_var[i][j].Wdew = 0.0;
            veg_var[i][j].Wdmax = 0.0;
            // Fluxes
            veg_var[i][j].canopyevap = 0.0;
            veg_var[i][j].throughfall = 0.0;
        }
        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                // Carbon states
                veg_var[i][j].AnnualNPP = 0.0;
                veg_var[i][j].AnnualNPPPrev = 0.0;
                veg_var[i][j].Ci = 0.0;
                veg_var[i][j].NPPfactor = 0.0;
                veg_var[i][j].rc = 0.0;
                for (k = 0; k < options.Ncanopy; k++) {
                    veg_var[i][j].CiLayer[k] = 0.0;
                    veg_var[i][j].NscaleFactor[k] = 0.0;
                    veg_var[i][j].rsLayer[k] = 0.0;
                }
                // Carbon fluxes
                veg_var[i][j].aPAR = 0.0;
                for (k = 0; k < options.Ncanopy; k++) {
                    veg_var[i][j].aPARLayer[k] = 0.0;
                }
                veg_var[i][j].GPP = 0.0;
                veg_var[i][j].Litterfall = 0.0;
                veg_var[i][j].NPP = 0.0;
                veg_var[i][j].Raut = 0.0;
                veg_var[i][j].Rdark = 0.0;
                veg_var[i][j].Rgrowth = 0.0;
                veg_var[i][j].Rmaint = 0.0;
                veg_var[i][j].Rphoto = 0.0;
            }
        }
    }
}
