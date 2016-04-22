/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initailizes the vegetation variable array.
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
 * @brief    This routine initailizes the vegetation variable array.
 *****************************************************************************/
void
initialize_veg(veg_var_struct **veg_var,
               size_t           Nveg)
{
    extern option_struct options;

    size_t               i, j, k;
    double               dummy_double;

    for (i = 0; i < Nveg; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Prognostic states
            veg_var[i][j].albedo = dummy_double;
            veg_var[i][j].displacement = dummy_double;
            veg_var[i][j].fcanopy = dummy_double;
            veg_var[i][j].LAI = dummy_double;
            veg_var[i][j].roughness = dummy_double;
            veg_var[i][j].Wdew = dummy_double;
            veg_var[i][j].Wdmax = dummy_double;
            // Fluxes
            veg_var[i][j].canopyevap = dummy_double;
            veg_var[i][j].throughfall = dummy_double;
        }
        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                // Carbon states
                veg_var[i][j].AnnualNPP = dummy_double;
                veg_var[i][j].AnnualNPPPrev = dummy_double;
                veg_var[i][j].Ci = dummy_double;
                veg_var[i][j].NPPfactor = dummy_double;
                veg_var[i][j].rc = dummy_double;
                for (k = 0; k < options.Ncanopy; k++) {
                    veg_var[i][j].CiLayer[k] = dummy_double;
                    veg_var[i][j].NscaleFactor[k] = dummy_double;
                    veg_var[i][j].rsLayer[k] = dummy_double;
                }
                // Carbon fluxes
                veg_var[i][j].aPAR = dummy_double;
                for (k = 0; k < options.Ncanopy; k++) {
                    veg_var[i][j].aPARLayer[k] = dummy_double;
                }
                veg_var[i][j].GPP = dummy_double;
                veg_var[i][j].Litterfall = dummy_double;
                veg_var[i][j].NPP = dummy_double;
                veg_var[i][j].Raut = dummy_double;
                veg_var[i][j].Rdark = dummy_double;
                veg_var[i][j].Rgrowth = dummy_double;
                veg_var[i][j].Rmaint = dummy_double;
                veg_var[i][j].Rphoto = dummy_double;
            }
        }
    }
}
