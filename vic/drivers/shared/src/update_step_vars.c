/******************************************************************************
* @section DESCRIPTION
*
* This subroutine updates data structures with values for the current
* time step.
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
******************************************************************************/

#include <vic_def.h>
#include <vic_run.h>

/******************************************************************************
* @brief        This subroutine controls the model core, it solves both the
*               energy and water balance models, as well as frozen soils.
******************************************************************************/
int
update_step_vars(all_vars_struct *all_vars,
                 veg_con_struct  *veg_con,
                 veg_hist_struct *veg_hist)
{
    extern option_struct options;

    unsigned short       iveg;
    size_t               Nveg;
    unsigned short       band;
    size_t               Nbands;
    veg_var_struct     **veg_var;

    /* set local pointers */
    veg_var = all_vars->veg_var;

    Nbands = options.SNOW_BAND;

    /* Set number of vegetation types */
    Nveg = veg_con[0].vegetat_type_num;

    /* Assign current veg characteristics */
    for (iveg = 0; iveg < Nveg; iveg++) {
        for (band = 0; band < Nbands; band++) {
            veg_var[iveg][band].vegcover = veg_hist[iveg].vegcover[0];
            veg_var[iveg][band].albedo = veg_hist[iveg].albedo[0];
            veg_var[iveg][band].LAI = veg_hist[iveg].LAI[0];
        }
    }

    return (0);
}
