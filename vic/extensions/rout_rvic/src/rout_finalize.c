/******************************************************************************
 * @section DESCRIPTION
 *
 * clean up functions for routing extension
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

#include <rout.h>

/******************************************************************************
 * @brief    Finalize RVIC by freeing memory.
 *****************************************************************************/
void
rout_finalize(void)
{
    extern rout_struct rout;

    free(rout.rout_param.source2outlet_ind);
    free(rout.rout_param.source_time_offset);
    free(rout.rout_param.source_x_ind);
    free(rout.rout_param.source_y_ind);
    free(rout.rout_param.source_lat);
    free(rout.rout_param.source_lon);
    free(rout.rout_param.source_VIC_index);
    free(rout.rout_param.outlet_lat);
    free(rout.rout_param.outlet_lon);
    free(rout.rout_param.outlet_VIC_index);
    free(rout.rout_param.unit_hydrograph);
    free(rout.rout_param.aggrunin);
    free(rout.discharge);
    free(rout.ring);
}
