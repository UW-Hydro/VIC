/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate vegetation parameters.
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
 * @brief    This subroutine estimates the displacement height of vegetation
 *           with a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_displacement(double height)
{
    extern parameters_struct param;

    double                   value;

    value = param.VEG_RATIO_DH_HEIGHT * height;

    return (value);
}

/******************************************************************************
 * @brief    This subroutine backs the vegetation height out of the given
 *           displacement using the reverse procedure from cal_veg_displacement.
 *****************************************************************************/
double
calc_veg_height(double displacement)
{
    extern parameters_struct param;

    double                   value;

    value = displacement / param.VEG_RATIO_DH_HEIGHT;

    return (value);
}

/******************************************************************************
 * @brief    This subroutine estimates the roughness height of vegetation with
 *           a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_roughness(double height)
{
    extern parameters_struct param;

    double                   value;

    value = param.VEG_RATIO_RL_HEIGHT * height;

    return (value);
}
