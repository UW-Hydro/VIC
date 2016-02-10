/******************************************************************************
* @section DESCRIPTION
*
* This subroutine computes dependent lake parameters from the specified
* lake parameters.
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
* @brief        This subroutine computes dependent lake parameters from the
*               specified lake parameters.
******************************************************************************/
void
compute_lake_params(lake_con_struct *lake_con,
                    soil_con_struct  soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;

    size_t                   i;
    double                   tempdz;
    double                   radius;
    double                   x;
    int                      ErrFlag;

    // miscellaneous lake parameters
    lake_con->bpercent = lake_con->rpercent;
    lake_con->maxdepth = lake_con->z[0];
    lake_con->basin[0] = lake_con->Cl[0] * soil_con.cell_area;

    if (!options.LAKE_PROFILE) {
        // generate lake depth-area relationship
        tempdz = (lake_con->maxdepth) / ((double) lake_con->numnod);
        radius = sqrt(lake_con->basin[0] / CONST_PI);

        for (i = 1; i <= lake_con->numnod; i++) {
            lake_con->z[i] = (lake_con->numnod - i) * tempdz;
            if (lake_con->z[i] < 0.0) {
                lake_con->z[i] = 0.0;
            }
            x =
                pow(lake_con->z[i] / lake_con->maxdepth,
                    param.LAKE_BETA) * radius;
            lake_con->basin[i] = CONST_PI * x * x;
            lake_con->Cl[i] = lake_con->basin[i] / soil_con.cell_area;
        }
    }
    else {
        // final point in depth-area relationship
        lake_con->z[lake_con->numnod] = 0;
        lake_con->Cl[lake_con->numnod] = 0;

        // depth-area relationship specified (for area fractions)
        // compute basin node surface areas
        for (i = 1; i <= lake_con->numnod; i++) {
            lake_con->basin[i] = lake_con->Cl[i] * soil_con.cell_area;
        }
    }

    // compute max volume
    lake_con->maxvolume = 0.0;
    for (i = 1; i <= lake_con->numnod; i++) {
        lake_con->maxvolume += (lake_con->basin[i] + lake_con->basin[i - 1]) *
                               (lake_con->z[i - 1] - lake_con->z[i]) / 2.;
    }

    // compute volume corresponding to mindepth
    ErrFlag = get_volume(*lake_con, lake_con->mindepth, &(lake_con->minvolume));
    if (ErrFlag == ERROR) {
        log_err("problem in get_volume(): depth %f volume %f",
                lake_con->mindepth, lake_con->minvolume);
    }
}
