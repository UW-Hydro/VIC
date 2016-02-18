/******************************************************************************
 * @section DESCRIPTION
 *
 * These routines compute lake dimensions.
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

#include <vic_def.h>
#include <vic_run.h>

/******************************************************************************
 * @brief    Function to compute surface area of liquid water in the lake,
 *           given the current depth of liquid water.
 *****************************************************************************/
int
get_sarea(lake_con_struct lake_con,
          double          depth,
          double         *sarea)
{
    size_t i;
    int    status;

    status = 0;
    *sarea = 0.0;

    if (depth > lake_con.z[0]) {
        *sarea = lake_con.basin[0];
    }
    else {
        for (i = 0; i < lake_con.numnod; i++) {
            if (depth <= lake_con.z[i] && depth > lake_con.z[i + 1]) {
                *sarea =
                    lake_con.basin[i +
                                   1] +
                    (depth -
                     lake_con.z[i +
                                1]) *
                    (lake_con.basin[i] -
                     lake_con.basin[i +
                                    1]) / (lake_con.z[i] - lake_con.z[i + 1]);
            }
        }
        if (*sarea == 0.0 && depth != 0.0) {
            status = ERROR;
        }
    }

    return status;
}

/******************************************************************************
 * @brief    Function to compute liquid water volume stored within the lake
 *           basin, given the current depth of liquid water.
 *****************************************************************************/
int
get_volume(lake_con_struct lake_con,
           double          depth,
           double         *volume)
{
    int    i;
    int    status;
    double m;

    status = 0;
    *volume = 0.0;

    if (depth > lake_con.z[0]) {
        status = 1;
        *volume = lake_con.maxvolume;
    }

    for (i = lake_con.numnod - 1; i >= 0; i--) {
        if (depth >= lake_con.z[i]) {
            *volume +=
                (lake_con.basin[i] +
                 lake_con.basin[i +
                                1]) * (lake_con.z[i] - lake_con.z[i + 1]) / 2.;
        }
        else if (depth < lake_con.z[i] && depth >= lake_con.z[i + 1]) {
            m =
                (lake_con.basin[i] -
                 lake_con.basin[i + 1]) / (lake_con.z[i] - lake_con.z[i + 1]);
            *volume +=
                (depth -
                 lake_con.z[i +
                            1]) *
                (m * (depth - lake_con.z[i + 1]) / 2. + lake_con.basin[i + 1]);
        }
    }

    if (*volume == 0.0 && depth != 0.0) {
        status = ERROR;
    }

    return status;
}

/******************************************************************************
 * @brief    Function to compute the depth of liquid water in the lake
 *           (distance between surface and deepest point), given volume of
 *           liquid water currently stored in lake.
 *****************************************************************************/
int
get_depth(lake_con_struct lake_con,
          double          volume,
          double         *depth)
{
    int    k;
    int    status;
    double m;
    double tempvolume;

    status = 0;

    if (volume < -1 * DBL_EPSILON) {
        volume = 0.0;
        status = 1;
    }

    if (volume >= lake_con.maxvolume) {
        *depth = lake_con.maxdepth;
        *depth += (volume - lake_con.maxvolume) / lake_con.basin[0];
    }
    else if (volume < DBL_EPSILON) {
        *depth = 0.0;
    }
    else {
        // Update lake depth
        *depth = 0.0;
        tempvolume = volume;
        for (k = lake_con.numnod - 1; k >= 0; k--) {
            if (tempvolume > ((lake_con.z[k] - lake_con.z[k + 1]) *
                              (lake_con.basin[k] +
                               lake_con.basin[k + 1]) / 2.)) {
                // current layer completely filled
                tempvolume -=
                    (lake_con.z[k] -
                     lake_con.z[k +
                                1]) *
                    (lake_con.basin[k] + lake_con.basin[k + 1]) / 2.;
                *depth += lake_con.z[k] - lake_con.z[k + 1];
            }
            else if (tempvolume > 0.0) {
                if (lake_con.basin[k] == lake_con.basin[k + 1]) {
                    *depth += tempvolume / lake_con.basin[k + 1];
                    tempvolume = 0.0;
                }
                else {
                    m =
                        (lake_con.basin[k] -
                         lake_con.basin[k +
                                        1]) /
                        (lake_con.z[k] - lake_con.z[k + 1]);
                    *depth +=
                        ((-1 *
                          lake_con.basin[k +
                                         1]) +
                         sqrt(lake_con.basin[k +
                                             1] *
                              lake_con.basin[k + 1] + 2. * m * tempvolume)) / m;
                    tempvolume = 0.0;
                }
            }
        }
        if (tempvolume / lake_con.basin[0] > DBL_EPSILON) {
            status = ERROR;
        }
    }

    if (*depth < 0.0 || (*depth == 0.0 && volume >= DBL_EPSILON)) {
        status = ERROR;
    }

    return status;
}
