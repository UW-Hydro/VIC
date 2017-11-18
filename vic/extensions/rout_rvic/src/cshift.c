/******************************************************************************
 * @section DESCRIPTION
 *
 * Function to shift columns or rows one position (used in convolution)
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
 * @brief   Function to shift columns or rows one position
 *****************************************************************************/
void
cshift(double *data,
       int     nx,
       int     ny,
       int     axis,
       int     direction)
{
    int    x, y;
    double b;

    if (axis == 0 && direction == 1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y);
            for (x = 0; x != nx - 1; x++) {
                *(data + y + ny * x) = *(data + y + ny * (x + 1));
            }
            *(data + y + ny * x) = b;
        }
    }

    if (axis == 0 && direction == -1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y + ny * (nx - 1));
            for (x = nx - 1; x >= 0; x--) {
                *(data + y + ny * (x + 1)) = *(data + y + ny * x);
            }
            *(data + y) = b;
        }
    }

    if (axis == 1 && direction == 1) {
        for (x = 0; x < nx; x++) {
            b = *(data + x * ny);
            for (y = 0; y != ny; y++) {
                *(data + y + ny * x) = *(data + y + 1 + ny * x);
            }
            *(data + y - 1 + ny * x) = b;
        }
    }

    if (axis == 1 && direction == -1) {
        for (x = 0; x < nx; x++) {
            b = *(data + ny - 1 + ny * x);
            for (y = ny - 2; y >= 0; y--) {
                *(data + y + 1 + ny * x) = *(data + y + ny * x);
            }
            *(data + x * ny) = b;
        }
    }
}
