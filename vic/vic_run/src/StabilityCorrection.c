/******************************************************************************
 * @section DESCRIPTION
 *
 * Stability correction
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
 * @brief    Calculate atmospheric stability correction for non-neutral
 *           conditions
 *****************************************************************************/
double
StabilityCorrection(double Z,
                    double d,
                    double TSurf,
                    double Tair,
                    double Wind,
                    double Z0)
{
    double Correction;          /* Correction to aerodynamic resistance */
    double Ri;                   /* Richardson's Number */
    double RiCr = 0.2;           /* Critical Richardson's Number */
    double RiLimit;              /* Upper limit for Richardson's Number */

    Correction = 1.0;

    /* Calculate the effect of the atmospheric stability using a Richardson
       Number approach */

    if (TSurf != Tair) {
        /* Non-neutral conditions */

        Ri = CONST_G * (Tair - TSurf) * (Z - d) /
             (((Tair +
                CONST_TKFRZ) + (TSurf + CONST_TKFRZ)) / 2.0 * Wind * Wind);

        RiLimit = (Tair + CONST_TKFRZ) /
                  (((Tair +
                     CONST_TKFRZ) +
                    (TSurf + CONST_TKFRZ)) / 2.0 * (log((Z - d) / Z0) + 5));

        if (Ri > RiLimit) {
            Ri = RiLimit;
        }

        if (Ri > 0.0) {
            Correction = (1 - Ri / RiCr) * (1 - Ri / RiCr);
        }
        else {
            if (Ri < -0.5) {
                Ri = -0.5;
            }

            Correction = sqrt(1 - 16 * Ri);
        }
    }

    return Correction;
}
