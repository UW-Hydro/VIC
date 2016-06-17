/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculates mass release of snow from canopy.
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
 * @brief    Calculates mass release of snow from canopy.
 *****************************************************************************/
void
MassRelease(double *InterceptedSnow,
            double *TempInterceptionStorage,
            double *ReleasedMass,
            double *Drip)
{
    extern parameters_struct param;

    double                   TempDrip;
    double                   TempReleasedMass;
    double                   Threshold;
    double                   MaxRelease;

    /* If the amount of snow in the canopy is greater than some minimum
       value, MIN_INTERCEPTION_STORAGE, then calculte mass release and Drip */

    if (*InterceptedSnow > param.VEG_MIN_INTERCEPTION_STORAGE) {
        Threshold = 0.10 * *InterceptedSnow;
        MaxRelease = 0.17 * *InterceptedSnow;

        /* If the amount of snow_melt after interception, snow_melt, is >= the
           theshhold then there is mass release.  If snow_melt is < the treshhold
           then there is no mass release but that water remains in
         * TempInterceptionStorage which will be augmented during the next
           compute period */

        if ((*TempInterceptionStorage) >= Threshold) {
            *Drip += Threshold;
            *InterceptedSnow -= Threshold;
            *TempInterceptionStorage -= Threshold;
            if (*InterceptedSnow < param.VEG_MIN_INTERCEPTION_STORAGE) {
                TempReleasedMass = 0.0;
            }
            else {
                TempReleasedMass =
                    min((*InterceptedSnow - param.VEG_MIN_INTERCEPTION_STORAGE),
                        MaxRelease);
            }
            *ReleasedMass += TempReleasedMass;
            *InterceptedSnow -= TempReleasedMass;
            MassRelease(InterceptedSnow, TempInterceptionStorage, ReleasedMass,
                        Drip);
        }
        else {
            TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
            *Drip += TempDrip;
            *InterceptedSnow -= TempDrip;
        }
    }

    /* (*InterceptedSnow < MIN_INTERCEPTION_STORAGE) If the amount of snow in
       the canopy is less than some minimum value, MIN_INTERCEPTION_STORAGE,
       then only melt can occur and there is no mass release. */

    else {
        TempDrip = min(*TempInterceptionStorage, *InterceptedSnow);
        *Drip += TempDrip;
        *InterceptedSnow -= TempDrip;
        *TempInterceptionStorage = 0.0;
    }
}
