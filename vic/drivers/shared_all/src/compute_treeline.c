/******************************************************************************
 * @section DESCRIPTION
 *
 * Compute treeline.
 *
 * This routine computes the annual average July temperature for the current
 * gridcell.  The temperature is than lapsed to determine the elevation at
 * which the annual average temperature is equal to 10C. Snow elevation bands
 * above this elevation are considered to be above the treeline.  When gridcell
 * data is output at the end of each time step, vegetation types with overstory
 * will be excluded from the variable averages of snow bands higher than the
 * treeline (e.g. decidous trees will be removed from high elevation snow bands,
 * while grass and shrubs will remain).  This is to serve as a preliminary fix
 * for high elevation "glaciers", a more permanent version would actually allow
 * for vegetation types to be excluded from various snow bands.
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
 * @brief    Compute treeline.
 *****************************************************************************/
void
compute_treeline(force_data_struct *force,
                 dmy_struct        *dmy,
                 double             avgJulyAirTemp,
                 double            *Tfactor,
                 bool              *AboveTreeLine)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern size_t              NF;

    double                     MonthSum;
    double                     AnnualSum;
    int                        MonthCnt;
    int                        AnnualCnt;
    unsigned int               rec;
    size_t                     band;
    size_t                     i;

    if (options.JULY_TAVG_SUPPLIED) {
        // use supplied average annual July air temperature
        AnnualSum = avgJulyAirTemp;
    }
    else {
        // compute average annual July air temperature from forcing
        AnnualSum = 0;
        AnnualCnt = 0;
        rec = 0;
        while (rec < global_param.nrecs) {
            if (dmy[rec].month == 7) {
                MonthSum = 0;
                MonthCnt = 0;
                while (dmy[rec].month == 7) {
                    for (i = 0; i < NF; i++) {
                        MonthSum += force[rec].air_temp[i];
                        MonthCnt++;
                    }
                    rec++;
                }
                if (MonthCnt > 0) {
                    // Sum monthly average July temperature
                    AnnualSum += MonthSum / (double) MonthCnt;
                    AnnualCnt++;
                }
            }
            rec++;
        }

        // Compute average annual July air temperature
        if (AnnualCnt > 0) {
            AnnualSum /= (double) AnnualCnt;
        }
    }

    // Lapse average annual July air temperature to 10C and determine elevation
    for (band = 0; band < options.SNOW_BAND; band++) {
        if (AnnualSum + Tfactor[band] <= 10.) {
            // Band is above treeline
            AboveTreeLine[band] = true;
        }
        else {
            AboveTreeLine[band] = false;
        }
    }
}
