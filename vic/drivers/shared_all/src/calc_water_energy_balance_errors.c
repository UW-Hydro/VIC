/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the overall model energy and water balance errors.
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
 * @brief    This subroutine computes the overall model water balance, and
 *           warns the model user if large errors are found.
 *****************************************************************************/
double
calc_water_balance_error(int    rec,
                         double inflow,
                         double outflow,
                         double storage)
{
    static double last_storage;
    static double cum_error;
    static double max_error;
    static int    Nrecs;

    double        error;

    if (rec < 0) {
        last_storage = storage;
        cum_error = 0.;
        max_error = 0.;
        Nrecs = -rec;

        return(0.0);
    }
    else {
        error = inflow - outflow - (storage - last_storage);
        cum_error += error;
        if (fabs(error) > fabs(max_error) && fabs(error) > 1e-5) {
            max_error = error;
            fprintf(LOG_DEST, "Maximum Moist Error:\t%i\t%.5f\t%.5f\n",
                    rec, error, cum_error);
        }
        if (rec == Nrecs - 1) {
            fprintf(LOG_DEST,
                    "Total Cumulative Water Error for Grid Cell = %.4f\n",
                    cum_error);
        }
        last_storage = storage;

        return(error);
    }
}

/******************************************************************************
 * @brief    This subroutine computes the overall model energy balance, and
 *           reports the maximum time step error above a thresehold to the
 *           user.  The total cumulative error for the grid cell is also
 *           computed and reported at the end of the model run.
 *****************************************************************************/
double
calc_energy_balance_error(int    rec,
                          double net_rad,
                          double latent,
                          double sensible,
                          double grnd_flux,
                          double snow_fluxes)
{
    static double cum_error;
    static double max_error;
    static int    Nrecs;

    double        error;

    if (rec < 0) {
        cum_error = 0;
        Nrecs = -rec;
        max_error = 0;
        error = 0.0;
    }
    else {
        error = net_rad - latent - sensible - grnd_flux + snow_fluxes;
        cum_error += error;
        if (fabs(error) > fabs(max_error) && fabs(error) > 0.001) {
            max_error = error;
            if (rec > 0) {
                fprintf(LOG_DEST, "Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
                        rec, error, cum_error / (double) rec);
            }
            else {
                fprintf(LOG_DEST, "Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
                        rec, error, cum_error);
            }
        }
        if (rec == Nrecs - 1) {
            fprintf(LOG_DEST,
                    "Total Cumulative Energy Error for Grid Cell = %.4f\n",
                    cum_error / (double) rec);
        }
    }

    return(error);
}
