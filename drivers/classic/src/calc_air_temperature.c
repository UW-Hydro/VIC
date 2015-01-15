/******************************************************************************
 * @section DESCRIPTION
 *
 * Estimate the daily temperature cycle from maximum and minimum daily
 * temperature measurements.
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
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    calculate the coefficients for the Hermite polynomials
 *****************************************************************************/
void
hermite(size_t  n,
        double *x,
        double *yc1,
        double *yc2,
        double *yc3,
        double *yc4)
{
    size_t i;
    double dx;
    double divdf1;
    double divdf3;

    for (i = 0; i < n - 1; i++) {
        dx = x[i + 1] - x[i];
        divdf1 = (yc1[i + 1] - yc1[i]) / dx;
        divdf3 = yc2[i] + yc2[i + 1] - 2 * divdf1;
        yc3[i] = (divdf1 - yc2[i] - divdf3) / dx;
        yc4[i] = divdf3 / (dx * dx);
    }
}

/******************************************************************************
 * @brief    use the Hermite polynomials, to find the interpolation function
 *           value at xbar
 *****************************************************************************/
double
hermint(double  xbar,
        int     n,
        double *x,
        double *yc1,
        double *yc2,
        double *yc3,
        double *yc4)
{
    int    klo, khi, k;
    double dx;
    double result;

    klo = 0;
    khi = n - 1;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (x[k] > xbar) {
            khi = k;
        }
        else {
            klo = k;
        }
    }

    dx = xbar - x[klo];
    result = yc1[klo] + dx * (yc2[klo] + dx * (yc3[klo] + dx * yc4[klo]));
    return result;
}

/******************************************************************************
 * @brief    calculate subdaily temperature.
 *****************************************************************************/
void
SubDailyT(size_t  stepsperday,
          size_t  ndays,
          double *TmaxSec,
          double *Tmax,
          double *TminSec,
          double *Tmin,
          double *Tair)
{
    double *x;
    double *Tyc1;
    double *yc2;
    double *yc3;
    double *yc4;
    double  sec;
    double  Dt;
    size_t  i;
    size_t  j;
    size_t  n;
    size_t  nsteps;

    nsteps = stepsperday * ndays;
    Dt = (double) (SEC_PER_DAY) / (double) stepsperday;

    n = ndays * 2 + 2;
    x = (double *) calloc(n, sizeof(double));
    Tyc1 = (double *) calloc(n, sizeof(double));
    yc2 = (double *) calloc(n, sizeof(double));
    yc3 = (double *) calloc(n, sizeof(double));
    yc4 = (double *) calloc(n, sizeof(double));
    if (x == NULL || Tyc1 == NULL || yc2 == NULL || yc3 == NULL || yc4 ==
        NULL) {
        log_err("Memory allocation failure in SubDailyT()");
    }

    /* First fill the x vector with the times for Tmin and Tmax, and fill the
       Tyc1 with the corresponding temperature and humidity values */
    for (i = 0, j = 1, sec = 0; i < ndays; i++, sec += SEC_PER_DAY) {
        if (TminSec[i] < TmaxSec[i]) {
            x[j] = TminSec[i] + sec;
            Tyc1[j++] = Tmin[i];
            x[j] = TmaxSec[i] + sec;
            Tyc1[j++] = Tmax[i];
        }
        else {
            x[j] = TmaxSec[i] + sec;
            Tyc1[j++] = Tmax[i];
            x[j] = TminSec[i] + sec;
            Tyc1[j++] = Tmin[i];
        }
    }

    /* To "tie" down the first and last values, repeat those */
    x[0] = x[2] - SEC_PER_DAY;
    Tyc1[0] = Tyc1[2];
    x[n - 1] = x[n - 3] + SEC_PER_DAY;
    Tyc1[n - 1] = Tyc1[n - 3];

    /* we want to preserve maxima and minima, so we require that the first
       derivative at these points is zero */
    for (i = 0; i < n; i++) {
        yc2[i] = 0.;
    }

    /* calculate the coefficients for the splines for the temperature */
    hermite(n, x, Tyc1, yc2, yc3, yc4);

    /* interpolate the temperatures */
    for (i = 0, sec = 0; i < nsteps; i++, sec += Dt) {
        Tair[i] = hermint(sec, n, x, Tyc1, yc2, yc3, yc4);
    }

    free(x);
    free(Tyc1);
    free(yc2);
    free(yc3);
    free(yc4);

    return;
}

/******************************************************************************
 * @brief    This function estimates the times of minimum and maximum
 *           temperature for each day of the simulation, based on the diurnal
 *           cycle of incoming solar radiation.
 *****************************************************************************/
void
set_max_min_sec(double *subdailyrad,
                size_t  ndays,
                size_t  stepsperday,
                double *tmaxsec,
                double *tminsec)
{
    double stepsize;
    size_t step;
    double risesec;
    double setsec;
    double sec;
    size_t i;

    stepsize = (double) (SEC_PER_DAY) / (double) stepsperday;

    for (i = 0; i < ndays; i++) {
        risesec = MISSING;
        setsec = MISSING;
        for (step = 0, sec = 0;
             sec < 12 * SEC_PER_HOUR;
             step++, sec += stepsize) {
            if (subdailyrad[i * stepsperday + step] > 0 &&
                (i * stepsperday + step == 0 ||
                 subdailyrad[i * stepsperday + step - 1] <= 0)) {
                risesec = sec;
            }
        }
        for (step = step + 1, sec = 12 * SEC_PER_HOUR;
             sec < SEC_PER_DAY;
             step++, sec += stepsize) {
            if (subdailyrad[i * stepsperday + step] <= 0 &&
                subdailyrad[i * stepsperday + step - 1] > 0) {
                setsec = sec;
            }
        }
        if (i == ndays - 1 && setsec == MISSING) {
            setsec = 23 * SEC_PER_HOUR;
        }
        if (risesec >= 0 && setsec >= 0) {
            tmaxsec[i] = 0.67 * (setsec - risesec) + risesec;
            if (risesec > stepsize) {
                tminsec[i] = risesec - stepsize;
            }
            else {
                tminsec[i] = 0.;
            }
        }
        else {
            /* arbitrarily set the min and max times to 2am and 2pm */
            tminsec[i] = 2 * SEC_PER_HOUR;
            tmaxsec[i] = 14 * SEC_PER_HOUR;
        }
    }
}
