/******************************************************************************
 * @section DESCRIPTION
 *
 * Newton-Raphson method to solve non-linear system adapted from "Numerical
 * Recipes"
 *
 * A relaxation factor "RELAX#" is added to help the Newton-Raphson trials to
 * converge during the initial formation of ice, where the shape of the
 * residual function becomes very difficult.
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
 * @brief    Newton-Raphson method to solve non-linear system adapted from
 *           "Numerical Recipes"
 *****************************************************************************/
int
newt_raph(void (*vecfunc)(double x[], double fvec[], int n, int init, ...),
          double x[],
          int n)
{
    extern parameters_struct param;

    int                      k, i, Error;
    double                   errx, errf, fvec[MAX_NODES], p[MAX_NODES];
    double                   a[MAX_NODES], b[MAX_NODES], c[MAX_NODES];

    Error = 0;

    for (k = 0; k < param.NEWT_RAPH_MAXTRIAL; k++) {
        // calculate function value for all nodes, i.e. focus = -1
        (*vecfunc)(x, fvec, n, 0, -1);

        // stop if TOLF is satisfied
        errf = 0.0;
        for (i = 0; i < n; i++) {
            errf += fabs(fvec[i]);
        }
        if (errf <= param.NEWT_RAPH_TOLF) {
            return (Error);
        }

        // calculate the Jacobian
        fdjac3(x, fvec, a, b, c, vecfunc, n);

        for (i = 0; i < n; i++) {
            p[i] = -fvec[i];
        }

        // only the tri-diagnal part of the Jacobian matrix is significant
        // and most of the off-belt entries are encessially zeros
        tridiag(a, b, c, p, n);

        errx = 0.0;
        for (i = 0; i < n; i++) {
            errx += fabs(p[i]);

            if (k > 10 && k <= 20 && x[i] < param.NEWT_RAPH_R_MAX && x[i] >
                param.NEWT_RAPH_R_MIN) {
                x[i] += p[i] * param.NEWT_RAPH_RELAX1;
            }
            else if (k > 20 && k <= 60 && x[i] < param.NEWT_RAPH_R_MAX && x[i] >
                     param.NEWT_RAPH_R_MIN) {
                x[i] += p[i] * param.NEWT_RAPH_RELAX2;
            }
            else if (k > 60 && x[i] < param.NEWT_RAPH_R_MAX && x[i] >
                     param.NEWT_RAPH_R_MIN) {
                x[i] += p[i] * param.NEWT_RAPH_RELAX3;
            }
            else {
                x[i] += p[i];
            }
        }

        // stop if TOLX is satisfied
        if (errx <= param.NEWT_RAPH_TOLX) {
            return (Error);
        }
    }
    Error = 1;

    return (Error);
}

/******************************************************************************
 * @brief    Forward difference approx to Jacobian, adapted from "Numerical
 *           Recipes"
 *****************************************************************************/
void
fdjac3(double x[],
       double fvec[],
       double a[],
       double b[],
       double c[],
       void (*vecfunc)(double x[], double fvec[], int n, int init, ...),
       int n)
{
    extern parameters_struct param;

    int                      j;
    double                   h, temp, f[MAX_NODES];

    for (j = 0; j < n; j++) {
        temp = x[j];
        h = param.NEWT_RAPH_EPS2 * fabs(temp);
        if (h == 0) {
            h = param.NEWT_RAPH_EPS2;
        }
        x[j] = temp + h;
        h = x[j] - temp;

        // only update column j-1, j and j+1, caused by change in x[j]
        (*vecfunc)(x, f, n, 0, j);

        x[j] = temp;

        b[j] = (f[j] - fvec[j]) / h;
        if (j != 0) {
            c[j - 1] = (f[j - 1] - fvec[j - 1]) / h;
        }
        if (j != n - 1) {
            a[j + 1] = (f[j + 1] - fvec[j + 1]) / h;
        }
    }
}

/******************************************************************************
 * @brief    function to solve tridiagonal linear system adapted from
 *           "Numerical Recipes"
 *****************************************************************************/
void
tridiag(double   a[],
        double   b[],
        double   c[],
        double   r[],
        unsigned n)
{
    unsigned i;
    int      j;
    double   factor;

    /* forward substitution */
    factor = b[0];
    b[0] = 1.0;
    c[0] = c[0] / factor;
    r[0] = r[0] / factor;

    for (i = 1; i < n; i++) {
        factor = a[i];
        a[i] = a[i] - b[i - 1] * factor;
        b[i] = b[i] - c[i - 1] * factor;
        r[i] = r[i] - r[i - 1] * factor;

        factor = b[i];
        b[i] = 1.0;
        c[i] = c[i] / factor;
        r[i] = r[i] / factor;
    }

    /* backward substitution */
    for (j = n - 2; j >= 0; j--) {
        factor = c[j];
        c[j] = c[j] - b[j + 1] * factor;
        r[j] = r[j] - r[j + 1] * factor;

        factor = b[j];
        r[j] = r[j] / factor;
    }
}
