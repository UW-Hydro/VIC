/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate energy of sublimation from blowing snow.
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
 * @brief    Calculate sublimation from blowing snow
 *****************************************************************************/
double
CalcBlowingSnow(double   Dt,
                double   Tair,
                unsigned LastSnow,
                double   SurfaceLiquidWater,
                double   Wind,
                double   Ls,
                double   AirDens,
                double   EactAir,
                double   ZO,
                double   Zrh,
                double   snowdepth,
                double   lag_one,
                double   sigma_slope,
                double   Tsnow,
                int      iveg,
                int      Nveg,
                double   fe,
                double   displacement,
                double   roughness,
                double  *TotalTransport)
{
    extern parameters_struct param;
    extern option_struct     options;

    /* Local variables: */
    double                   Age;
    double                   U10, Uo, prob_occurence;
    double                   es, Ros, F;
    double                   SubFlux;
    double                   Diffusivity;
    double                   ushear;
    double                   Tk;
    double                   utshear;
    int                      p;
    double                   upper, lower, Total;
    double                   area;
    double                   sigma_w;
    double                   Zo_salt;
    double                   ratio, wind10;
    double                   Uveg, hv, Nd;
    double                   Transport;

    /*******************************************************************/
    /* Calculate some general variables, that don't depend on wind speed. */

    /* Age in hours */
    Age = LastSnow * Dt / SEC_PER_HOUR;

    /* Saturation density of water vapor, Liston A-8 */
    es = svp(Tair);

    Tk = Tair + CONST_TKFRZ;

    Ros = CONST_EPS * es / (CONST_RDAIR * Tk);

    /* Diffusivity in m2/s, Liston eq. A-7 */
    Diffusivity = (2.06e-5) * pow(Tk / 273., 1.75);

    // Essery et al. 1999, eq. 6 (m*s/kg)
    F = (Ls / (param.BLOWING_KA * Tk)) * (Ls * Tk / CONST_RDAIR - 1.);
    F += 1. / (Diffusivity * Ros);

    /* grid cell 10 m wind speed = 50th percentile wind */
    /* Wind speed at 2 m above snow was passed to this function. */

    wind10 = Wind * log(10. / ZO) / log((2 + ZO) / ZO);

    /* Check for bare soil case. */
    if (iveg == Nveg) {
        fe = 1500;
        sigma_slope = .0002;
    }
    // sigma_w/uo:
    ratio = (2.44 - (0.43) * lag_one) * sigma_slope;

    sigma_w = wind10 * ratio;
    Uo = wind10;

    /*********** Parameters for roughness above snow. *****************/
    hv = (3. / 2.) * displacement;
    Nd = (4. / 3.) * (roughness / displacement);

    /*******************************************************************/
    /** Begin loop through wind probability function.                  */

    Total = 0.0;
    *TotalTransport = 0.0;
    area = 1. / (double) param.BLOWING_NUMINCS;

    if (snowdepth > 0.0) {
        if (options.BLOWING_SPATIAL_WIND && sigma_w != 0.) {
            for (p = 0; p < param.BLOWING_NUMINCS; p++) {
                SubFlux = lower = upper = 0.0;
                /* Find the limits of integration. */
                if (p == 0) {
                    lower = -9999;
                    upper = Uo + sigma_w * log(2. * (p + 1) * area);
                }
                else if (p > 0 && p < param.BLOWING_NUMINCS / 2) {
                    lower = Uo + sigma_w * log(2. * (p) * area);
                    upper = Uo + sigma_w * log(2. * (p + 1) * area);
                }
                else if (p < (param.BLOWING_NUMINCS - 1) && p >=
                         (double) param.BLOWING_NUMINCS / 2) {
                    lower = Uo - sigma_w * log(2. - 2. * (p * area));
                    upper = Uo - sigma_w * log(2. - 2. * ((p + 1.) * area));
                }
                else if (p == param.BLOWING_NUMINCS - 1) {
                    lower = Uo - sigma_w * log(2. - 2. * (p * area));
                    upper = 9999;
                }

                if (lower > upper) { /* Could happen if lower > Uo*2 */
                    lower = upper;
                    log_err("Error with probability boundaries");
                }


                /* Find expected value of wind speed for the interval. */
                U10 = Uo;
                if (lower >= Uo) {
                    U10 = -0.5 *
                          ((upper +
                            sigma_w) * exp((-1. / sigma_w) * (upper - Uo)) -
                           (lower +
                            sigma_w) *
                           exp((-1. / sigma_w) * (lower - Uo))) / area;
                }
                else if (upper <= Uo) {
                    U10 = 0.5 *
                          ((upper -
                            sigma_w) * exp((1. / sigma_w) * (upper - Uo)) -
                           (lower -
                            sigma_w) *
                           exp((1. / sigma_w) * (lower - Uo))) / area;
                }
                else {
                    log_err("Problem with probability ranges: Increment = %d, "
                            "integration limits = %f - %f", p, upper, lower);
                }

                if (U10 < 0.4) {
                    U10 = .4;
                }

                if (U10 > 25.) {
                    U10 = 25.;
                }
                /*******************************************************************/
                /* Calculate parameters for probability of blowing snow occurence. */
                /* ( Li and Pomeroy 1997) */

                if (snowdepth < hv) {
                    Uveg = U10 / sqrt(1. + 170 * Nd * (hv - snowdepth));
                }
                else {
                    Uveg = U10;
                }


                prob_occurence = get_prob(Tair, Age, SurfaceLiquidWater, Uveg);

                /*******************************************************************/
                /* Calculate threshold shear stress. Send 0 for constant or  */
                /* 1 for variable threshold after Li and Pomeroy (1997)      */

                utshear =
                    get_thresh(Tair, SurfaceLiquidWater, ZO);

                /* Iterate to find actual shear stress during saltation. */

                shear_stress(U10, ZO, &ushear, &Zo_salt, utshear);

                if (ushear > utshear) {
                    SubFlux = CalcSubFlux(EactAir, es, Zrh, AirDens, utshear,
                                          ushear, fe, Tsnow,
                                          Tair, U10, Zo_salt, F, &Transport);
                }
                else {
                    SubFlux = 0.0;
                    Transport = 0.0;
                }

                Total += (1. / (double) param.BLOWING_NUMINCS) * SubFlux *
                         prob_occurence;
                *TotalTransport += (1. / (double) param.BLOWING_NUMINCS) *
                                   Transport * prob_occurence;
            }
        }
        else {
            U10 = Uo;
            /*******************************************************************/
            /* Calculate parameters for probability of blowing snow occurence. */
            /* ( Li and Pomeroy 1997) */

            if (snowdepth < hv) {
                Uveg = U10 / sqrt(1. + 170 * Nd * (hv - snowdepth));
            }
            else {
                Uveg = U10;
            }

            prob_occurence = get_prob(Tair, Age, SurfaceLiquidWater, Uveg);

            /*******************************************************************/
            /* Calculate threshold shear stress. Send 0 for constant or  */
            /* 1 for variable threshold after Li and Pomeroy (1997)      */

            utshear = get_thresh(Tair, SurfaceLiquidWater, ZO);

            /* Iterate to find actual shear stress during saltation. */

            shear_stress(Uo, ZO, &ushear, &Zo_salt, utshear);

            if (ushear > utshear) {
                SubFlux = CalcSubFlux(EactAir, es, Zrh, AirDens, utshear,
                                      ushear, fe, Tsnow,
                                      Tair, Uo, Zo_salt, F, &Transport);
            }
            else {
                SubFlux = 0.0;
                Transport = 0.0;
            }
            Total = SubFlux * prob_occurence;
            *TotalTransport = Transport * prob_occurence;
        }
    }

    if (Total < -.00005) {
        Total = -.00005;
    }

    return Total;
}

/******************************************************************************
 * @brief    Integration is performed by Romberg's method:  Numerical Recipes
 *           in C Section 4.3
 *****************************************************************************/
double
qromb(double (*funcd)(),
      double   es,
      double   Wind,
      double   AirDens,
      double   ZO,
      double   EactAir,
      double   F,
      double   hsalt,
      double   phi_r,
      double   ushear,
      double   Zrh,
      double   a,
      double   b)
{
    extern parameters_struct param;

    double                   ss, dss;
    double                   s[param.BLOWING_MAX_ITER + 1];
    double                   h[param.BLOWING_MAX_ITER + 2];
    int                      j;

    h[1] = 1.0;
    for (j = 1; j <= param.BLOWING_MAX_ITER; j++) {
        s[j] = trapzd(funcd, es, Wind, AirDens, ZO, EactAir, F, hsalt, phi_r,
                      ushear, Zrh, a, b, j);
        if (j >= param.BLOWING_K) {
            polint(&h[j - param.BLOWING_K], &s[j - param.BLOWING_K],
                   param.BLOWING_K, 0.0, &ss, &dss);
            if (fabs(dss) <= DBL_EPSILON * fabs(ss)) {
                return ss;
            }
        }
        h[j + 1] = 0.25 * h[j];
    }
    log_err("Too many steps");
    return 0.; // To avoid warnings.
}

/******************************************************************************
 * @brief    Interpolate a set of N points by fitting a polynomial of degree N-1
 *****************************************************************************/
void
polint(double  xa[],
       double  ya[],
       int     n,
       double  x,
       double *y,
       double *dy)
{
    int     i, m, ns;
    double  den, dif, dift, ho, hp, w;
    double *c = NULL;
    double *d = NULL;

    ns = 1;
    dif = fabs(x - xa[1]);
    c = (double *)malloc((size_t) ((n + 1) * sizeof(double)));
    check_alloc_status(c, "Memory allocation error.");
    d = (double *)malloc((size_t) ((n + 1) * sizeof(double)));
    check_alloc_status(d, "Memory allocation error.");


    for (i = 1; i <= n; i++) {
        if ((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            if ((den = ho - hp) == 0.0) {
                log_err("interpolation error");
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
    free(d);
    free(c);
}

/******************************************************************************
 * @brief    Compute the nth stage of refinement of an extended trapezoidal rule.
 *****************************************************************************/
double
trapzd(double (*funcd)(),
       double   es,
       double   Wind,
       double   AirDens,
       double   ZO,
       double   EactAir,
       double   F,
       double   hsalt,
       double   phi_r,
       double   ushear,
       double   Zrh,
       double   a,
       double   b,
       int      n)
{
    double x, tnm, sum, del;
    int    it, j;

    // TODO: remove use of static variables (see GH #735), for now:
    // make static variables thread safe
    static double s;
    #pragma omp threadprivate(s)

    if (n == 1) {
        return (s = 0.5 *
                    (b -
                     a) *
                    ((*funcd)(a, es, Wind, AirDens, ZO, EactAir, F, hsalt,
                              phi_r, ushear, Zrh) +
                     (*funcd)(b, es, Wind, AirDens, ZO, EactAir, F, hsalt,
                              phi_r, ushear, Zrh)));
    }
    else {
        for (it = 1, j = 1; j < n - 1; j++) {
            it <<= 1;
        }
        tnm = it;
        del = (b - a) / tnm;
        x = a + 0.5 * del;
        for (sum = 0.0, j = 1; j <= it; j++, x += del) {
            sum +=
                (*funcd)(x, es, Wind, AirDens, ZO, EactAir, F, hsalt, phi_r,
                         ushear, Zrh);
        }
        s = 0.5 * (s + (b - a) * sum / tnm);
        return s;
    }
    return 0.; // To avoid warnings.
}

/******************************************************************************
 * @brief    Newton-Raphson method.
 *****************************************************************************/
double
rtnewt(double x1,
       double x2,
       double acc,
       double Ur,
       double Zr)
{
    extern parameters_struct param;

    int                      j;
    double                   df, dx, dxold, f, fh, fl;
    double                   temp, xh, xl, rts;

    get_shear(x1, &fl, &df, Ur, Zr);
    get_shear(x2, &fh, &df, Ur, Zr);

    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        log_err("Root must be bracketed");
    }

    if (fl == 0.0) {
        return x1;
    }
    if (fh == 0.0) {
        return x2;
    }
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2);
    dxold = fabs(x2 - x1);
    dx = dxold;
    get_shear(rts, &f, &df, Ur, Zr);
    for (j = 1; j <= param.BLOWING_MAX_ITER; j++) {
        if ((((rts - xh) * df - f) * ((rts - x1) * df - f) > 0.0) ||
            (fabs(2.0 * f) > fabs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) {
                return rts;
            }
        }
        else {
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) {
                return rts;
            }
        }
        if (fabs(dx) < acc) {
            return rts;
        }
        // if(rts < .025) rts=.025;
        get_shear(rts, &f, &df, Ur, Zr);
        if (f < 0.0) {
            xl = rts;
        }
        else {
            xh = rts;
        }
    }
    log_err("Maximum number of iterations exceeded");
    return 0.; // To avoid warnings.
}

/******************************************************************************
 * @brief    This routine resets the values of all output variables to 0.
 *****************************************************************************/
void
get_shear(double  x,
          double *f,
          double *df,
          double  Ur,
          double  Zr)
{
    *f =
        log(2. * CONST_G * Zr / .12) + log(1 / (x * x)) - CONST_KARMAN * Ur / x;
    *df = CONST_KARMAN * Ur / (x * x) - 2. / x;
}

/******************************************************************************
 * @brief    Calculate the sublimation rate for a given height above the
 *           boundary layer.
 *****************************************************************************/
double
sub_with_height(double z,
                double es,
                double Wind,
                double AirDens,
                double ZO,
                double EactAir,
                double F,
                double hsalt,
                double phi_r,
                double ushear,
                double Zrh)
{
    extern parameters_struct param;

    /* Local variables */
    double                   Rrz, ALPHAz, Mz;
    double                   Rmean, terminal_v, fluctuat_v;
    double                   Vtz, Re, Nu;
    double                   sigz, dMdt;
    double                   temp;
    double                   psi_t, phi_t;


    // Calculate sublimation loss rate (1/s)
    Rrz = 4.6e-5 * pow(z, -.258);
    ALPHAz = 4.08 + 12.6 * z;
    Mz =
        (4. /
         3.) * CONST_PI * CONST_RHOICE * Rrz * Rrz * Rrz *
        (1. + (3. / ALPHAz) + (2. / (ALPHAz * ALPHAz)));

    Rmean = pow((3. * Mz) / (4. * CONST_PI * CONST_RHOICE), 1. / 3.);

    // Pomeroy and Male 1986
    terminal_v = 1.1e7 * pow(Rmean, 1.8);

    // Pomeroy (1988)
    fluctuat_v = 0.005 * pow(Wind, 1.36);

    // Ventilation velocity for turbulent suspension Lee (1975)
    Vtz = terminal_v + 3. * fluctuat_v * cos(CONST_PI / 4.);

    Re = 2. * Rmean * Vtz / param.BLOWING_KIN_VIS;
    Nu = 1.79 + 0.606 * pow(Re, 0.5);

    // LCB: found error in rh calc, 1/20/04, check impact
    sigz = ((EactAir / es) - 1.) * (1.019 + .027 * log(z));

    dMdt = 2 * CONST_PI * Rmean * sigz * Nu / F;
    // sublimation loss rate coefficient (1/s)

    psi_t = dMdt / Mz;

    // Concentration of turbulent suspended snow Kind (1992)

    temp = (0.5 * ushear * ushear) / (Wind * param.BLOWING_SETTLING);
    phi_t = phi_r *
            ((temp +
              1.) *
             pow((z / hsalt),
                 (-1. *
                  param.BLOWING_SETTLING) / (CONST_KARMAN * ushear)) - temp);

    return psi_t * phi_t;
}

/******************************************************************************
 * @brief    Calculate parameters for probability of blowing snow occurence.
 *
 * @note     see Li and Pomeroy 1997
 *****************************************************************************/
double
get_prob(double Tair,
         double Age,
         double SurfaceLiquidWater,
         double U10)
{
    extern option_struct options;

    double               mean_u_occurence;
    double               sigma_occurence;
    double               prob_occurence;

    if (options.BLOWING_CALC_PROB) {
        if (SurfaceLiquidWater < 0.001) {
            mean_u_occurence = 11.2 + 0.365 * Tair + 0.00706 * Tair * Tair +
                               0.9 * log(Age);
            sigma_occurence = 4.3 + 0.145 * Tair + 0.00196 * Tair * Tair;

            prob_occurence = 1. /
                             (1. +
                              exp(sqrt(CONST_PI) *
                                  (mean_u_occurence - U10) / sigma_occurence));
        }
        else {
            mean_u_occurence = 21.;
            sigma_occurence = 7.;

            prob_occurence = 1. /
                             (1. +
                              exp(sqrt(CONST_PI) *
                                  (mean_u_occurence - U10) / sigma_occurence));
        }

        if (prob_occurence < 0.0) {
            prob_occurence = 0.0;
        }
        if (prob_occurence > 1.0) {
            prob_occurence = 1.0;
        }
    }
    else {
        prob_occurence = 1.;
    }

    return prob_occurence;
}

/******************************************************************************
 * @brief    Calculate threshold shear stress.
 *****************************************************************************/
double
get_thresh(double Tair,
           double SurfaceLiquidWater,
           double Zo_salt)
{
    double                   ut10;
    double                   utshear;

    extern parameters_struct param;
    extern option_struct     options;

    if (SurfaceLiquidWater < 0.001) {
        // Threshold wind speed after Li and Pomeroy (1997)
        ut10 = 9.43 + .18 * Tair + .0033 * Tair * Tair;
    }
    else {
        // Threshold wind speed after Li and Pomeroy (1997)
        ut10 = 9.9;
    }

    if (options.BLOWING_VAR_THRESHOLD) {
        // Variable threshold, Li and Pomeroy 1997
        utshear = CONST_KARMAN * ut10 / log(10. / Zo_salt);
    }
    // Constant threshold, i.e. Liston and Sturm
    else {
        utshear = param.BLOWING_UTHRESH;
    }

    return utshear;
}

/******************************************************************************
 * @brief    Iterate to find actual shear stress during saltation.
 *****************************************************************************/
void
shear_stress(double  U10,
             double  ZO,
             double *ushear,
             double *Zo_salt,
             double  utshear)
{
    double umin, umax, xacc;
    double fl, fh, df;

    /* Find min & max shear stress to bracket value. */
    umin = utshear;
    umax = CONST_KARMAN * U10;
    xacc = 0.10 * umin;

    /* Check to see if value is bracketed. */
    get_shear(umin, &fl, &df, U10, 10.);
    get_shear(umax, &fh, &df, U10, 10.);

    if (fl < 0.0 && fh < 0.0) {
        log_err("Solution surpasses upper boundary."
                "fl(%f)=%f, fh(%f)=%f", umin, fl, umax, fh);
    }

    if (fl > 0.0 && fh > 0.0) {
        *Zo_salt = ZO;
        *ushear = CONST_KARMAN * U10 / log(10. / ZO);
    }
    else {
        /* Iterate to find actual shear stress. */
        *ushear = rtnewt(umin, umax, xacc, U10, 10.);
        *Zo_salt = 0.12 * (*ushear) * (*ushear) / (2. * CONST_G);
    }
}

/******************************************************************************
 * @brief    Calculate the sublimation flux.
 *****************************************************************************/
double
CalcSubFlux(double  EactAir,
            double  es,
            double  Zrh,
            double  AirDens,
            double  utshear,
            double  ushear,
            double  fe,
            double  Tsnow,
            double  Tair,
            double  U10,
            double  Zo_salt,
            double  F,
            double *Transport)
{
    extern parameters_struct param;
    extern option_struct     options;

    double                   b, undersat_2;
    double                   SubFlux;
    double                   Qsalt, hsalt;
    double                   phi_s, psi_s;
    double                   T, ztop;
    double                   particle;
    double                   saltation_transport;
    double                   suspension_transport;

    SubFlux = 0.0;
    particle = utshear * 2.8;
    // SBSM:
    if (options.BLOWING_SIMPLE) {
        b = .25;
        if (EactAir >= es) {
            undersat_2 = 0.0;
        }
        else {
            undersat_2 =
                ((EactAir / es) - 1.) * (1. - .027 * log(Zrh) + 0.027 * log(2));
        }
        SubFlux = b * undersat_2 * pow(U10, 5.) / F;
    }
    else {
        // Sublimation flux (kg/m2*s) = mass-concentration * sublimation rate * height
        // for both the saltation layer and the suspension layer

        // Saltation layer is assumed constant with height
        // Maximum saltation transport rate (kg/m*s)
        // Liston and Sturm 1998, eq. 6
        Qsalt = (param.BLOWING_CSALT * AirDens / CONST_G) *
                (utshear / ushear) * (ushear * ushear - utshear * utshear);
        if (options.BLOWING_FETCH) {
            Qsalt *= (1. + (500. / (3. * fe)) * (exp(-3. * fe / 500.) - 1.));
        }

        // Pomeroy and Male (1992)
        hsalt = 0.08436 * pow(ushear, 1.27);

        // Saltation layer mass concentration (kg/m3)
        phi_s = Qsalt / (hsalt * particle);

        T = 0.5 * (ushear * ushear) / (U10 * param.BLOWING_SETTLING);
        ztop = hsalt *
               pow(T / (T + 1.),
                   (CONST_KARMAN * ushear) / (-1. * param.BLOWING_SETTLING));

        if (EactAir >= es) {
            SubFlux = 0.0;
        }
        else {
            // Sublimation loss-rate for the saltation layer (s-1)
            psi_s = sub_with_height(hsalt / 2., es, U10, AirDens, Zo_salt,
                                    EactAir, F, hsalt,
                                    phi_s, ushear, Zrh);

            // Sublimation from the saltation layer in kg/m2*s
            SubFlux = phi_s * psi_s * hsalt;

            // Suspension layer must be integrated
            SubFlux += qromb(sub_with_height, es, U10, AirDens, Zo_salt,
                             EactAir, F, hsalt,
                             phi_s, ushear, Zrh, hsalt, ztop);
        }

        // Transport out of the domain by saltation Qs(fe) (kg/m*s), eq 10 Liston and Sturm
        saltation_transport = Qsalt * (1 - exp(-3. * fe / 500.));

        // Transport in the suspension layer
        suspension_transport = qromb(transport_with_height, es, U10, AirDens,
                                     Zo_salt,
                                     EactAir, F, hsalt, phi_s, ushear, Zrh,
                                     hsalt, ztop);

        // Transport at the downstream edge of the fetch in kg/m*s
        *Transport = (suspension_transport + saltation_transport);
        if (options.BLOWING_FETCH) {
            *Transport /= fe;
        }
    }

    return SubFlux;
}

/******************************************************************************
 * @brief    Calculate the transport rate for a given height above the boundary
 *           layer.
 *****************************************************************************/
double
transport_with_height(double z,
                      double es,
                      double Wind,
                      double AirDens,
                      double ZO,
                      double EactAir,
                      double F,
                      double hsalt,
                      double phi_r,
                      double ushear,
                      double Zrh)
{
    extern parameters_struct param;

    /* Local variables */
    double                   u_z;
    double                   temp;
    double                   phi_t;

    // Find wind speed at current height

    u_z = ushear * log(z / ZO) / CONST_KARMAN;

    // Concentration of turbulent suspended snow Kind (1992)

    temp = (0.5 * ushear * ushear) / (Wind * param.BLOWING_SETTLING);
    phi_t = phi_r *
            ((temp +
              1.) *
             pow((z / hsalt),
                 (-1. *
                  param.BLOWING_SETTLING) / (CONST_KARMAN * ushear)) - temp);

    return u_z * phi_t;
}
