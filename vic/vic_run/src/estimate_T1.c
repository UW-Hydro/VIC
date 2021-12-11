/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine uses Xu Liangs 3-layer energy balance formulation to estimate
 * the temperature between the first and second layers.  Formerly calculated
 * independently in each of the surface energy balance equation routines.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    3-layer energy balance formulation to estimate the temperature
 *           between the first and second layers.
 *****************************************************************************/
double
estimate_T1(double Ts,
            double T1_old,
            double T2,
            double D1,
            double D2,
            double kappa1,
            double kappa2,
            double Cs2,
            double dp,
            double delta_t)
{
    double C1;
    double C2;
    double C3;
    double T1;

    C1 = Cs2 * dp / D2 * (1. - exp(-D2 / dp));
    C2 = -(1. - exp(D1 / dp)) * exp(-D2 / dp);
    C3 = kappa1 / D1 - kappa2 / D1 + kappa2 / D1 * exp(-D1 / dp);

    T1 = (kappa1 / 2. / D1 / D2 * (Ts) + C1 / delta_t * T1_old +
          (2. * C2 - 1. + exp(-D1 / dp)) * kappa2 / 2. / D1 / D2 * T2) /
         (C1 / delta_t + kappa2 / D1 / D2 * C2 + C3 / 2. / D2);

    return(T1);
}
