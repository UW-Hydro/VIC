/******************************************************************************
* @section DESCRIPTION
*
* Calculate values related to the saturated vapor pressure.
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief        This routine computes the saturated vapor pressure
*
* @note         Handbook of Hydrology eqn 4.2.2.
******************************************************************************/
double
svp(double temp)
{
    extern parameters_struct param;

    double                   SVP;

    SVP = param.SVP_A * exp((param.SVP_B * temp) / (param.SVP_C + temp));

    if (temp < 0) {
        SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;
    }

    return (SVP * PA_PER_KPA);
}

/******************************************************************************
* @brief        This routine computes the gradient of d(svp)/dT
*
* @note         Handbook of Hydrology eqn 4.2.3
******************************************************************************/
double
svp_slope(double temp)
{
    extern parameters_struct param;

    return (param.SVP_B * param.SVP_C) / ((param.SVP_C + temp) *
                                          (param.SVP_C + temp)) * svp(temp);
}
