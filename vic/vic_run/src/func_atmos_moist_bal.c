/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine solves the atmospheric exchange moisture balance.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This routine solves the atmospheric exchange moisture balance.
 *****************************************************************************/
double
func_atmos_moist_bal(double  VPcanopy,
                     va_list ap)
{
    double  InLatentHeat;
    double  Lv;
    double  Ra;
    double  atmos_density;
    double  gamma;
    double  vp; // atmospheric vapor pressure

    double *LatentHeat;

    // internal routine variables
    double  Error;

    // extract variables from va_arg
    InLatentHeat = (double)  va_arg(ap, double);
    Lv = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    gamma = (double)  va_arg(ap, double);
    vp = (double)  va_arg(ap, double);

    LatentHeat = (double *) va_arg(ap, double *);

    // compute sensible heat flux between canopy and atmosphere
    (*LatentHeat) = Lv * calc_sensible_heat(atmos_density, vp, VPcanopy,
                                            gamma * Ra);

    // compute energy balance error
    Error = InLatentHeat - (*LatentHeat);

    return (Error);
}
