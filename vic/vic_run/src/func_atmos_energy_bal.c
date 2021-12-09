/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine solves the atmospheric exchange energy balance.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This routine solves the atmospheric exchange energy balance.
 *****************************************************************************/
double
func_atmos_energy_bal(double  Tcanopy,
                      va_list ap)
{
    double  Ra;
    double  Tair;
    double  atmos_density;
    double  InSensible;

    double *SensibleHeat;

    // internal routine variables
    double  Error;

    // extract variables from va_arg
    Ra = (double)  va_arg(ap, double);
    Tair = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    InSensible = (double)  va_arg(ap, double);
    SensibleHeat = (double *)va_arg(ap, double *);

    // compute sensible heat flux between canopy and atmosphere
    (*SensibleHeat) = calc_sensible_heat(atmos_density, Tair, Tcanopy, Ra);

    // compute energy balance error
    Error = InSensible - (*SensibleHeat);

    return (Error);
}
