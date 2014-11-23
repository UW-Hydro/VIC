#include <vic_def.h>
#include <vic_run.h>

double
func_atmos_energy_bal(double  Tcanopy,
                      va_list ap)
{
/**********************************************************************
   func_atmos_energy_bal.c      Keith Cherkauer        February 6, 2001

   This routine solves the atmospheric exchange energy balance.

**********************************************************************/

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
    (*SensibleHeat) = atmos_density * Cp * (Tair - Tcanopy) / Ra;

    // compute energy balance error
    // Error = NetRadiation + LatentHeat + (*SensibleHeat);
    Error = InSensible - (*SensibleHeat);

    return (Error);
}
