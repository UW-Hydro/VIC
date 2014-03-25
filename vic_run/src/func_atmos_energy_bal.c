#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double func_atmos_energy_bal(double Tcanopy, va_list ap) {
/**********************************************************************
  func_atmos_energy_bal.c      Keith Cherkauer        February 6, 2001

  This routine solves the atmospheric exchange energy balance.

**********************************************************************/

  double  LatentHeat;
  double  NetRadiation;
  double  Ra;
  double  Tair;
  double  atmos_density;
  double  InSensible;

  double *SensibleHeat;
 
  // internal routine variables
  double  Error;

  // extract variables from va_arg
  LatentHeat    = (double)  va_arg(ap, double);
  NetRadiation  = (double)  va_arg(ap, double);
  Ra            = (double)  va_arg(ap, double);
  Tair          = (double)  va_arg(ap, double);
  atmos_density = (double)  va_arg(ap, double);
  InSensible    = (double)  va_arg(ap, double);

  SensibleHeat  = (double *)va_arg(ap, double *);

  // compute sensible heat flux between canopy and atmosphere
  (*SensibleHeat) = atmos_density * Cp * (Tair - Tcanopy) / Ra;

  // compute energy balance error
  //Error = NetRadiation + LatentHeat + (*SensibleHeat);
  Error = InSensible - (*SensibleHeat);

  return ( Error );

}
