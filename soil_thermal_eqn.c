#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: soil_thermal_eqn.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

double soil_thermal_eqn(double T, va_list ap) {

  double value;

  double TL;
  double TU;
  double T0;
  double moist;
  double max_moist;
#if QUICK_FS
  double **ufwc_table;
#else
  double bubble;
  double expt;
#endif
  double ice0;
  double gamma;
  double fprime;
  double A;
  double B;
  double C;
  double D;
  double E;
  double ice;

  TL         = (double) va_arg(ap, double);
  TU         = (double) va_arg(ap, double);
  T0         = (double) va_arg(ap, double);
  moist      = (double) va_arg(ap, double);
  max_moist  = (double) va_arg(ap, double);
#if QUICK_FS
  ufwc_table = (double **) va_arg(ap, double **);
#else
  bubble     = (double) va_arg(ap, double);
  expt       = (double) va_arg(ap, double);
#endif
  ice0       = (double) va_arg(ap, double);
  gamma      = (double) va_arg(ap, double);
  fprime     = (double) va_arg(ap, double);
  A          = (double) va_arg(ap, double);
  B          = (double) va_arg(ap, double);
  C          = (double) va_arg(ap, double);
  D          = (double) va_arg(ap, double);
  E          = (double) va_arg(ap, double);

  if(T<0.) {
#if QUICK_FS
    ice = moist - maximum_unfrozen_water_quick(T, max_moist,
					       ufwc_table);
#else
    ice = moist - maximum_unfrozen_water(T,max_moist,bubble,expt);
#endif
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
  value = T*E - A*(TL-TU) - B*(TL+TU-gamma*fprime) - C - D*(ice-ice0);

  return(value);

}
