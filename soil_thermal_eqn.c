#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double soil_thermal_eqn(double T, va_list ap) {

  double value;

  double TL;
  double TU;
  double T0;
  double kappa;
  double kappaL;
  double kappaU;
  double Cs;
  double dt;
  double moist;
  double max_moist;
  double bubble;
  double expt;
  double ice0;
  double alpha;
  double beta;
  double gamma;
  double fprime;
  double ice;

  TL = (double) va_arg(ap, double);
  TU = (double) va_arg(ap, double);
  T0 = (double) va_arg(ap, double);
  kappa = (double) va_arg(ap, double);
  kappaL = (double) va_arg(ap, double);
  kappaU = (double) va_arg(ap, double);
  Cs = (double) va_arg(ap, double);
  dt = (double) va_arg(ap, double);
  moist = (double) va_arg(ap, double);
  max_moist = (double) va_arg(ap, double);
  bubble = (double) va_arg(ap, double);
  expt = (double) va_arg(ap, double);
  ice0 = (double) va_arg(ap, double);
  alpha = (double) va_arg(ap, double);
  beta = (double) va_arg(ap, double);
  gamma = (double) va_arg(ap, double);
  fprime = (double) va_arg(ap, double);
  

  if(T<0.) {
    ice = moist - maximum_unfrozen_water(T,max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
  value = T*(alpha*beta*Cs + 4.*kappa*alpha*dt) 
        - beta*dt*(kappaL-kappaU)*(TL-TU) 
        - 2.*alpha*dt*kappa*(TL+TU-gamma*fprime) - alpha*beta*Cs*T0
        - alpha*beta*ice_density*Lf*(ice-ice0);

  return(value);

}
