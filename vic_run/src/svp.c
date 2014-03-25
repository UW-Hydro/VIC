#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double svp(double temp)
/**********************************************************************
  This routine computes the saturated vapor pressure using Handbook
  of Hydrology eqn 4.2.2

  Pressure in Pa

**********************************************************************/
{
  double SVP;
  
  SVP = A_SVP * exp((B_SVP * temp)/(C_SVP+temp));

  if(temp<0) SVP *= 1.0 + .00972 * temp + .000042 * temp * temp;

  return (SVP*1000.);
}

double svp_slope(double temp)
/**********************************************************************
  This routine computes the gradient of d(svp)/dT using Handbook
  of Hydrology eqn 4.2.3

  returned value in Pa
**********************************************************************/
{
  return (B_SVP * C_SVP) / ((C_SVP + temp) * (C_SVP + temp)) * svp(temp);
}


 
