#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double svp(double temp)
/**********************************************************************
  This routine computes the saturated vapor pressure using Handbook
  of Hydrology eqn 4.2.2

  Pressure in kPa

**********************************************************************/
{
  return (A_SVP * exp((B_SVP * temp)/(C_SVP+temp)));
}

double svp_slope(double temp)
/**********************************************************************
  This routine computes the gradient of d(svp)/dT using Handbook
  of Hydrology eqn 4.2.3
**********************************************************************/
{
  return (B_SVP * C_SVP) / pow((C_SVP + temp),2.0) * svp(temp);
}


 
