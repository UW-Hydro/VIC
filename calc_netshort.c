#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double calc_netshort(double trans, int day, double lat)
/**********************************************************************
	calc_netshort		Dag Lohmann	January 1996

  This routine calculates the net short wave radiation per time step.

**********************************************************************/
{
  double netshort, albedo;

  netshort = 1000.0/SEC_PER_DAY * 
             fltrad((double) 0.0, (double) 0.0, lat, 
             day, trans);
  albedo = 0.2;
  netshort *= (1.0 - albedo);
  return netshort;
}
