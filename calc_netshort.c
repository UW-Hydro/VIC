#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double calc_netshort(double  trans, 
		     int     day, 
		     double  lat, 
		     double *day_len)
/**********************************************************************
	calc_netshort		Dag Lohmann	January 1996

  This routine calculates the net short wave radiation per time step.

  09/04/98  Modified to get the length of day in hours from fltrad. KAC

**********************************************************************/
{
  double netshort, albedo;

  netshort = 1000.0/SEC_PER_DAY * 
             fltrad((double) 0.0, (double) 0.0, lat, 
             day, trans, day_len);
  albedo = 0.2;
  netshort *= (1.0 - albedo);
  return netshort;
}
