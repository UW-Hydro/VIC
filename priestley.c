#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double priestley(double tair, double rad)
/**********************************************************************
	priestley	Dag Lohmann		January 1996

  This routine computes potential evaporation using the Priestly-Taylor
  method.

  REFERENCE:	Handbook of Hydrology pg 4.16-17

  UNITS: 	slope	kPa/C

**********************************************************************/
{
  double slope, pot;

  slope = svp_slope(tair); 
  /* calculate Priestley-Taylor potential evaporation, the radiation 
     returned by fltrad is in kJ/day/m^2 */
    
  pot = ALPHA_PT * slope/(slope + GAMMA_PT) * 0.9 * rad/ LV_PT *
          SEC_PER_DAY;
  return pot;
}
