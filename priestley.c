#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double priestley(double tair, 
		 double rad)
/**********************************************************************
	priestley	Dag Lohmann		January 1996

  This routine computes potential evaporation using the Priestly-Taylor
  method.

  REFERENCE:	Kimball, J. S., et. al., "An Improved method for estimating
                surface humidity from daily minimum temperature",
		Agriculture and Froest Meteorology, 85(1997) 87-98.

  double tair    C       Current air temperature
  double rad     kJ/day/m^2   Net radiation at the surface

  returns potential evapotranspiration in kg-m^2/s

**********************************************************************/
{

  double slope, pot;

  slope = svp_slope(tair); /* gradient in kPa/C */
    
  pot = ALPHA_PT * slope/(slope + GAMMA_PT) * 0.9 * rad / LV_PT;
  /** The 0.9 is because ground heat flux (G) is assumed to be 10% of the
      net radiation (Rn), so the radiation term (Rn - G) can be 
      approximated as (0.9 * Rn) **/

  return pot;

}
