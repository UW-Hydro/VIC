#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void correct_precip(double *rain, 
		    double *rainonly, 
		    double  wind,
		    double  wind_h,
		    double  roughness) {
/**********************************************************************
	correct_precip	Keith Cherkauer		May 21, 1997

  This routine corrects preciptation measurements for gauge catch
  deficiencies.  Correction values read from Bras Figure 4.16.

  NOTE: Should locate better reference with fitted equations, or
  at least data with which to fit an equation.

  Modifications:
  05-12-2000 Modified to use a logorithmic wind profile to bring
             observed wind speeds to a height of 2m             KAC

**********************************************************************/

  double snow;
  double Frain,Fsnow;

  wind = log( ( 2. + roughness ) / roughness ) / log( wind_h / roughness );

  if(wind<=8.) {
    Frain = 0.9959 + 0.03296 * wind - 7.748e-6 * exp(wind);
    Fsnow = 0.9987 + 0.09080 * wind - 5.803e-6 * exp(wind);
  }
  else {
    Frain = 0.9959 + 0.03296 * 9. - 7.748e-6 * exp(8.);
    Fsnow = 0.9987 + 0.09080 * 9. - 5.803e-6 * exp(8.);
  }

  snow = *rain - *rainonly;

  snow *= Fsnow;
  *rainonly *= Frain;
  *rain = snow + *rainonly;

}
