#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: correct_precip.c,v 4.2 2000/05/30 22:35:57 vicadmin Exp $";

#define GAUGE_HEIGHT 1.0  /* precipitation gauge height (m) */

void correct_precip(double *gauge_correction, 
		    double  wind,
		    double  wind_h,
		    double  roughness,
		    double  snow_roughness) {
/**********************************************************************
	correct_precip	Keith Cherkauer		May 21, 1997

  This routine corrects preciptation measurements for gauge catch
  deficiencies.  Correction values read from Bras Figure 4.16.

  NOTE: Should locate better reference with fitted equations, or
  at least data with which to fit an equation.

  Modifications:
  05-12-2000 Modified to use a logorithmic wind profile to bring
             observed wind speeds to a height of 2m             KAC
  05-20-2000 Modified to use the WMO correction equations for
             standard NWS shielded 8 inch gauges.  Equations
             presented in: Yang, D., B. E. Goodison, J. R. Metcalfe,
	     V. S. Golubev, R. Bates, T. Pangburn, and C. L. Hanson,
	     "Accuracy of NWS 8" Standard Nonrecording Precipitation
	     Gauge: Results and Application of WMO Intercomparison,"
	     J. Atmos. Oceanic Tech, 15(1), 54-68, 1998.       KAC

**********************************************************************/

  double gauge_wind;

  gauge_wind = wind * (log( ( GAUGE_HEIGHT + roughness ) / roughness ) 
		       / log( wind_h / roughness ));

  gauge_correction[RAIN] = 100. / exp(4.606 - 0.041 
				      * pow ( gauge_wind, 0.69 ) );

  gauge_wind = wind * (log( ( GAUGE_HEIGHT + snow_roughness ) 
			    / snow_roughness ) / 
		       log( wind_h / snow_roughness ));

  gauge_correction[SNOW] = 100. / exp(4.606 - 0.036 
				      * pow ( gauge_wind, 1.75 ) );

}

#undef GAUGE_HEIGHT
