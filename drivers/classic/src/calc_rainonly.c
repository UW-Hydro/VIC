#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

#define MIN_PREC     1.e-5      /* smallest amount of precipitation that
				   is allowed to fall as snow or rain in
				   a mixed precipitation event */

static char vcid[] = "$Id$";

double calc_rainonly(double air_temp,
		     double prec,
		     double MAX_SNOW_TEMP,
		     double MIN_RAIN_TEMP,
		     double mu) {
/**********************************************************************
  calc_rainonly.c	Keith Cherkauer		March 7, 1998

  Determines from the air temperature what fraction of incoming
  precipitation is frozen and unfrozen (snow and rain).

  Modifications:
  09-22-98 Modified to filter out very small fractions of snow
           or rain in mixed precipitation.  Minimum value MIN_PREC
	   is adjusted to account for the size of mu (minimum
	   is based of fractional precipitation with mu=1, since
	   snow cannot be solved for when mu<1).                  KAC
  10-May-04 Changed test
		if ( MAX_SNOW_TEMP < MIN_RAIN_TEMP )
	    to
		if ( MAX_SNOW_TEMP <= MIN_RAIN_TEMP )
            to avoid possibility of dividing by zero.		TJB
  10-May-04 Changed test
		else if(air_temp > MAX_SNOW_TEMP)
	    to
		else if(air_temp >= MAX_SNOW_TEMP)
	    to fix situation in which, if air_temp = MAX_SNOW_TEMP,
	    rainfall (rainonly) was set to 0 and snowfall was set
	    to 100% of precip, causing function to fail.	TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
           main subroutine, rather than ending the simulation.   GCT/KAC
**********************************************************************/

  double rainonly;

  rainonly = 0.;
  if ( MAX_SNOW_TEMP <= MIN_RAIN_TEMP ) {
    fprintf( stderr, "ERROR: MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
    return (ERROR);
  }
  if(air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
    rainonly = (air_temp - MIN_RAIN_TEMP)
        / (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
  }
  else if(air_temp >= MAX_SNOW_TEMP) {
    rainonly = prec;
  }

  if(mu < 1) rainonly = prec;

  return(rainonly);

}

#undef MIN_RAIN
