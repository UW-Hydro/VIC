#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

double calc_air_temperature(double *tmax, double *tmin, int hour) {
/**********************************************************************
  calc_air_temperature.c	Keith Cherkauer		March 7, 1998

  This subroutine is based on equations from the NWS snow melt model,
  which estimate air temperature based on minimum and maximum daily
  air temperatures for the 6th, 12th, 18th, and 24th hours of the day.

**********************************************************************/

  double air_temp;

  if(hour == 6)
    air_temp = 0.95 * tmin[1] + 0.05 * tmax[0];
  else if(hour == 12)
    air_temp = 0.40 * tmin[1] + 0.60 * tmax[1];
  else if(hour == 18)
    air_temp = 0.925 * tmax[1] + 0.025 * tmin[1] + 0.05 * tmin[2];
  else if(hour == 24)
    air_temp = 0.33 * tmax[1] + 0.67 * tmin[2];
  else if(hour == 0)
    air_temp = 0.33 * tmax[0] + 0.67 * tmin[1];

  return (air_temp);

}
