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

  Modified:
  6/12/98  

**********************************************************************/

  int    A, B;
  double air_temp;
  double air_temp_A;
  double air_temp_B;
  double a0[] = { 0.330, 0.260, 0.200, 0.140, 0.100, 0.070, 0.050, 0.000, 
                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
  double a1[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.100,
                 0.180, 0.270, 0.380, 0.510, 0.600, 0.700, 0.790, 0.850, 
                 0.900, 0.940, 0.925, 0.890, 0.800, 0.700, 0.600, 0.450 };
  double b1[] = { 0.670, 0.740, 0.800, 0.860, 0.900, 0.930, 0.950, 0.900,
                 0.820, 0.730, 0.620, 0.490, 0.400, 0.300, 0.210, 0.150, 
                 0.100, 0.060, 0.025, 0.000, 0.000, 0.000, 0.000, 0.000 };
  double b2[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
                 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                 0.000, 0.000, 0.050, 0.110, 0.200, 0.300, 0.400, 0.550 };

  air_temp = a0[hour] * tmax[0] + a1[hour] * tmax[1] 
           + b1[hour] * tmin[1] + b2[hour] * tmin[2];

  /** Check for Air Temperature Less Than Tmin at end of Current Day **/
  B = 23;
  air_temp_B = a0[B] * tmax[0] + a1[B] * tmax[1] 
             + b1[B] * tmin[1] + b2[B] * tmin[2];    
  if(hour>19 && air_temp_B<tmin[1]) {
    A = 19;
    air_temp_A = a0[A] * tmax[0] + a1[A] * tmax[1] 
               + b1[A] * tmin[1] + b2[A] * tmin[2];
    air_temp = (hour - B)/(A-B) * (air_temp_A - tmin[1]) + tmin[1];
  }
   
  if(air_temp>tmax[1]) {
    air_temp = tmax[1];
    fprintf(stderr,"WARNING: Daily temperature cyle exceeded maximum temperature\n");
  }
  if(air_temp<tmin[1]) {
    air_temp = tmin[1];
    fprintf(stderr,"WARNING: Daily temperature cyle fell below minimum temperature\n");
  }

  return (air_temp);

}
