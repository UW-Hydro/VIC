#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void calc_long_shortwave(double *shortwave,
                         double *longwave,
                         double *tskc,
                         double air_temp,
                         double vp,
                         double theta_l,
                         double theta_s,
                         double phi,
                         double jdate,
                         double hour,
                         char   HAVE_SHORTWAVE,
                         char   HAVE_LONGWAVE,
                         char   HAVE_TSKC) {
/**********************************************************************
  calc_long_shortwave.c         Keith Cherkauer         March 7, 1998
 
  This routine computes long and short wave radiation based on
  geographic location, time and cloud cover.  Equations are based on
  Bras 21-47.

      shortwave        incoming shortwave radiation (W/m^2)
      longwave         incoming longwave radiation (W/m^2)
      tskc             fraction of sky covered by clouds (fract)
      air_temp         air temperature (C)
      vp               vapor pressure (kPa)
      theta_l          defined longitude of time zone (degree)
      theta_s          longitude of gird cell (degree)
      phi              latitude of grid cell (degree)
      jdate            day in year of current time step (day)
      hour             hour of current time step (hour)
      HAVE_SHORTWAVE   if TRUE do not calculate shortwave
      HAVE_LONGWAVE    if TRUE do not calculate longwave
      HAVE_TSKC        if TRUE do not calculate cloud coverage
 
  Modifications:
  11-30-98 TSKC is a fractional value from 0 to 1, but was being 
           used as if it was a percent (0 to 100).  Based on Gregs
           testing have fixed the problem.                     KAC

**********************************************************************/
 
  static double last_tskc;

  double declination;
  double i_var;
  double tau;
  double sin_alpha;
  double radius;
  double I0;
  double m;
  double Ic;
  int    sum_exceed=0;
  int    sum_zero=0;

  hour = hour - SOLARTIMEOFFSET;   /* assume shortwave
                              measurements made during previous hour */
  declination = 23.45 * PI / 180.0 * cos(2.0 * PI / 365.0
              * (172.0 - jdate));
  if(fabs(theta_l)==theta_l) i_var = 1.0;
  else i_var = -1.0;
 
  /** Check if sun is east or west of cell longitude **/
  if(((float)hour > (12.+(theta_l-theta_s)*24./360.)
      && (float)hour < (24.+(theta_l-theta_s)*24./360.))
      || ((float)hour < (0. +(theta_l-theta_s)*24./360.)))
    tau = (hour - 12.0 - (i_var/15.0 * (fabs(theta_s)
        - fabs(theta_l)))) * 15.0;
  else tau = (hour + 12.0 - (i_var/15.0 * (fabs(theta_s)
           - fabs(theta_l)))) * 15.0;
  sin_alpha = (sin(declination) * sin(phi*PI/180.0) + cos(declination)
            * cos(phi*PI/180.0) * cos(tau*PI/180.0));
  radius = 1.0+0.017*cos((double)(2*PI/365*(186-jdate)));
  I0 = 1353.0 * sin_alpha / radius / radius;
 
  if(!HAVE_SHORTWAVE || (HAVE_SHORTWAVE && *shortwave<0.0)) *shortwave=0.;
  if(I0>0.0) {
    m = pow((sin_alpha + 0.15*pow((asin(sin_alpha) + 3.885),-1.253)),-1.);
    Ic = I0 * exp(-2.0 * (0.128 - 0.054 * log10(m)) * m);
    if(HAVE_SHORTWAVE) Ic = I0;      /** so far the above eqns appear
                                        to over correct **/
    if(!HAVE_SHORTWAVE && HAVE_TSKC) {
      /** Need to Calculate Shortwave Radiation **/
      *shortwave = (1.0 - 0.65 * ((*tskc) * (*tskc)) / 1.00) * Ic;
      if(*shortwave < 0.0) {
        sum_zero++;
        *shortwave = 0.0;
      }
      if(!HAVE_LONGWAVE)
        *longwave = (1.0 + 0.17 * (*tskc)
		     * (*tskc)) * (0.740+0.0049* (vp)*10.0)
	  * STEFAN_B * (air_temp+KELVIN) * (air_temp+KELVIN) 
	  * (air_temp+KELVIN) * (air_temp+KELVIN) / LWAVE_COR;
    }

    else if(HAVE_SHORTWAVE && !HAVE_LONGWAVE) {
 
      /** Shortwave Measured, Cloud Cover Needed **/
      if(*shortwave < Ic) {
        *tskc = sqrt((1.0 - *shortwave / Ic) / 0.65);
        last_tskc = *tskc;
        *longwave = (1.0 + 0.17 * (*tskc) * (*tskc))
	  * (0.740+0.0049*vp*10.0) * STEFAN_B
	  * (air_temp+KELVIN) * (air_temp+KELVIN) 
	  * (air_temp+KELVIN) * (air_temp+KELVIN)  / LWAVE_COR;
      }
      else {
        /** Measured shortwave exceeds estimated **/
        sum_exceed++;
        *tskc = last_tskc;
        *longwave = (1.0 + 0.17*(*tskc)*(*tskc))
	  * (0.740+0.0049*vp*10.0) * STEFAN_B
	  * (air_temp+KELVIN) * (air_temp+KELVIN) 
	  * (air_temp+KELVIN) * (air_temp+KELVIN)  / LWAVE_COR;
      }
    }
    else {
      nrerror("ERROR: To compute long and shortwave radiation, need TSKC (cloud cover fraction), or measured shortwave");
    }
  }
  else {
    if(!HAVE_SHORTWAVE) *shortwave = 0.0;
    else if(HAVE_SHORTWAVE && !HAVE_TSKC) *tskc = last_tskc;
    if(!HAVE_LONGWAVE)
      *longwave = (1.0 + 0.17*(*tskc)*(*tskc))
	* (0.740+0.0049*vp*10.0) * STEFAN_B
	* (air_temp+KELVIN) * (air_temp+KELVIN) 
	* (air_temp+KELVIN) * (air_temp+KELVIN)  / LWAVE_COR;
  }

}

void get_rise_and_set_hours(int    *rise_hour,
			    int    *set_hour,
			    double  theta_l,
			    double  theta_s,
			    double  phi,
			    double  jdate) {
/**********************************************************************
  get_rise_and_set_hours.c         Keith Cherkauer    September 8, 1998
 
  This routine computes the hours at which the sun rises and sets for
  the current day.  Equations are based on Bras 21-47.

      rise_hour        hour the sun rises above the horizon
      set_hour         hour the sun sets below the horizon
      theta_l          defined longitude of time zone (degree)
      theta_s          longitude of gird cell (degree)
      phi              latitude of grid cell (degree)
      jdate            day in year of current time step (day)
 
**********************************************************************/
 
  double declination;
  double i_var;
  double tau;
  double sin_alpha;
  double radius;
  double I0;
  int    hour;

  rise_hour[0] = -1;
  set_hour[0]  = -1;
  for(hour=0;hour<24;hour++) {
    hour = hour - SOLARTIMEOFFSET;   /* assume shortwave measurements made 
					during previous hour */
    declination = 23.45 * PI / 180.0 * cos(2.0 * PI / 365.0
					   * (172.0 - jdate));
    if(fabs(theta_l)==theta_l) i_var = 1.0;
    else i_var = -1.0;
    
    /** Check if sun is east or west of cell longitude **/
    if(((float)hour > (12.+(theta_l-theta_s)*24./360.)
	&& (float)hour < (24.+(theta_l-theta_s)*24./360.))
       || ((float)hour < (0. +(theta_l-theta_s)*24./360.)))
      tau = (hour - 12.0 - (i_var/15.0 * (fabs(theta_s)
					  - fabs(theta_l)))) * 15.0;
    else tau = (hour + 12.0 - (i_var/15.0 * (fabs(theta_s)
					     - fabs(theta_l)))) * 15.0;
    sin_alpha = (sin(declination) * sin(phi*PI/180.0) + cos(declination)
		 * cos(phi*PI/180.0) * cos(tau*PI/180.0));
    radius = 1.0+0.017*cos((double)(2*PI/365*(186-jdate)));
    I0 = 1353.0 * sin_alpha / radius / radius;

    if(rise_hour[0] < 0 && I0 > 0) rise_hour[0] = hour;
    if(set_hour[0] < 0 && rise_hour[0] > 0 && I0 < 0) set_hour[0] = hour;

  }
}
