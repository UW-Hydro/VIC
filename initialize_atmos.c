#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
#define MAXSTR 512

void initialize_atmos(atmos_data_struct *temp,
                      dmy_struct *dmy,
                      double theta_l,
                      double theta_s,
                      double phi,
                      double MAX_SNOW_TEMP,
                      double MIN_RAIN_TEMP,
                      int nrecs,
                      int dt)
/**********************************************************************
	initialize_atmos	Keith Cherkauer		February 3, 1997

  This routine initializes atmospheric variables for hourly time
  step data read in from surface airways data files, or SHAW model
  type input files.  Shortwave and longwave radiation calculations
  are made using equations from Bras' Hydrology.

  WARNING: This subroutine is site specific.  Location parameters
    must be changed before compilation.

  UNITS: mks
	energy - W/m^2

  MODIFIES:
	pressure		(kPa)
	vapor pressure		(kPa)
	density			(kg/m^3)
	surface albedo		(fraction)
	specific humidity	(fraction)
	tmax and tmin		(C)
	longwave radiation	(W/m^2)
	shortwave radiation	(W/m^2)

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;
 
  int    i, j, rec;
  int    date_cnt;
  double jdate;
  double declination, sin_alpha, tau;
  double I0, Ic, m, f, radius;
  double tmax, tmin;
  double i_var;
  double hour;
  double last_tskc=0;
  double sat_vp;

  if(!param_set.PREC)
    nrerror("ERROR: Precipitation must be given to the model, check input files\n");

  if(!param_set.AIR_TEMP && options.FULL_ENERGY)
    nrerror("ERROR: Unable to initialize met data without air temperature.");

  /** Initialize Variables **/
  date_cnt=1;
  for (rec = 0; rec < nrecs; rec++) {

    if(param_set.TMAX && param_set.TMIN && !param_set.AIR_TEMP)
      temp[rec].air_temp = (temp[rec].tmax + temp[rec].tmin) / 2.;

    if(options.FULL_ENERGY) {
      /** Set snow/rain division **/
      temp[rec].rainonly = calc_rainonly(temp[rec].air_temp, temp[rec].prec,
					 MAX_SNOW_TEMP,MIN_RAIN_TEMP);

      /** Initialize snow melt **/
      temp[rec].melt = 0.0;             /** melt will be calculated **/
    }

    /** Set wind speed **/
    if(!param_set.WIND) temp[rec].wind = 1.5;

    /** Set bare soil albedo **/
    if(!param_set.ALBEDO) temp[rec].albedo = 0.20;

    /** Set atmospheric pressure (kPa) **/
    if(!param_set.PRESSURE) temp[rec].pressure = 95.5;

    /** Set atmospheric density (kg/m^3) HBH 4.5 **/
    if(!param_set.DENSITY) temp[rec].density = 3.486*temp[rec].pressure/(275.0
        + temp[rec].air_temp);

    /** Correct precipitation for gage undercatch **/
    if(options.CORRPREC && temp[rec].prec>0 && param_set.WIND)
      correct_precip(&temp[rec].prec, &temp[rec].rainonly, temp[rec].wind);

    if(options.FULL_ENERGY) {
  
      /** Compute humidity and vapor pressure parameters **/
      if(!param_set.SVP) sat_vp = svp(temp[rec].air_temp);
      if(!param_set.VPD && !param_set.VP && param_set.REL_HUMID) {
        temp[rec].vpd = sat_vp * (1. - temp[rec].rel_humid / 100.);
        temp[rec].vp = temp[rec].rel_humid / 100. * sat_vp;
        temp[rec].spec_humid = 0.622 * temp[rec].vp / (temp[rec].pressure
                             - 0.378 * temp[rec].vp);
      }
      else if(!param_set.VPD && !param_set.VP && param_set.SPEC_HUMID) {
        temp[rec].vp = temp[rec].spec_humid * temp[rec].pressure / 0.622;
        temp[rec].rel_humid = 100. * temp[rec].vp / sat_vp;
        if(temp[rec].rel_humid>100.) temp[rec].rel_humid=100.;
        temp[rec].vpd = sat_vp * (1. - temp[rec].rel_humid / 100.);
      }
      else if(!param_set.VPD && param_set.VP) {
        temp[rec].vpd = sat_vp - temp[rec].vp;
        if(!param_set.REL_HUMID)
          temp[rec].rel_humid = 100. * temp[rec].vp / sat_vp;
        if(!param_set.SPEC_HUMID)
          temp[rec].spec_humid = 0.622 * temp[rec].vp / (temp[rec].pressure
                               - 0.378 * temp[rec].vp);
      }
      else
        nrerror("ERROR: VIC does not support this combination and vapor pressure and humidity");
    }

    /**************************************************
      Calculate Minimum and Maximum Daily Temperatures
    **************************************************/

    if(!param_set.TMAX || !param_set.TMIN) {
      if(date_cnt==24/dt) {
        tmax=tmin=temp[rec].air_temp;
        for(j=1;j<24/dt;j++) {
          if(temp[rec-j].air_temp>tmax) tmax=temp[rec-j].air_temp;
          if(temp[rec-j].air_temp<tmin) tmin=temp[rec-j].air_temp;
        }
        for(j=0;j<24/dt;j++) {
          temp[rec-j].tmax = tmax;
          temp[rec-j].tmin = tmin;
        }
        date_cnt=1;
      }
      date_cnt++;
    }

    /**************************************************
      Calculate Incoming Long and Short Wave Radiation (W/m^2)
    **************************************************/

    if(options.FULL_ENERGY) {
      if(!param_set.SHORTWAVE || !param_set.LONGWAVE)
        calc_long_shortwave(&temp[rec].shortwave,
             &temp[rec].longwave,&temp[rec].tskc,temp[rec].air_temp,
             temp[rec].vp,theta_l,theta_s,phi,
             (double)dmy[rec].day_in_year,(double)dmy[rec].hour,
             param_set.SHORTWAVE,param_set.LONGWAVE,param_set.TSKC);
 
    }
  }
}
#undef MAXSTR
