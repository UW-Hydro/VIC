#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
#define MAXSTR 512

static char vcid[] = "$Id$";

void initialize_atmos(atmos_data_struct *temp,
                      dmy_struct        *dmy,
                      double             theta_l,
                      double             theta_s,
                      double             phi,
		      double             elevation,
                      double             MAX_SNOW_TEMP,
                      double             MIN_RAIN_TEMP,
		      double            *Tfactor,
                      int                nrecs,
                      int                dt)
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

  Modifications:
  11-18-98  Removed variable array yearly_epot, since yearly potential
            evaporation is no longer used for estimating the dew
            point temperature from daily minimum temperature.   KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;
 
  /*FILE *ftmp;*/

  int     rec;
  int     band;
  int     tmax_hour[3];
  int     tmin_hour[3];
  int     Nyears;
  int     Ndays;
  int     Nhours;
  int     day;
  int     hour;
  double  tmp_tmax, tmp_tmin;
  double  sat_vp;
  double  min_Tfactor;
  double  tmax[3];
  double  tmin[3];
  double *Tair;
  double  Tmin;
  double  priest;
  double  deltat;
  double  shortwave;
  double  tair;
  double  trans_clear;
  double  qdp;
  double *yearly_prec;
  double *day_len_hr;

  double *storerad;
  double  sumrad;

  /*ftmp=open_file("humid_compare","w");*/

  if(!param_set.PREC)
    nrerror("ERROR: Precipitation must be given to the model, check input files\n");

  /***********************
    Initialize Variables 
  ***********************/

  min_Tfactor = Tfactor[0];
  for(band=1;band<options.SNOW_BAND;band++)
    if(Tfactor[band]<min_Tfactor) min_Tfactor = Tfactor[band];
  
  Nyears = dmy[nrecs-1].year - dmy[0].year + 1;
  Ndays  = (nrecs * dt) / 24;
  Nhours = 24 / dt;
  
  yearly_prec = (double*) calloc(Nyears, sizeof(double));
  day_len_hr  = (double*) calloc(Ndays, sizeof(double));

  /************************************************************
    Determine the hours of sunrise and sunset for all records 
  ************************************************************/

  for (day = 0; day < Ndays; day++) {
    get_rise_and_set_hours(&temp[day*Nhours].rise_hour,
			   &temp[day*Nhours].set_hour,
			   theta_l,theta_s,phi,
			   (double)dmy[day*Nhours].day_in_year);
    for(hour=1;hour<Nhours;hour++) {
      temp[day*Nhours+hour].rise_hour = temp[day*Nhours].rise_hour;
      temp[day*Nhours+hour].set_hour  = temp[day*Nhours].set_hour;
    }
  }

  /****************************************************
    Make sure air temperature is set for every record 
  ****************************************************/

  if(param_set.TMAX && param_set.TMIN 
     && !param_set.AIR_TEMP && dt == 24) {
    /** Calculate Average Daily Temeprature from Tmax and Tmin **/
    for (rec = 0; rec < nrecs; rec++) 
      temp[rec].air_temp = (temp[rec].tmax + temp[rec].tmin) / 2.;
  }
  else if(param_set.TMAX && param_set.TMIN 
	  && !param_set.AIR_TEMP && dt<24) {
    /** Estimate Subdaily Temperature Cycle from Tmax and Tmin **/
    Tair = (double *)calloc(Nhours,sizeof(double));
    for (day = 0; day < Ndays; day++) {
      store_max_min_temp(&temp[day*Nhours],tmax,tmax_hour,tmin,
			 tmin_hour,day*Nhours,nrecs,Nhours);
      Tmin = HourlyT(1,tmax_hour,tmax,tmin_hour,tmin,Tair);
      for(hour=0;hour<Nhours;hour++) 
	temp[day*Nhours+hour].air_temp = Tair[hour];
    }
    free((char *)Tair);
  }
  else if(!param_set.TMAX && !param_set.TMIN 
	  && param_set.AIR_TEMP && dt<24) {
    /** Compute Daily Maximum and Minimum Temperatures if Not Provided **/
    for(day=0;day<Ndays;day++) {
      tmp_tmax=tmp_tmin=temp[day*Nhours].air_temp;
      for(hour=1;hour<Nhours;hour++) {
	if(temp[day*Nhours+hour].air_temp>tmp_tmax) 
	  tmp_tmax = temp[day*Nhours+hour].air_temp;
	if(temp[day*Nhours+hour].air_temp<tmp_tmin) 
	  tmp_tmin = temp[day*Nhours+hour].air_temp;
      }
      for(hour=0;hour<Nhours;hour++) {
	temp[day*Nhours+hour].tmax = tmp_tmax;
	temp[day*Nhours+hour].tmin = tmp_tmin;
      }
    }
  }
  else if(!param_set.TMAX && !param_set.TMIN && !param_set.AIR_TEMP) {
    nrerror("ERROR: Must define temperature (either per time step or daily maximum and minimum) in order to run the model!");
  }

  for (rec = 0; rec < nrecs; rec++) {

    /********************************************************
      Determine Rain/Snow fraction of current precipitation
    ********************************************************/

    if(options.SNOW_MODEL) {
      /** Set snow/rain division for highest elevation band **/
	
      if(!options.FULL_ENERGY)
	temp[rec].rainonly 
	  = calc_rainonly(temp[rec].tmin+min_Tfactor, temp[rec].prec,
			  MAX_SNOW_TEMP,MIN_RAIN_TEMP,1.);
      else
	temp[rec].rainonly 
	  = calc_rainonly(temp[rec].air_temp+min_Tfactor, temp[rec].prec,
			  MAX_SNOW_TEMP,MIN_RAIN_TEMP,1.);

      /** Initialize snow melt **/
      temp[rec].melt = 0.0;             /** melt will be calculated **/
    }

    /*****************
      Set wind speed 
    *****************/

    if(!param_set.WIND) temp[rec].wind = 1.5;

    /***********************
      Set bare soil albedo 
    ***********************/

    if(!param_set.ALBEDO) temp[rec].albedo = 0.20;

    /*********************************
      Set atmospheric pressure (kPa) 
    *********************************/

    if(!param_set.PRESSURE) temp[rec].pressure = 95.5;

    /******************************************* 
      Set atmospheric density (kg/m^3) HBH 4.5 
    *******************************************/

    if(!param_set.DENSITY) temp[rec].density = 3.486*temp[rec].pressure/(275.0
        + temp[rec].air_temp);

    /******************************************** 
      Correct precipitation for gage undercatch 
    ********************************************/

    if(options.CORRPREC && temp[rec].prec>0 && param_set.WIND)
      correct_precip(&temp[rec].prec, &temp[rec].rainonly, temp[rec].wind);

  }

  /*****************************************************************
    Calculate Daily and Annual Precipitation and Evaporation
  *****************************************************************/

  storerad = (double *)calloc(Ndays,sizeof(double));

  for(day=0;day<Ndays;day++) {
    rec = day * Nhours;
    if ( day == (Ndays-1) )
      deltat = temp[rec].tmax - temp[rec].tmin;
    else
      deltat = temp[rec].tmax - (temp[rec].tmin 
				  + temp[rec+Nhours].tmin)/2.0;
    deltat = (deltat < 0) ? -deltat : deltat;
    deltat = (deltat==0) ?  (temp[rec].tmax 
			     - temp[rec+Nhours].tmin) : deltat;
    temp[rec].trans = calc_trans(deltat, elevation);
    for(hour=1;hour<Nhours;hour++) temp[rec+hour].trans = temp[rec].trans;
    shortwave = calc_netshort(temp[rec].trans, dmy[rec].day_in_year, 
			      phi,&day_len_hr[day]);

    storerad[day] = shortwave;
    
    tair = (temp[rec].tmax + temp[rec].tmin) / 2.0;
    priest = priestley(tair, shortwave);
    for(hour=0;hour<Nhours;hour++) {
      temp[rec+hour].priest = priest;
      yearly_prec[dmy[rec].year-dmy[0].year] += temp[rec+hour].prec;
    }
  }

  trans_clear = A1_TRANS + A2_TRANS * elevation;

  /*************************************************
    Compute humidity and vapor pressure parameters 
  *************************************************/

  if(!param_set.VPD && !param_set.VP && param_set.REL_HUMID) {
    /** Estimate vapor pressure from relative humidity **/
    for(rec=0;rec<nrecs;rec++) {
      sat_vp = svp(temp[rec].air_temp);
      temp[rec].vpd = sat_vp * (1. - temp[rec].rel_humid / 100.);
      temp[rec].vp = temp[rec].rel_humid / 100. * sat_vp;
      temp[rec].spec_humid = 0.622 * temp[rec].vp / (temp[rec].pressure
						     - 0.378 * temp[rec].vp);
    }
  }
  else if(!param_set.VPD && !param_set.VP && param_set.SPEC_HUMID) {
    /** Estimate vapor pressure from specific humidity **/
    for(rec=0;rec<nrecs;rec++) {
      sat_vp = svp(temp[rec].air_temp);
      temp[rec].vp = temp[rec].spec_humid * temp[rec].pressure / 0.622;
      temp[rec].rel_humid = 100. * temp[rec].vp / sat_vp;
      if(temp[rec].rel_humid>100.) temp[rec].rel_humid=100.;
      temp[rec].vpd = sat_vp * (1. - temp[rec].rel_humid / 100.);
    }
  }
  else if(!param_set.VPD && param_set.VP) {
    /** estimate vapor pressure deficit from vapor pressure **/
    for(rec=0;rec<nrecs;rec++) {
      sat_vp = svp(temp[rec].air_temp);
      temp[rec].vpd = sat_vp - temp[rec].vp;
      if(!param_set.REL_HUMID)
	temp[rec].rel_humid = 100. * temp[rec].vp / sat_vp;
      if(!param_set.SPEC_HUMID)
	temp[rec].spec_humid = 0.622 * temp[rec].vp / (temp[rec].pressure
						       - 0.378 * temp[rec].vp);
    }
  }
  else { 
    /** Estimate vapor pressure from daily tmax and tmin **/
    for(day=0;day<Ndays;day++) {
      tair = (temp[day*Nhours].tmax + temp[day*Nhours].tmin) / 2.0;
      qdp = svp(estimate_dew_point(dmy,yearly_prec,
				   temp[day*Nhours].tmin,
				   temp[day*Nhours].tmax,
				   temp[day*Nhours].priest,
				   day_len_hr[day]*3600.,day*Nhours));
      for(hour=0;hour<Nhours;hour++) {
	rec = day*Nhours+hour;


	/*fprintf(ftmp,"%lf\t%lf\t",temp[rec].rel_humid,temp[rec].vp);*/


	temp[rec].rel_humid = qdp / svp(tair) * 100.;
	temp[rec].vp        = (temp[rec].rel_humid 
				* svp(temp[rec].air_temp)) / 100.;
	temp[rec].vpd       = (svp(temp[rec].air_temp) - temp[rec].vp);


	/*fprintf(ftmp,"%lf\t%lf\n",temp[rec].rel_humid,temp[rec].vp);*/


      }
    }
  } 

  /*******************************
    Estimate radiation forcings 
  *******************************/

  for(rec=0;rec<nrecs;rec++) {
    if(dt<24) {
      if(!param_set.SHORTWAVE || !param_set.LONGWAVE)
	/** Estimate subdaily radiation **/
        calc_long_shortwave(&temp[rec].shortwave, &temp[rec].longwave,
			    &temp[rec].tskc, temp[rec].air_temp,
			    temp[rec].vp, theta_l, theta_s, phi,
			    (double)dmy[rec].day_in_year, 
			    (double)dmy[rec].hour,
			    param_set.SHORTWAVE, param_set.LONGWAVE,
			    param_set.TSKC);
    }
    else {
      if(!param_set.SHORTWAVE || !param_set.LONGWAVE) {
	/** Estimate Daily Radiation **/
	temp[rec].shortwave 
	  = in_shortwave(phi, dmy[rec].day_in_year, temp[rec].trans);
	temp[rec].longwave 
	  = -net_out_longwave(temp[rec].trans, trans_clear, 
			      temp[rec].air_temp, temp[rec].vp*1000.0, 
			      &temp[rec].tskc);
	temp[rec].rad = temp[rec].shortwave + temp[rec].longwave;
      }
    }
  }
  
  for(day=0;day<Ndays;day++) {
    sumrad = 0.;
    for(hour=0;hour<Nhours;hour++) {
      sumrad += temp[day*Nhours+hour].shortwave * 3600.;
    }

    if(storerad[day]>sumrad/3600./24.) {
      fprintf(stderr,"ERROR: Estimated net shortwave radiation exceeds total incoming shortwave radiation.\n");
      fprintf(stdout,"%i\t%lf\t%lf\t%lf\t%lf\n",day,storerad[day],sumrad/3600./24.,sumrad/1000.,sumrad);
      exit(0);
    }
  }

  free((char *)yearly_prec);
  free((char *)day_len_hr);
  
  /*fclose(ftmp);*/

  }
#undef MAXSTR
