#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void rad_and_vpd (atmos_data_struct *atmos, 
		  soil_con_struct    soil_con,
		  int                nrecs, 
		  dmy_struct        *dmy,
		  double           **yearly_prec,
		  double           **yearly_epot)
/**********************************************************************
  rad_and_vpd           Dag Lohmann and Bart Nijssen            1995

  This subroutine calculates net radation and the vapor pressure 
  deficit at the current time step for the current grid cell.  The
  vapor pressure subroutine uses a regression curve and subroutine
  developed by John Kimball 1995 (references noted in aurad.c).  Since
  this method derives questionable shortwave values, a different method
  to calculate the incoming shortwave is used for calculations of 
  net radiation.

  UNITS:	kPa
	
  Modifications:
  5/18/96	Removed first shortwave variable from global usage,
	and added comments.					KAC
  9/4/98  Modified to apply Kimball et. al.'s published regression
          for estimating the daily dew point temperature:
          Kimball, J. S., S. W. Running, R. Nemani, An Improved
          method for estimating surface humidity from Daily 
          Minimum Temperature, Agricultural and Forest Hydrology, 
          Vol. 85, no. 1-2, June 1997, pg 87-98.               KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;

  int    i, no_of_years;
  double deltat, tair, trans_clear, qdp;
  double shortwave;
  double day_len_hr;

  no_of_years = dmy[nrecs-1].year - dmy[0].year + 1;

  yearly_prec[0] = (double*) calloc(no_of_years, sizeof(double));
  yearly_epot[0] = (double*) calloc(no_of_years, sizeof(double));

  /** Compute Priestly Evaporation **/
  if(!options.FULL_ENERGY) {
    for (i = 0; i < nrecs; i++) {
      if ( i == (nrecs-1) )
	deltat = atmos[i].tmax - atmos[i].tmin;
      else
	deltat = atmos[i].tmax - (atmos[i].tmin + atmos[i+1].tmin)/2.0;
      deltat = (deltat < 0) ? -deltat : deltat;
      deltat = (deltat==0) ?  (atmos[i].tmax - atmos[i+1].tmin) : deltat;
      atmos[i].trans = calc_trans(deltat, soil_con.elevation);
      shortwave = calc_netshort(atmos[i].trans, dmy[i].day_in_year, 
				soil_con.lat,&day_len_hr);
      tair = (atmos[i].tmax + atmos[i].tmin) / 2.0;
      atmos[i].priest = priestley(tair, shortwave);
    }
  }
  
  for (i = 0 ; i < nrecs; i++) {
    yearly_prec[0][dmy[i].year-dmy[0].year] += atmos[i].prec;
    yearly_epot[0][dmy[i].year-dmy[0].year] += atmos[i].priest;
  }
  trans_clear = A1_TRANS + A2_TRANS * soil_con.elevation;
  
  for ( i = 0 ; i < nrecs ; i++) {
    tair = (atmos[i].tmax + atmos[i].tmin) / 2.0;
    qdp = svp(estimate_dew_point(dmy,yearly_epot[0],yearly_prec[0],
				 atmos[i].tmin,atmos[i].tmax,atmos[i].priest,
				 day_len_hr*3600.,i));
    if(!param_set.SPEC_HUMID) {
      atmos[i].vpd = (svp(tair) - qdp);
      atmos[i].vp = qdp;
      atmos[i].rel_humid = qdp / svp(tair) * 100.;
    }

    if(!options.FULL_ENERGY) {
      atmos[i].shortwave = /* (1.0 - atmos[i].albedo) *  */
	in_shortwave(soil_con.lat,dmy[i].day_in_year,atmos[i].trans);
      atmos[i].longwave = -net_out_longwave(atmos[i].trans, trans_clear, tair,
					    qdp*1000.0,&atmos[i].tskc);
      atmos[i].rad = atmos[i].shortwave + atmos[i].longwave;
    }
  }
}

double estimate_dew_point(dmy_struct *dmy,
			  double     *yearly_epot,
			  double     *yearly_prec,
			  double      tmin,
			  double      tmax,
			  double      priest,
			  double      day_len,
			  int         day) {
/*******************************************************************
  As described in Kimball et. al. 1997.  See reference above.  KAC  
*******************************************************************/

  double Tdew;
  double EF;

  EF = ((priest / RHO_W) * day_len) / yearly_prec[dmy[day].year-dmy[0].year];
  Tdew = (tmin + KELVIN) * ( -0.127 + 1.121 * (1.003 - 1.444 * EF
					       + 12.312 * EF * EF
					       -32.766 * EF * EF * EF) 
			     + 0.0006 * (tmax - tmin) );

  return(Tdew-KELVIN);
}
