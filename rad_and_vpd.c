#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void rad_and_vpd (atmos_data_struct *atmos, soil_con_struct soil_con,
   int nrecs, dmy_struct *dmy)
/**********************************************************************
  	rad_and_vpd

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

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  int    i, no_of_years;
  double deltat, tair, trans_clear, qdp;
  double *yearly_prec, *yearly_epot;
  double shortwave;

  no_of_years = dmy[nrecs-1].year - dmy[0].year + 1;

  yearly_prec = (double*) calloc(no_of_years, sizeof(double));
  yearly_epot = (double*) calloc(no_of_years, sizeof(double));

  /** Compute Priestly Evaporation **/
  if(!options.FULL_ENERGY) {
    for (i = 0; i < nrecs; i++) {
      if ( i == (nrecs-1) )
        deltat = atmos[i].tmax - atmos[i].tmin;
      else
        deltat = atmos[i].tmax - (atmos[i].tmin + atmos[i+1].tmin)/2.0;
      deltat = (deltat < 0) ? -deltat : deltat;
      atmos[i].trans = calc_trans(deltat, soil_con.elevation);
      shortwave = calc_netshort(atmos[i].trans, dmy[i].day_in_year, 
          soil_con.lat);
      tair = (atmos[i].tmax + atmos[i].tmin) / 2.0;
      atmos[i].priest = priestley(tair, shortwave);
    }
  }

  for (i = 0 ; i < nrecs; i++) {
    yearly_prec[dmy[i].year-dmy[0].year] += atmos[i].prec;
    yearly_epot[dmy[i].year-dmy[0].year] += atmos[i].priest;
  }
  trans_clear = A1_TRANS + A2_TRANS * soil_con.elevation;

  for ( i = 0 ; i < nrecs ; i++) {
    tair = (atmos[i].tmax + atmos[i].tmin) / 2.0;
    if (yearly_epot[dmy[i].year-dmy[0].year]/
      yearly_prec[dmy[i].year-dmy[0].year] < HUMID_RATIO) {
      qdp = svp(atmos[i].tmin);
    }
    else { 
      qdp = svp(atmos[i].tmin) * 0.2 *
               (pow((atmos[i].priest/yearly_prec[dmy[i].year-dmy[0].year]),
               -0.1) - 0.83) - 0.016 * tair + 0.867;
    }
    atmos[i].vpd = (svp(tair) - qdp);
    atmos[i].vp = qdp;
    if(!options.FULL_ENERGY) {
      atmos[i].shortwave = (1.0 - atmos[i].albedo) * 
          in_shortwave(soil_con.lat,dmy[i].day_in_year,atmos[i].trans);
      atmos[i].longwave = net_out_longwave(atmos[i].trans, trans_clear, tair,
          qdp*1000.0,&atmos[i].tskc);
      atmos[i].rad = atmos[i].shortwave - atmos[i].longwave;
    }
  }
  free((char *)yearly_prec);
  free((char *)yearly_epot);
}
