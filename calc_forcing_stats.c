#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#if OUTPUT_FORCE_STATS

double get_mean(double *, int, double);
double get_stdev(double *, int, double, double);
double get_sum(double *, int, double);
double get_min(double *, int, double);
double get_max(double *, int, double);

void calc_forcing_stats(int                Nrecs,
/* 			int                cell, */
/* 			double             lat, */
/* 			double             lng, */
			atmos_data_struct *atmos) {
/**********************************************************************
  calc_forcing_stats.c     Keith Cherkauer          November 16, 2000

  This routine finds the maximum, minimum and mean values for each 
  data type.  Results are output to stdout for inclusion in screen or 
  log file output.  These statistics are meant only to help the user
  identify possible problems in their input forcing data and are not
  an exhaustive study of that data.

**********************************************************************/

  extern int                 NR, NF;

  int rec, i;
  double *values;

  values = (double *) calloc ( Nrecs, sizeof(double) );

  printf("Variable\tMean\tStd. Dev.\tSum\tMaximum\tMinimum\n");

  /** Air Temperature **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].air_temp[NR];
  printf("air temp (C):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Density **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].density[NR];
  printf("Density (kg/m^3):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Longwave **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].longwave[NR];
  printf("Longwave (W/m^2):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Precipitation **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].prec[NR];
  printf("Precip (mm):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Pressure **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].pressure[NR];
  printf("Pressure (Pa):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Shortwave **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].shortwave[NR];
  printf("Shortwave (W/m^2):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Vapor Pressure **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].vp[NR];
  printf("vp (Pa):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Vapor Pressure Deficit **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].vpd[NR];
  printf("vpd (Pa):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  /** Wind Speed **/
  for ( rec = 0; rec < Nrecs; rec ++ )
    values[rec] = atmos[rec].wind[NR];
  printf("wind speed (m/s):\t%f\t%f\t%f\t%f\t%f\n", 
	 get_mean(values, Nrecs, -999), 
	 get_stdev(values, Nrecs, get_mean(values, Nrecs, -999), -999),
	 get_sum(values, Nrecs, -999), get_max(values, Nrecs, -999), 
	 get_min(values, Nrecs, -999) );

  fflush(stdout);

}

double get_mean(double *values, int N, double NoData) {

  int index;
  int cnt;
  double mean=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      mean+=values[index];
      cnt++;
    }
  }

  if(cnt>0) mean /= (double)cnt;
  else mean = NoData;

  return(mean);

}

double get_stdev(double *values, int N, double mean, double NoData) {

  int index;
  int cnt;
  double stdev=0.0;

  cnt = 0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      stdev+=pow(values[index]-mean,2.0);
      cnt++;
    }
  }

  if(cnt>0) stdev = sqrt((stdev)/(double)(cnt-1));
  else stdev = NoData;

  return(stdev);

}

double get_sum(double *values, int N, double NoData) {

  int index;
  double sum;

  sum=0.0;
  for(index=0;index<N;index++) {
    if(values[index]!=NoData) {
      sum += values[index];
    }
  }
  
  return (sum);

}

double get_min(double *values, int N, double NoData) {

  int index;
  double min;

  min=values[0];
  for(index=1;index<N;index++)
    if(values[index]!=NoData)
      if(min>values[index]) min=values[index];
  
  return (min);

}

double get_max(double *values, int N, double NoData) {

  int index;
  double max;

  max=values[0];
  for(index=1;index<N;index++)
    if(values[index]!=NoData)
      if(max<values[index]) max=values[index];
  
  return (max);

}

#endif
