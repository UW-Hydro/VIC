 
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>
 
void read_atmosdata(atmos_data_struct *temp,
		    FILE              *snowf,
		    int               *nrecs,
		    int                dt,
		    int                file_dt)
/**********************************************************************
	read_atmosdata	Dag Lohmann	

  This routine reads in atmospheric data values from gridded daily
  precipitation station data records.

  Modifications:
    7-25-96  Routine modified to read data into structure created
	     for time steps of less than 1 day.			KAC

**********************************************************************/
{
  extern param_set_struct param_set;

  int    i, j, n, rec, maxline = 210;
  int    store_rec;
  char   str[210];
  double junk;

  n = 0;
  while (fgets(str,maxline,snowf) != '\0') n++;
  printf("nrecs = %d\n",n*24/dt);
  if(n==0)
    nrerror("No data in Daily Precipitation forcing file.  Model stopping...");

  rewind(snowf);

  rec=0;
  while ( !feof(snowf) && (rec < *nrecs) ) {
    fgets(str, maxline, snowf);
    sscanf(str, "%lf %lf %lf", &junk, &temp[rec].tmax, &temp[rec].tmin);
    temp[rec].prec = (double)dt/(double)file_dt * junk;
    if(temp[rec].prec < 0.) {
      fprintf(stderr,"ERROR: Negative precipitation for record #%i.\n",rec);
      exit(0);
    }
    if(temp[rec].tmax < temp[rec].tmin) {
      fprintf(stderr,"WARNING: Maximum temperature (%.3lf) less than minimum (%.3lf) for record #%i.\n\tSwitching values.\n",temp[rec].tmax,temp[rec].tmin,rec);
      junk           = temp[rec].tmax;
      temp[rec].tmax = temp[rec].tmin;
      temp[rec].tmin = junk;
    }
    if(rec>0) {
      if(temp[rec-1].tmax<temp[rec].tmin) temp[rec].tmin=temp[rec-1].tmax;
      if(temp[rec-1].tmin>temp[rec].tmax) temp[rec].tmax=temp[rec-1].tmin;
    }
    rec++;

    if(file_dt < dt) {
      /** Time Step in Forcing File Finer than Used by Model: 
	  Skip Records **/
      for(i=0;i<dt/file_dt-1;i++) fgets(str,maxline,snowf);
    }
    else if(file_dt > dt) {
      /** Time step used by model finer than that used in forcing file:
	  Repeat Data Into Extra Columns **/
      store_rec = rec - 1;
      for(i=1;i<file_dt/dt;i++) {
	temp[rec].tmax = temp[store_rec].tmax;
	temp[rec].tmin = temp[store_rec].tmin;
	temp[rec].prec = temp[store_rec].prec;
	rec++;
      }
    }

  }

  if(rec < *nrecs) {
    fprintf(stderr,"WARNING: Not enough records in the Daily Precipitation forcing file to run the number of records defined in the global file.  Check forcing file time step (%i), and global file.  Number of records being modified to stop model when available data has run out.\n",file_dt);
    *nrecs = rec;
  }

  param_set.TMAX = param_set.TMIN = param_set.PREC = TRUE;

}
