 
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>
 
void read_atmosdata(atmos_data_struct *temp,
		    FILE              *snowf,
		    int               *nrecs,
		    int                dt)
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
  char   str[210];
  double junk;

  n = 0;
  while (fgets(str,maxline,snowf) != '\0') n++;
  printf("nrecs = %d\n",n*24/dt);
  if(*nrecs>n*24/dt) {
    fprintf(stderr,"WARNING: daily precipitation file does not have as many records as defined in the global parameter file, truncating run to %i records.\n",n*24/dt);
    *nrecs=n*24/dt;
  }
  if(*nrecs<n) {
    fprintf(stderr,"WARNING: daily precipitation file has more records then were defined in the global parameter file, run will stop after %i records.\n",*nrecs);
    n = *nrecs;
  }

  rewind(snowf);

  rec=0;
  for(i=0;i<n;i++) {
    fgets(str, maxline, snowf);
    for(j=0;j<24/dt;j++) {
      sscanf(str, "%lf %lf %lf", &junk, &temp[rec].tmax, &temp[rec].tmin);
      temp[rec].prec = (double)dt/24. * junk;
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
      rec++;
    }

  }

  param_set.TMAX = param_set.TMIN = param_set.PREC = TRUE;

}
