#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void read_rosemount(atmos_data_struct *temp,
		    FILE              *snowf,
		    int               *nrecs,
		    int                dt)
/**********************************************************************
	read_rosemount	Keith Cherkauer		January 7, 1997

  This subroutine reads in hourly meteorolgical data from the
  Rosemount weather file (SHAW model hourly met data format).

  NOTE: Use snow file flag to identify this file.

**********************************************************************/
{
  extern debug_struct debug;
  extern param_set_struct param_set;

  int    i, n, rec, maxline = 210;
  int    FIRST = TRUE;
  int    day, year, hour;
  char   str[210];
  double junk;

  n = 0;
  while (fgets(str,maxline,snowf) != '\0') n++;
  printf("nrecs = %d\n",n);
  if(*nrecs>n) {
    fprintf(stderr,"WARNING: SAWD file does not have as many records as defined in the global parameter file, truncating run to %i records.\n",n);
    *nrecs=n;
  }
  if(*nrecs<n) {
    fprintf(stderr,"WARNING: SAWD file has more records then were defined in the global parameter file, run will stop after %i records.\n",*nrecs);
    n = *nrecs;
  }

  rewind(snowf);

  rec=0;
  for(i=0;i<n;i++) {
    fgets(str, maxline, snowf);
    sscanf(str, "%*i %i",&hour);
    if(!(FIRST && hour!=0) || !FIRST) {
      if(FIRST) FIRST=FALSE;
      sscanf(str, "%i %*i %i %lf %lf %lf %lf %lf %lf", &day, &year, 
          &temp[rec].air_temp, &temp[rec].wind, &temp[rec].rel_humid, 
          &temp[rec].prec, &junk, &temp[rec].shortwave);

      temp[rec].wind *= 1609.347/3600.0; /** convert miles to meter per sec **/
      temp[rec].prec *= 25.4;		/** convert inches to mm **/
      temp[rec].pressure = 92.5;	/** air pressure (kPa) assumed constant
                                            for for Rosemount elevation **/
      if(temp[rec].rel_humid>100) temp[rec].rel_humid=100;
      else if(temp[rec].rel_humid<0) temp[rec].rel_humid=0;
      if(temp[rec].air_temp>150 || temp[rec].air_temp<-150) {
        fprintf(stderr,"ERROR: Invalid Air Temperature %lf\n", 
            temp[rec].air_temp);
        exit(0);
      }
      if(temp[rec].wind<0) temp[rec].wind=0;
      rec++;
    }

  }

  param_set.WIND = param_set.AIR_TEMP = param_set.SHORTWAVE = TRUE;
  param_set.REL_HUMID = param_set.PREC = TRUE;

}
