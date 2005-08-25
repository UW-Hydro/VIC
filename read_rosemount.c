#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void read_rosemount(atmos_data_struct *temp,
		    FILE              *snowf,
		    int               *nrecs,
		    int                starthour,
		    int                dt,
		    int                file_dt,
		    int                fileskip)
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
  int    store_rec;
  int    skip_bytes;
  char   str[210];
  double junk;

  /** locate starting record **/
  skip_bytes = (int)((float)(dt * fileskip)) / (float)file_dt - 1;
  if((dt * fileskip) % (24 / file_dt) > 0) 
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");
  for(i=0;i<skip_bytes;i++) {
    fgets(str, maxline, snowf);
  }

  /** read forcing data **/
  rec=0;
  while ( !feof(snowf) && (rec < *nrecs) ) {
    fgets(str, maxline, snowf);
    sscanf(str, "%*d %d",&hour);
    if(!(FIRST && hour!=starthour) || !FIRST) {
      if(FIRST) FIRST=FALSE;
      sscanf(str, "%d %*d %d %lf %lf %lf %lf %lf %lf", &day, &year, 
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
	  temp[rec].air_temp  = temp[store_rec].air_temp;
	  temp[rec].wind      = temp[store_rec].wind;
	  temp[rec].rel_humid = temp[store_rec].rel_humid;
	  temp[rec].prec      = temp[store_rec].prec;
	  temp[rec].shortwave = temp[store_rec].shortwave;
	  rec++;
	}
      }
    }
  }

  param_set.WIND = param_set.AIR_TEMP = param_set.SHORTWAVE = TRUE;
  param_set.REL_HUMID = param_set.PREC = TRUE;

}
