#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void read_sawd_binary(atmos_data_struct *temp,
		      FILE              *sawdf,
		      int               *nrecs,
		      int                dt,
		      int                file_dt,
		      int                fileskip)
/**********************************************************************
  read_sawd_binary	Keith Cherkauer		April 26, 1998

  This routine reads in atmospheric data values from surface
  airways binary gridded hourly data files.

				Input	        Output
				Units	        Units
	air_temp	-	C*100	        C
	pressure	-	Pa*100	        kPa
	wind		-	m/s*1000        m/s
	rel_humid	-	%*100		%
	tskc		-	%*10	        fract

**********************************************************************/
{
  extern debug_struct debug;
  extern param_set_struct param_set;
 
  int       i, n, rec;
  int       fixcnt;
  int       store_rec;
  int       skip_bytes;
  char      tmpmem[20];
  short int values[5];

  /** locate starting record **/
  skip_bytes = (int)((float)(dt * fileskip)) / (float)file_dt - 1;
  if((dt * fileskip) % (24 / file_dt) > 0) 
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");
  fseek(sawdf,skip_bytes*12,SEEK_SET);
  if(feof(sawdf))
    nrerror("No data for the specified time period in the Daily Precipitation forcing file.  Model stopping...");

  /** Read forcing data **/
  fixcnt = 0;
  rec = 0;
  while ( !feof(sawdf) && (rec < *nrecs) ) {
    fread(values,sizeof(short int),5,sawdf);
    temp[rec].air_temp = (double)values[0]/100.;
    temp[rec].pressure = (double)values[1]/100.;
    if(temp[rec].pressure == 0. && rec>0) {
      temp[rec].pressure = temp[rec-1].pressure;
      fixcnt++;
    }
    temp[rec].wind = (double)values[2]/1000.;
    temp[rec].rel_humid = (double)values[3]/100.;
    if(temp[rec].rel_humid == 0. && rec>0) {
      temp[rec].rel_humid = temp[rec-1].rel_humid;
      fixcnt++;
    }
    if(temp[rec].rel_humid > 100.) {
      temp[rec].rel_humid = 100.;
      fixcnt++;
    }
    temp[rec].tskc = (double)values[4] / 1000.;
    rec++;

    if(file_dt < dt) {
      /** Time Step in Forcing File Finer than Used by Model: 
	  Skip Records **/
      for(i=0;i<dt/file_dt-1;i++) fread(tmpmem,sizeof(char),10,sawdf);
    }
    else if(file_dt > dt) {
      /** Time step used by model finer than that used in forcing file:
	  Repeat Data Into Extra Columns **/
      store_rec = rec - 1;
      for(i=1;i<file_dt/dt;i++) {
	temp[rec].air_temp   = temp[store_rec].air_temp;
	temp[rec].wind       = temp[store_rec].wind;
	temp[rec].pressure   = temp[store_rec].pressure;
	temp[rec].rel_humid  = temp[store_rec].rel_humid;
	temp[rec].tskc       = temp[store_rec].tskc;
	rec++;
      }
    }

  }

  if(fixcnt>0) {
    fprintf(stderr,"WARNING: Had to fix %i values in sawd file.\n",fixcnt);
  }

  if(rec < *nrecs) {
    fprintf(stderr,"WARNING: Not enough records in the SAWD Binary forcing file to run the number of records defined in the global file.  Check forcing file time step (%i), and global file.  Number of records being modified to stop model when available data has run out.\n",file_dt);
    *nrecs = rec;
  }

  param_set.WIND = param_set.AIR_TEMP = param_set.PRESSURE = TRUE;
  param_set.TSKC = TRUE;
  param_set.REL_HUMID = TRUE;

}
