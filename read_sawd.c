#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void read_sawd(atmos_data_struct *temp,
               FILE              *sawdf,
               int               *nrecs,
               int                dt,
	       int                file_dt,
	       int                fileskip)
/**********************************************************************
	read_sawd	Keith Cherkauer		July 25, 1996

  This routine reads in atmospheric data values from surface
  airways gridded hourly data files.  If time step is less than
  hourly, extra data is skipped.

				Input	Output
				Units	Units
	air_temp	-	C	C
	rel_humid	-	%	%
	tskc		-	%	fraction
	wind		-	m/s	m/s
	pressure	-	Pa	kPa

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;
 
  int    i, n, rec, maxline = 210;
  int    fixcnt, headcnt;
  int    store_rec;
  int    skip_bytes;
  char   str[210],jnkstr[MAXSTRING];

  /** Count Records **/
  n = 0;
  while (fgets(str,maxline,sawdf) != '\0') n++;
  printf("nrecs = %d\n",n);

  rewind(sawdf);

  /** Check for Header, and Skip **/
  fgets(jnkstr,MAXSTRING,sawdf);
  if((jnkstr[0]<48 || jnkstr[0]>57) && jnkstr[0]!=46) {
    fgets(jnkstr,MAXSTRING,sawdf);
    fprintf(stderr,"SAWD ... skipping header\n");
  }

  /** find simulation start record **/
  skip_bytes = (int)((float)(dt * fileskip)) / (float)file_dt - 1;
  if((dt * fileskip) % (24 / file_dt) > 0) 
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");
  for(i=0;i<skip_bytes;i++) {
    fgets(str, maxline, sawdf);
  }

  /** read forcing data **/
  fixcnt = 0;
  rec = 0;
  while ( !feof(sawdf) && (rec < *nrecs) ) {
    fscanf(sawdf,"%*s");
    fscanf(sawdf,"%s",str);
    temp[rec].air_temp = atof(str);
    fscanf(sawdf,"%s",str);
    temp[rec].pressure = atof(str) / 1000.;
    if(temp[rec].pressure <= 0. && rec>0) {
      temp[rec].pressure = temp[rec-1].pressure;
      fixcnt++;
    }
    fscanf(sawdf,"%s",str);
    temp[rec].wind = atof(str);
    fscanf(sawdf,"%s",str);
    temp[rec].rel_humid = atof(str);
    if(temp[rec].rel_humid == 0. && rec>0) {
      temp[rec].rel_humid = temp[rec-1].rel_humid;
      fixcnt++;
    }
    fscanf(sawdf,"%s",str);
    temp[rec].tskc = atof(str) / 100.;
    rec++;

    if(file_dt < dt) {
      /** Time Step in Forcing File Finer than Used by Model: 
	  Skip Records **/
      for(i=0;i<dt/file_dt-1;i++) fgets(str,maxline,sawdf);
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
    fprintf(stderr,"WARNING: Not enough records in the SAWD ASCII forcing file to run the number of records defined in the global file.  Check forcing file time step (%i), and global file.  Number of records being modified to stop model when available data has run out.\n",file_dt);
    *nrecs = rec;
  }

  param_set.WIND = param_set.AIR_TEMP = param_set.PRESSURE = TRUE;
  param_set.TSKC = TRUE;
  param_set.REL_HUMID = TRUE;

}
