#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void read_snowmodel(atmos_data_struct *temp,
                    FILE              *snowf,
                    int                nrecs,
                    int                dt,
		    int                file_dt,
		    int                fileskip)
/**********************************************************************
	read_snowmodel	Keith Cherkauer		July 25, 1996

  This routine reads in atmospheric data values from the output
  files for the NWS snow melt model.

  Modifications:
    7-25-96  Routine modified to read data into structure created
	     for time steps of less than 1 day.			KAC

**********************************************************************/
{
  extern param_set_struct param_set;

  int    i, n, rec, maxline = 210;
  int    store_rec;
  int    skip_bytes;
  char   str[210];

  /** locate starting record **/
  skip_bytes = (int)((float)(dt * fileskip)) / (float)file_dt - 1;
  if((dt * fileskip) % (24 / file_dt) > 0) 
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");
  for(i=0;i<skip_bytes;i++) {
    fgets(str, maxline, snowf);
  }

  /** read forcing data **/
  rec = 0;
  while ( !feof(snowf) && (rec < nrecs) ) {
    fscanf(snowf,"%s",str);
    temp[rec].melt = (double)dt/(double)file_dt * atof(str);
    fscanf(snowf,"%s",str);
    temp[rec].prec = (double)dt/(double)file_dt * atof(str);    
    fscanf(snowf,"%s",str);
    temp[rec].air_temp = atof(str);
    fscanf(snowf,"%*s");
    fscanf(snowf,"%*s");
    fscanf(snowf,"%*s");
    fscanf(snowf,"%s",str);
    temp[rec].rainonly = (double)dt/(double)file_dt * atof(str);
    fscanf(snowf,"%*s");
    fscanf(snowf,"%*s"); 
    fscanf(snowf,"%*s");

    fprintf(stderr,"WARNING: Reading Albedo from NWS Snow Model - Check VERSION.  This was not a feature of the original snow model, so make sure your files have albedo in the last column (total of 11 columns).\n");

    fscanf(snowf,"%s",str);
    temp[rec].albedo = atof(str);
    param_set.ALBEDO = TRUE;
    rec ++;

    if(file_dt < dt) {
      /** Time Step in Forcing File Finer than Used by Model: 
	  Skip Records **/
      for(i=0;i<dt/file_dt-1;i++) fgets(str,maxline,snowf);
    }
    else if(file_dt > dt) {
      /** Time step used by model finer than that used in forcing file:
	  Repeat Data Into Extra Columns **/
      store_rec = rec-1;
      for(i=1;i<file_dt/dt;i++) {
	temp[rec].prec     = temp[store_rec].prec;
	temp[rec].melt     = temp[store_rec].melt;
	temp[rec].air_temp = temp[store_rec].air_temp;
	temp[rec].rainonly = temp[store_rec].rainonly;
	rec++;
      }
    }

  }

  if(rec < nrecs) {
    fprintf(stderr,"WARNING: Not enough records in the NWS Snow Model Output forcing file to run the number of records defined in the global file.  Check forcing file time step (%i), and global file.  Unable to run model without enough snow inputs.\n",file_dt);
    exit(0);
  }

  param_set.PREC = TRUE;
  param_set.MELT = TRUE;
  param_set.AIR_TEMP = TRUE;

}
