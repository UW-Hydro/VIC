#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void read_PILPS2c(atmos_data_struct *temp,
                  FILE              *PILPS2c,
                  int               *nrecs,
                  int                dt,
		  int                file_dt,
		  int                fileskip)
/**********************************************************************
	read_PILPS2c	Dag Lohmann		Feb. 12, 1998

  This routine reads in atmospheric data values from the PILPS2c
  data files.
				Input	Output
				Units	Units
	shortwave               W/m2    W/m2
        longwave                W/m2    W/m2
	prec		-	mm	mm	
	air_temp	-	C	C
	wind		-	m/s	m/s
	pressure	-	mbar	kPa
	spec_humid	-	kg/kg	kg/kg

**********************************************************************/
{
  extern debug_struct debug;
  extern param_set_struct param_set;
  
  int    i, n, rec, maxline = 210;
  int    year, month, day, hour;
  int    store_rec;
  int    skip_bytes;
  char   str[210];
  
  /** locate starting record **/
  skip_bytes = (int)((float)(dt * fileskip)) / (float)file_dt - 1;
  if((dt * fileskip) % (24 / file_dt) > 0) 
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");
  for(i=0;i<skip_bytes;i++) {
    fgets(str, maxline, PILPS2c);
  }

  /** read forcing data **/
  rec = 0;
  while ( !feof(PILPS2c) ) {
    fscanf(PILPS2c,"%d",&year);
    fscanf(PILPS2c,"%d",&month);
    fscanf(PILPS2c,"%d",&day);
    fscanf(PILPS2c,"%d",&hour);
    fscanf(PILPS2c,"%lf",&temp[rec].shortwave);
    fscanf(PILPS2c,"%lf",&temp[rec].longwave);
    fscanf(PILPS2c,"%lf",&temp[rec].prec);
    fscanf(PILPS2c,"%lf",&temp[rec].air_temp);
    fscanf(PILPS2c,"%lf",&temp[rec].wind);
    fscanf(PILPS2c,"%lf",&temp[rec].wind);
    fscanf(PILPS2c,"%lf",&temp[rec].pressure);
    fscanf(PILPS2c,"%lf",&temp[rec].spec_humid);
    temp[rec].pressure =  temp[rec].pressure / 10.; 
    rec++;

    if(file_dt < dt) {
      /** Time Step in Forcing File Finer than Used by Model: 
	  Skip Records **/
      for(i=0;i<dt/file_dt-1;i++) fgets(str,maxline,PILPS2c);
    }
    else if(file_dt > dt) {
      /** Time step used by model finer than that used in forcing file:
	  Repeat Data Into Extra Columns **/
      store_rec = rec;
      for(i=1;i<dt/file_dt;i++) {
	temp[rec].shortwave  = temp[store_rec].shortwave;
	temp[rec].longwave   = temp[store_rec].longwave;
	temp[rec].prec       = temp[store_rec].prec;
	temp[rec].air_temp   = temp[store_rec].air_temp;
	temp[rec].wind       = temp[store_rec].wind;
	temp[rec].pressure   = temp[store_rec].pressure;
	temp[rec].spec_humid = temp[store_rec].spec_humid;
	rec++;
      }
    }
  }

  if(rec < *nrecs) {
    fprintf(stderr,"WARNING: Not enough records in the PILPS2c forcing file to run the number of records defined in the global file.  Check forcing file time step (%i), and global file.  Number of records being modified to stop model when available data has run out.\n",file_dt);
    *nrecs = rec;
  }

  param_set.SHORTWAVE=TRUE;
  param_set.LONGWAVE=TRUE;
  param_set.PREC=TRUE;
  param_set.AIR_TEMP=TRUE;
  param_set.WIND=TRUE;
  param_set.PRESSURE=TRUE;
  param_set.SPEC_HUMID=TRUE;
}
