#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void read_PILPS2c(atmos_data_struct *temp,
                  FILE              *PILPS2c,
                  int               *nrecs,
                  int                dt)
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
  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;
  
  int    i, j, n, rec, maxline = 210;
  int    year, month, day, hour;
  char   str[210],jnkstr[MAXSTRING];
  char   errstr[MAXSTRING];
  
  /** Count Records **/
  n = 0;
  while (fgets(str,maxline,PILPS2c) != '\0') 
    n++;
  printf("nrecs = %d\n",n);
  if(*nrecs>n) {
    fprintf(stderr,"WARNING: SAWD file does not have as many records as defined in the global parameter file, truncating run to %i records.\n",n);
    *nrecs=n;
  }
  if(*nrecs<n) {
    fprintf(stderr,"WARNING: SAWD file has more records then were defined in the global parameter file, run will stop after %i records.\n",*nrecs);
    n = *nrecs;
  }
  
  rewind(PILPS2c);
  
  rec = 0;
  for (i = 0; i < n; i++) {
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
  }

  param_set.SHORTWAVE=TRUE;
  param_set.LONGWAVE=TRUE;
  param_set.PREC=TRUE;
  param_set.AIR_TEMP=TRUE;
  param_set.WIND=TRUE;
  param_set.PRESSURE=TRUE;
  param_set.SPEC_HUMID=TRUE;
}
