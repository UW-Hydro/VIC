 
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void read_snowmodel(atmos_data_struct *temp,
                    FILE              *atmosf,
                    int                nrecs,
                    int                dt)
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

  int    i, j, n, rec, maxline = 210;
  char   str[210];
  double junk;
  double pan_evap;

  n = 0;
  while (fgets(str,maxline,atmosf) != '\0') n++;
  printf("nrecs = %d\n",n);

  rewind(atmosf);

  rec = 0;

  /** Read NWS Snow Melt Output File **/
  for (i = 0; i < n; i++) {
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) temp[rec+j].melt = (double)(dt)/24 * atof(str);
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) temp[rec+j].prec = (double)(dt)/24 * atof(str);    
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) temp[rec+j].air_temp = atof(str);
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) pan_evap = atof(str);
    fscanf(atmosf,"%s",str);
    fscanf(atmosf,"%s",str);
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) temp[rec+j].rainonly = atof(str);
    fscanf(atmosf,"%s",str);
    fscanf(atmosf,"%s",str);       
    fscanf(atmosf,"%s",str);


/** Albedo only calculated using modified snow model code **/
    if(i==0)
      fprintf(stderr,"WARNING: Reading Albedo from NWS Snow Model - Check VERSION\n");
    fscanf(atmosf,"%s",str);
    for(j=0;j<24/dt;j++) temp[rec+j].albedo = atof(str);
    param_set.ALBEDO = TRUE;
/**
    if(i==0)
      fprintf(stderr,"WARNING: Not Reading Albedo from NWS Snow Model - Check VERSION\n");


    for(j=0;j<24/dt;j++) temp[rec+j].albedo = 0.5;
**/
    rec += 24/dt;
  }
  if(rec != nrecs) {
    fprintf(stderr,"\nERROR: number of records from snow melt model (%i) does not match those from precipitation station data (%i).\n",rec,nrecs);
    exit(0);
  }

  param_set.PREC = TRUE;
  param_set.AIR_TEMP = TRUE;

}
