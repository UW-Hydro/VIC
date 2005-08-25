#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

snow_data_struct read_initial_snow(FILE *initsnow,
				   int cellnum,
				   int veg,
				   int band)
/******************************************************************************
  read_initial_snow	Keith Cherkauer	      January 11, 1999

  This routine reads parameters for an initial snowpack - if present.

  Parameters Read from File:
  TYPE   NAME                    UNITS   DESCRIPTION
  int    cellnum                 N/A     cell number for current data
  char   snow                    N/A     1 = snow is present, 0 = no snow present
  int    last_snow               (dt)    number of time steps since last snowfall
  double swq                     mm      snow water equivalent
  double surf_temp               C       surface temperature of snowpack
  double density                 kg.m^3  density of snowpack

  Parameters Computed from those in the File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double coverage                fract   fraction of grid cell covered in snow

  Modifications:

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char             ErrStr[MAXSTRING];
  char             tmpstr[MAXSTRING];
  int              index, tmpveg, tmpband;
  snow_data_struct temp; 

  rewind(initsnow);

  fscanf(initsnow, "%s", tmpstr);
  if(tmpstr[0] != '#') {
    index = atoi(tmpstr);
    fscanf(initsnow, "%d %d", &tmpveg, &tmpband);
  }
  else index = tmpveg = tmpband = -999;
  while(index != cellnum && veg != tmpveg && band != tmpband) {
    fgets(tmpstr,MAXSTRING,initsnow);
    fscanf(initsnow, "%s", tmpstr);
    if(tmpstr[0] != '#') {
      index = atoi(tmpstr);
      fscanf(initsnow, "%d %d", &tmpveg, &tmpband);
    }
    else index = tmpveg = tmpband = -999;
  }

  if(!feof(initsnow)) {
    fscanf(initsnow, "%s",  &tmpstr);
    if(tmpstr[0]=='1') temp.snow = TRUE;
    else temp.snow = FALSE;
    fscanf(initsnow, "%d",  &temp.last_snow);
    fscanf(initsnow, "%lf", &temp.swq);
    fscanf(initsnow, "%lf", &temp.surf_temp);
    fscanf(initsnow, "%lf", &temp.density);

    if(temp.swq < MAX_FULL_COVERAGE_SWQ)
      temp.coverage = 1. / MAX_FULL_COVERAGE_SWQ * temp.swq;
    else if(temp.swq > 0.) temp.coverage = 1.0;
    else temp.coverage = 0.0;
  }
  else {
    /** Assume no snow cover **/
    fprintf(stderr,"No snow cover defined for cell number %d.\nModel assuming no snow cover.\n",cellnum);
    temp.snow      = 0;
    temp.last_snow = 0;
    temp.swq       = 0.0;
    temp.surf_temp = 0.0;
    temp.density   = 0.0;
    temp.coverage  = 0.0;
  }
    
  return temp;

} 
