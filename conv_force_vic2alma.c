#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void conv_force_vic2alma(atmos_data_struct *atmos, int rec, int j, int dt, atmos_data_alma_struct *atmos_alma)
/**********************************************************************
	conv_force_vic2alma	Ted Bohn	March 24, 2005

  This routine converts standard VIC forcing variables into ALMA-
  compliant forcing variables and stores them in the structure atmos_alma.
  
  modifications:
  2005-04-17 Rain and Snow fractions now vary linearly between MIN_RAIN_TEMP
	     and MAX_SNOW_TEMP.						TJB
  2005-04-29 Fixed typo in calculation of Qair.				TJB

**********************************************************************/
{
  extern global_param_struct global_param;

  double MIN_RAIN_TEMP;
  double MAX_SNOW_TEMP;

  MIN_RAIN_TEMP = global_param.MIN_RAIN_TEMP;
  MAX_SNOW_TEMP = global_param.MAX_SNOW_TEMP;

  atmos_alma->SWdown = atmos[rec].shortwave[j];
  atmos_alma->LWdown = atmos[rec].longwave[j];
  atmos_alma->Tair = atmos[rec].air_temp[j] + KELVIN;
  atmos_alma->Qair = 0.622 * atmos[rec].vp[j] / atmos[rec].pressure[j];
  atmos_alma->Psurf = atmos[rec].pressure[j];
  if (atmos[rec].air_temp[j] >= MAX_SNOW_TEMP) {
    atmos_alma->Rainf = atmos[rec].prec[j] / (dt * 3600);
    atmos_alma->Snowf = 0.0;
  }
  else if (atmos[rec].air_temp[j] > MIN_RAIN_TEMP) {
    atmos_alma->Rainf = ((atmos[rec].air_temp[j]-MIN_RAIN_TEMP) / (MAX_SNOW_TEMP-MIN_RAIN_TEMP))
                        * atmos[rec].prec[j] / (dt * 3600);
    atmos_alma->Snowf = atmos[rec].prec[j] / (dt * 3600) - atmos_alma->Rainf;
  }
  else {
    atmos_alma->Rainf = 0.0;
    atmos_alma->Snowf = atmos[rec].prec[j] / (dt * 3600);
  }
  atmos_alma->Wind = atmos[rec].wind[j];

}

