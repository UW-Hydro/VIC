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

**********************************************************************/
{

  atmos_alma->SWdown = atmos[rec].shortwave[j];
  atmos_alma->LWdown = atmos[rec].longwave[j];
  atmos_alma->Tair = atmos[rec].air_temp[j] + KELVIN;
  atmos_alma->Qair = atmos[rec].vp[j] / atmos[rec].pressure[j];
  atmos_alma->Psurf = atmos[rec].pressure[j];
  // Setting upper limit of air temperature for snowfall at 275.65K (2.5 C),
  // based on reference in CLM model to Fig. 1, Plate 3-1, of Snow Hydrology (1956).
  if (atmos[rec].air_temp[j] > 2.5) {
    atmos_alma->Rainf = atmos[rec].prec[j] / (dt * 3600);
    atmos_alma->Snowf = 0.0;
  }
  else {
    atmos_alma->Rainf = 0.0;
    atmos_alma->Snowf = atmos[rec].prec[j] / (dt * 3600);
  }
  atmos_alma->Wind = atmos[rec].wind[j];

}

