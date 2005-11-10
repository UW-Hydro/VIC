#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

dist_prcp_struct make_dist_prcp(int  nveg)

/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.  The value of mu changes with
  the intensity of incoming precipitation, and is set in the routine
  dist_prec.

Modifications

  2005-11-09 (Port from 4.1.0) Updated arglist to make_energy_bal() as 
             part of fix for QUICK_FLUX state file compatibility. GCT
**********************************************************************/
{
  extern option_struct options;

  dist_prcp_struct temp;
  int              i;

  temp.mu     = (double *)calloc(nveg+1,sizeof(double));
  for ( i = 0; i < nveg + 1; i++ ) temp.mu[i] = 1;
  temp.snow   = make_snow_data(nveg+1);
  temp.energy = make_energy_bal(nveg+1);
  for(i=0;i<2;i++) {
    temp.veg_var[i]  = make_veg_var(nveg);
    temp.cell[i]     = make_cell_data(nveg+1,options.Nlayer);
  }

  return (temp);

}
