#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
dist_prcp_struct make_dist_prcp(int nveg)
/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.  The value of mu changes with
  the intensity of incoming precipitation, and is set in the routine
  dist_prec.

**********************************************************************/
{
  extern option_struct options;

  dist_prcp_struct temp;

  temp.mu = 1.;
  if(options.DIST_PRCP)
    temp.dist = (prcp_var_struct *) calloc(2,sizeof(prcp_var_struct));
  else
    temp.dist = (prcp_var_struct *) calloc(1,sizeof(prcp_var_struct));

  return (temp);

}
