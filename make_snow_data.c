#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

snow_data_struct **make_snow_data(int nveg)
/**********************************************************************
	make_snow_data	Keith Cherkauer		January 22, 1997

  This routine makes an array of snow cover data structures, one 
  for each vegetation type plus bare soil.

  modifications:
  07-09-98 modified to make te make a two dimensional array which 
           also accounts for a variable number of snow elevation
           bands                                               KAC

**********************************************************************/
{
  extern option_struct options;

  int                i;
  snow_data_struct **temp;

  temp = (snow_data_struct **) calloc(nveg, 
				      sizeof(snow_data_struct *));

  for(i=0;i<nveg;i++) {
    temp[i] = (snow_data_struct *) calloc(options.SNOW_BAND, 
					  sizeof(snow_data_struct));
  }
    
  return temp;
}
