#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
snow_data_struct *make_snow_data(int nveg)
/**********************************************************************
	make_snow_data	Keith Cherkauer		January 22, 1997

  This routine makes an array of snow cover data structures, one 
  for each vegetation type plus bare soil.

**********************************************************************/
{
  int i;
  snow_data_struct *temp;

  temp = (snow_data_struct*) calloc(nveg, 
                                  sizeof(snow_data_struct));

  /** Initialize all records to unfrozen conditions */
  for(i=0;i<nveg;i++) temp[i].snow = FALSE;

  return temp;
}
