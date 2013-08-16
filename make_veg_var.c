#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id: make_veg_var.c,v 3.1 1999/02/16 18:02:07 vicadmin Exp $";

veg_var_struct **make_veg_var(int veg_type_num)
/**********************************************************************
	make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

  Modifications:
  07-13-98 modified to add structure definitions for all defined 
           elevation bands                                       KAC

**********************************************************************/
{
  extern option_struct options;
  
  int              i;
  veg_var_struct **temp;

  temp = (veg_var_struct **) calloc(veg_type_num, 
				    sizeof(veg_var_struct *));
  for(i=0;i<veg_type_num;i++)
    temp[i] = (veg_var_struct *) calloc(options.SNOW_BAND, 
					sizeof(veg_var_struct));
  return temp;
}
