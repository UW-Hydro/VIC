#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
veg_var_struct *make_veg_var(int veg_type_num)
/**********************************************************************
	make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

**********************************************************************/
{
  veg_var_struct *temp;

  temp = (veg_var_struct*) calloc(veg_type_num, 
                                  sizeof(veg_var_struct));
  return temp;
}
