#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void initialize_veg (veg_var_struct *veg_var,
		     veg_con_struct *veg_con,
                     global_param_struct gp)
/**********************************************************************
	initialize_veg		Dag Lohmann		January 1996

  This routine initailizes the vegetation variable array.

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;

  int i;

  for ( i = 0 ; i < veg_con[0].vegetat_type_num ; i++) {
    veg_var[i].Wdew = 0.0;
    veg_var[i].throughfall = 0.0;
  }
}
