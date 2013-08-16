#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: initialize_veg.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void initialize_veg(veg_var_struct      **veg_var,
		    veg_con_struct       *veg_con,
		    global_param_struct   *gp)
/**********************************************************************
  initialize_veg		Dag Lohmann	 January 1996

  This routine initailizes the vegetation variable array.

  Modifications:
  07-13-98 modified to initialize vegetation structure for all 
           defined elevation bands                                 KAC

**********************************************************************/
{
  extern option_struct   options;

  int i, j;

  for ( i = 0 ; i < veg_con[0].vegetat_type_num ; i++) {
    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      veg_var[i][j].Wdew = 0.0;
      veg_var[i][j].throughfall = 0.0;
    }
  }
}
