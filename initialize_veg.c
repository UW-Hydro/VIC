#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_veg(veg_var_struct      **veg_var,
		    veg_con_struct       *veg_con,
		    global_param_struct   *gp,
		    int                    Nveg)
/**********************************************************************
  initialize_veg		Dag Lohmann	 January 1996

  This routine initailizes the vegetation variable array.

  Modifications:
  07-13-98 modified to initialize vegetation structure for all 
           defined elevation bands                                 KAC
  11-18-02 modified to get the maximum number of vegetation types
           passed to it.  This allows the maximum number of vegetation
           types to include the wetland vegetation fraction when the 
           lake model is active.                                  LCB

**********************************************************************/
{
  extern option_struct   options;

  int i, j;

  for ( i = 0 ; i < Nveg ; i++) {
    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      veg_var[i][j].Wdew = 0.0;
      veg_var[i][j].throughfall = 0.0;
    }
  }
}
