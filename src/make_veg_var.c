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
  2013-Jul-25 Added photosynthesis terms.				TJB
  2013-Jul-25 Added looping over water table (zwt) distribution.	TJB
**********************************************************************/
{
  extern option_struct options;
  
  int              i, j, zwtidx;
  veg_var_struct **temp;

  temp = (veg_var_struct **) calloc(veg_type_num, sizeof(veg_var_struct *));
  for(i=0;i<veg_type_num;i++) {
    temp[i] = (veg_var_struct *) calloc(options.SNOW_BAND, sizeof(veg_var_struct));

    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      if (options.CARBON) {
        temp[i][j].NscaleFactor = (double *)calloc(options.Ncanopy,sizeof(double));
        temp[i][j].aPARLayer = (double *)calloc(options.Ncanopy,sizeof(double));
        temp[i][j].CiLayer = (double *)calloc(options.Ncanopy,sizeof(double));
        temp[i][j].rsLayer = (double *)calloc(options.Ncanopy,sizeof(double));
      }
      if (options.DIST_ZWT) {
        temp[i][j].rc_dist_zwt = (double *)calloc(options.Nzwt,sizeof(double));
        if (options.CARBON) {
          temp[i][j].NPPfactor_dist_zwt = (double *)calloc(options.Nzwt,sizeof(double));
          temp[i][j].rsLayer_dist_zwt = (double **)calloc(options.Nzwt,sizeof(double*));
          for(zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {
            temp[i][j].rsLayer_dist_zwt[zwtidx] = (double *)calloc(options.Ncanopy,sizeof(double));
          }
        }
      }
    }

  }

  return temp;
}
