#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>
 
veg_var_struct **make_veg_var(int veg_type_num)
/**********************************************************************
	make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

  Modifications:
  07-13-98 modified to add structure definitions for all defined 
           elevation bands                                       KAC
  2013-Jul-25 Added photosynthesis terms.				TJB
**********************************************************************/
{
  extern option_struct options;
  
  int              i, j;
  veg_var_struct **temp;

  temp = (veg_var_struct **) calloc(veg_type_num, sizeof(veg_var_struct *));
  for(i=0;i<veg_type_num;i++) {
    temp[i] = (veg_var_struct *) calloc(options.SNOW_BAND, sizeof(veg_var_struct));

    if (options.CARBON) {
      for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
         temp[i][j].NscaleFactor = (double *)calloc(options.Ncanopy,sizeof(double));
         temp[i][j].aPARLayer = (double *)calloc(options.Ncanopy,sizeof(double));
         temp[i][j].CiLayer = (double *)calloc(options.Ncanopy,sizeof(double));
         temp[i][j].rsLayer = (double *)calloc(options.Ncanopy,sizeof(double));
      }
    }

  }

  return temp;
}
