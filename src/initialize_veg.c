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
  2013-Jul-25 Added photosynthesis terms.				TJB
  2013-Jul-25 Added AnnualNPP.						TJB
  2014-Apr-25 Added LAI, Wdmax, and albedo.				TJB
  2014-Apr-25 Added vegcover.						TJB

**********************************************************************/
{
  extern option_struct   options;

  int i, j, k;

  for ( i = 0 ; i < Nveg ; i++) {
    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      veg_var[i][j].Wdew = 0.0;
      veg_var[i][j].throughfall = 0.0;
      veg_var[i][j].LAI = 0.0;
      veg_var[i][j].Wdmax = 0.0;
      veg_var[i][j].vegcover = 0.0;
      veg_var[i][j].albedo = 0.0;
      veg_var[i][j].crop_frac = 0.0;
    }
    if (options.CARBON) {
      for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
        veg_var[i][j].aPAR = 0.0;
        veg_var[i][j].Ci = 0.0;
        veg_var[i][j].rc = 0.0;
        veg_var[i][j].GPP = 0.0;
        veg_var[i][j].Rphoto = 0.0;
        veg_var[i][j].Rdark = 0.0;
        veg_var[i][j].Rmaint = 0.0;
        veg_var[i][j].Rgrowth = 0.0;
        veg_var[i][j].Raut = 0.0;
        veg_var[i][j].NPP = 0.0;
        veg_var[i][j].AnnualNPP = 0.0;
        veg_var[i][j].AnnualNPPPrev = 0.0;
        for ( k = 0 ; k < options.Ncanopy ; k++ ) {
          veg_var[i][j].NscaleFactor[k] = 0.0;
          veg_var[i][j].aPARLayer[k] = 0.0;
          veg_var[i][j].CiLayer[k] = 0.0;
          veg_var[i][j].rsLayer[k] = 0.0;
        }
      }
    }
    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      veg_var[i][j].irrig = 0.0;
      veg_var[i][j].irr_apply = FALSE;
    }
  }
}
