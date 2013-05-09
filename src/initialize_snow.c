#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_snow (snow_data_struct **snow, 
		      int                veg_num,
		      int                cellnum)
/**********************************************************************
	initialize_snow		Keith Cherkauer		January 22, 1997

  This routine initializes the snow variable arrays for each new
  grid cell.

  VARIABLES INITIALIZED:
    snow[i][j].snow;	          TRUE = snow, FALSE = no snow 
    snow[i][j].last_snow;         time steps since last snowfall 
    snow[i][j].snow_canopy;       amount of snow on canopy (m) 
    snow[i][j].swq;               snow water equivalent of the entire pack (m) 
    snow[i][j].surf_water;        liquid water content of the surface 
                                  layer (m) 
    snow[i][j].pack_water;        liquid water content of the snow pack (m) 
    snow[i][j].surf_temp;         depth averaged temperature of the snow pack
                                  surface layer (C) 
    snow[i][j].pack_temp;         depth averaged temperature of the snow pack
                                  (C) 
    snow[i][j].vapor_flux;        depth of water evaporation, sublimation, or 
                                  condensation from snow pack (m) 
    snow[i][j].canopy_vapor_flux; depth of water evaporation, sublimation, or 
                                  condensation from intercepted snow (m) 
    snow[i][j].albedo;            snow surface albedo (fraction) 
    snow[i][j].coldcontent;       cold content of snow pack 
    snow[i][j].mass_error;        snow mass balance error 
    snow[i][j].density;	          snow density (kg/m^3) 
    snow[i][j].depth;	          snow depth (m) 
    snow[i][j].tmp_int_storage;   temporary canopy storage, used in 
                                  snow_canopy 
    snow[i][j].Qnet;              Net energy error in snow model 
    snow[i][j].band_elev;         median elevation of the current snow band 
    snow[i][j].prec_frac;         fracton of precipitation that falls in the 
	  		          current snow band 

  modifications:
  07-09-98 modified to initialize snow variables for each defined
           snow elevation band.                                   KAC
  01-11-99 modified to read new initial snow conditions file format KAC
  04-17-00 removed call for read_initial_snow properties file, the
           file read is now incorporated into a single model state
           file.                                                  KAC
  xx-xx-01 modified to initialize spatially distributed snow variables. KAC
  11-18-02 modified to initalize blowing_snow variable.			LCB
  2006-Oct-16 Removed unused init_snow file.				TJB
  2009-Sep-19 Initializing last_snow to MISSING.			TJB
  2009-Sep-28 Added initialization of some terms that previously had
	      not been initialized.					TJB
  2010-Apr-24 Added initialization of surf_temp_fbcount and fbflag.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
**********************************************************************/
{
  extern option_struct options;
  int i, j;
  int startlayer;

  if(options.FROZEN_SOIL) startlayer=2;
  else startlayer=0;

  for ( i = 0 ; i <= veg_num ; i++ ) {
    for ( j = 0 ; j < options.SNOW_BAND ; j++ ) {
      // State vars
      snow[i][j].albedo            = 0.0;
      snow[i][j].canopy_albedo     = 0.0;
      snow[i][j].coldcontent       = 0.0;
      snow[i][j].coverage          = 0.0;
      snow[i][j].density           = 0.0;
      snow[i][j].depth             = 0.0;
      snow[i][j].last_snow         = MISSING;
      snow[i][j].max_snow_depth    = 0.0;
      snow[i][j].MELTING           = FALSE;
      snow[i][j].pack_temp         = 0.0;
      snow[i][j].pack_water        = 0.0;
      snow[i][j].snow              = FALSE;
      snow[i][j].snow_canopy       = 0.0;
      snow[i][j].store_coverage    = 0.0;
      snow[i][j].store_snow        = FALSE;
      snow[i][j].store_swq         = 0.0;
      snow[i][j].surf_temp         = 0.0;
      snow[i][j].surf_temp_fbflag  = 0;
      snow[i][j].surf_temp_fbcount = 0;
      snow[i][j].surf_water        = 0.0;
      snow[i][j].swq               = 0.0;
      snow[i][j].snow_distrib_slope= 0.0;
      snow[i][j].tmp_int_storage   = 0.0;
      // Fluxes
      snow[i][j].blowing_flux      = 0.0;
      snow[i][j].canopy_vapor_flux = 0.0;
      snow[i][j].mass_error        = 0.0;
      snow[i][j].melt              = 0.0;
      snow[i][j].Qnet              = 0.0;
      snow[i][j].surface_flux      = 0.0;
      snow[i][j].transport         = 0.0;
      snow[i][j].vapor_flux        = 0.0;
    }
  }
}
