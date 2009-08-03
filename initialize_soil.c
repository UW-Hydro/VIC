#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_soil (cell_data_struct **cell, 
                      soil_con_struct   *soil_con,
                      veg_con_struct    *veg_con,
		      int                veg_num)
/**********************************************************************
	initialize_soil		Keith Cherkauer		July 31, 1996

  This routine initializes the soil variable arrays for each new
  grid cell.

  modifications:
  11-18-02 Modified to initialize wetland soil moisture.          LCB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Aug-10 Added features for EXCESS_ICE option.			JCA
  2009-Mar-16 Modified to use min_liq (minimum allowable liquid water
	      content) instead of resid_moist.  For unfrozen soil,
	      min_liq = resid_moist.					TJB
  2009-Jul-31 Replaced extra lake/wetland veg tile with reference to
	      veg_con[j].LAKE.						TJB
**********************************************************************/
{
  extern option_struct options;

  int j, band, index, frost_area;
  
  for ( j = 0 ; j <= veg_num ; j++) {
    for(band=0;band<options.SNOW_BAND;band++) {
      for(index=0;index<options.Nlayer;index++) {
	cell[j][band].layer[index].moist = soil_con->init_moist[index];
        if (options.LAKES && veg_con[j].LAKE)
#if EXCESS_ICE
          cell[j][0].layer[index].moist = soil_con->effective_porosity[index]*soil_con->depth[index]*1000.;
#else
          cell[j][0].layer[index].moist = soil_con->porosity[index]*soil_con->depth[index]*1000.;
#endif
#if SPATIAL_FROST
        for(frost_area=0;frost_area<FROST_SUBAREAS;frost_area++)
	  cell[j][band].layer[index].min_liq[frost_area] = soil_con->resid_moist[index];
#else
	cell[j][band].layer[index].min_liq = soil_con->resid_moist[index];
#endif
      }
    }
  }

}
