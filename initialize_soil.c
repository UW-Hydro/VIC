#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_soil (cell_data_struct **cell, 
                      soil_con_struct   *soil_con,
		      int                veg_num)
/**********************************************************************
	initialize_soil		Keith Cherkauer		July 31, 1996

  This routine initializes the soil variable arrays for each new
  grid cell.

  modifications:
  11-18-02 Modified to initialize wetland soil moisture.          LCB
  2006-Nov-07 Removed LAKE_MODEL option. TJB
  2007-Aug-10 Added features for EXCESS_ICE option.  JCA
**********************************************************************/
{
  extern option_struct options;

  int j, band, index;
  
  for ( j = 0 ; j <= veg_num ; j++)
    for(band=0;band<options.SNOW_BAND;band++) 
      for(index=0;index<options.Nlayer;index++)
	cell[j][band].layer[index].moist = soil_con->init_moist[index];

  if ( options.LAKES ) {
    for(index=0;index<options.Nlayer;index++){
#if EXCESS_ICE
      cell[veg_num][0].layer[index].moist = soil_con->effective_porosity[index]*soil_con->depth[index]*1000.; 
#else
      cell[veg_num][0].layer[index].moist = soil_con->porosity[index]*soil_con->depth[index]*1000.;
#endif
    }
  }

}
