#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: initialize_soil.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void initialize_soil (cell_data_struct **cell, 
                      soil_con_struct   *soil_con,
		      int                veg_num)
/**********************************************************************
	initialize_soil		Keith Cherkauer		July 31, 1996

  This routine initializes the soil variable arrays for each new
  grid cell.

**********************************************************************/
{
  extern option_struct options;

  int j, band, index;
  
  for ( j = 0 ; j <= veg_num ; j++)
    for(band=0;band<options.SNOW_BAND;band++) 
      for(index=0;index<options.Nlayer;index++)
	cell[j][band].layer[index].moist = soil_con->init_moist[index];
}
