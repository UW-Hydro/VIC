#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void initialize_soil (cell_data_struct *cell, 
                      soil_con_struct soil_con,
		      int veg_num)
/**********************************************************************
	initialize_soil		Keith Cherkauer		July 31, 1996

  This routine initializes the soil variable arrays for each new
  grid cell.

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  int i, index;

  for ( i = 0 ; i <= veg_num ; i++) {
    for(index=0;index<options.Nlayer;index++) {
      cell[i].layer[index].moist = soil_con.init_moist[index];
    }
  }
}
