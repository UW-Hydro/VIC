#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id: make_cell_data.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

cell_data_struct **make_cell_data(int veg_type_num, int Nlayer)
/**********************************************************************
	make_cell_data	Keith Cherkauer		July 9, 1997

  This subroutine makes an array of type cell, which contains soil
  column variables for a single grid cell.

**********************************************************************/
{
  extern option_struct options;

  int i;
  cell_data_struct **temp;

  temp = (cell_data_struct**) calloc(veg_type_num, 
                                  sizeof(cell_data_struct*));
  for(i=0;i<veg_type_num;i++) {
    temp[i] = (cell_data_struct*) calloc(options.SNOW_BAND, 
					 sizeof(cell_data_struct));
/*     for(j=0;j<options.SNOW_BAND;j++) { */
/*       temp[i][j].layer  */
/* 	= (layer_data_struct*)calloc(Nlayer, */
/* 				     sizeof(layer_data_struct)); */
      
/*     } */
  }
  return temp;
}
