#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void free_vegcon(veg_con_struct **veg_con)
/**********************************************************************
  free_vegcon.c	            Keith Cherkauer	  September 25, 1998

  This subroutine frees all components of the veg_con structure.

**********************************************************************/
{
 
  int i;
  
  for(i=0;i<veg_con[0][0].vegetat_type_num;i++) { 
    free((char *)veg_con[0][i].zone_depth);
    free((char *)veg_con[0][i].zone_fract);
  }
  free((char *)veg_con[0]);

}
