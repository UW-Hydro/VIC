#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

energy_bal_struct **make_energy_bal(int nveg, int *Nnodes)
/**********************************************************************
	make_energy_bal	Keith Cherkauer		May 26, 1996

  This routine makes an array of frozen soil data structures, one 
  for each vegetation type and bare soil.

**********************************************************************/
{
  extern option_struct options;

  int i, j;
  energy_bal_struct **temp;

  temp = (energy_bal_struct**) calloc(nveg, 
				      sizeof(energy_bal_struct*));

  /** Initialize all records to unfrozen conditions */
  for(i = 0; i < nveg; i++) {
    temp[i] = (energy_bal_struct*) calloc(options.SNOW_BAND, 
					  sizeof(energy_bal_struct));
    for(j = 0; j < options.SNOW_BAND; j++) {
      temp[i][j].frozen = FALSE;
      if(options.QUICK_FLUX) {
	if(options.FULL_ENERGY) {
	  *Nnodes        = 3;
	}
	else {
	  *Nnodes        = 1;
	}
      }
    }
  }

  return temp;
}
