#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
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
  for(i=0;i<nveg;i++) {
    temp[i] = (energy_bal_struct*) calloc(options.SNOW_BAND, 
					  sizeof(energy_bal_struct));
    for(j=0;j<options.SNOW_BAND;j++) {
      temp[i][j].frozen = FALSE;
      if(options.FROZEN_SOIL) {
	temp[i][j].T   = (double *)calloc( *Nnodes, sizeof(double));
	temp[i][j].dz  = (double *)calloc( *Nnodes, sizeof(double));
	temp[i][j].ice = (double *)calloc( *Nnodes, sizeof(double));
      }
      else if(options.FULL_ENERGY) {
	*Nnodes        = 3;
	temp[i][j].T   = (double *)calloc(3, sizeof(double));
	temp[i][j].dz  = (double *)calloc(3, sizeof(double));
	temp[i][j].ice = NULL;
      }
      else {
	*Nnodes        = 0;
	temp[i][j].T   = NULL;
	temp[i][j].dz  = (double *)calloc(options.Nlayer, sizeof(double));
	temp[i][j].ice = NULL;
      }
    }
  }

  return temp;
}
