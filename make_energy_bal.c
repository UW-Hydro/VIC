#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
energy_bal_struct *make_energy_bal(int nveg, int *Ulayer, int *Llayer)
/**********************************************************************
	make_energy_bal	Keith Cherkauer		May 26, 1996

  This routine makes an array of frozen soil data structures, one 
  for each vegetation type and bare soil.

**********************************************************************/
{
  extern option_struct options;

  int i;
  energy_bal_struct *temp;

  temp = (energy_bal_struct*) calloc(nveg, 
                                  sizeof(energy_bal_struct));

  /** Initialize all records to unfrozen conditions */
  for(i=0;i<nveg;i++) {
    temp[i].frozen = FALSE;
    if(options.FROZEN_SOIL) {
      temp[i].T = (double *)calloc( *Ulayer+ *Llayer+2, sizeof(double));
      temp[i].dz = (double *)calloc( *Ulayer+ *Llayer+2, sizeof(double));
      temp[i].ice = (double *)calloc( *Ulayer+ *Llayer+2, sizeof(double));
    }
    else if(options.FULL_ENERGY) {
      *Ulayer = 1;
      *Llayer = 1;
      temp[i].T = (double *)calloc(4, sizeof(double));
      temp[i].dz = (double *)calloc(4, sizeof(double));
      temp[i].ice = NULL;
    }
    else {
      temp[i].T = NULL;
      temp[i].dz = (double *)calloc(options.Nlayer, sizeof(double));
      temp[i].ice = NULL;
    }
  }

  return temp;
}
