#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void free_dist_prcp(dist_prcp_struct *prcp, int Nveg)
/**********************************************************************
	free_dist_prcp	Keith Cherkauer		March 1998

  This routine frees all memory allocated down the distributed 
  precipitation data structure.  This include all grid cell specific
  variables (soil, vegetation, energy, snow).

**********************************************************************/
{
  extern option_struct options;

  int Ndist;
  int i, j;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;

  for(i=0;i<Ndist;i++) {
    for(j=0;j<=Nveg;j++) {
      free((char *)prcp[0].dist[i].cell[j].layer);
      free((char *)prcp[0].dist[i].energy[j].dz);
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        free((char *)prcp[0].dist[i].energy[j].T);
        if(options.FROZEN_SOIL) free((char *)prcp[0].dist[i].energy[j].ice);
      }
    }
    free((char *)prcp[0].dist[i].cell);
    free((char *)prcp[0].dist[i].veg_var);
    if(options.FULL_ENERGY || options.SNOW_MODEL) {
      free((char *)prcp[0].dist[i].snow);
      free((char *)prcp[0].dist[i].energy);
    }
  }
  free((char *)prcp[0].dist);

}
