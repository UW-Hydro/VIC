#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
void free_dist_prcp(dist_prcp_struct *prcp, int Nveg)
/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.

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
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
        free((char *)prcp[0].dist[i].energy[i].T);
        free((char *)prcp[0].dist[i].energy[i].dz);
        if(options.FROZEN_SOIL) free((char *)prcp[0].dist[i].energy[i].ice);
      }
    }
    free((char *)prcp[0].dist[i].cell);
    free((char *)prcp[0].dist[i].veg_var);
    free((char *)prcp[0].dist[i].snow);
    free((char *)prcp[0].dist[i].energy);
  }
  free((char *)prcp[0].dist);
  free((char *)prcp);

}
