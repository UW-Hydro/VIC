#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void free_dist_prcp(dist_prcp_struct *prcp, 
		    int               Nveg)
/**********************************************************************
	free_dist_prcp	Keith Cherkauer		March 1998

  This routine frees all memory allocated down the distributed 
  precipitation data structure.  This include all grid cell specific
  variables (soil, vegetation, energy, snow).

  modifications:
  06-24-98 modified to account for redesign of distributed precipitation
           data structures                                          KAC
  2007-Apr-21 Replaced loop over Nveg to loop over Nitems, so that lake-
	      specific veg tiles could be freed.			TJB
  2009-Jul-31 Removed extra veg tile for lake/wetland.			TJB

**********************************************************************/
{
  extern option_struct options;

  int Ndist;
  int i, j, Nitems;

  Ndist = 2;
  Nitems = Nveg + 1;

  for(i=0;i<Ndist;i++) {
    for(j=0;j<Nitems;j++) {
      free((char *)prcp[0].cell[i][j]);
    }
    free((char *)prcp[0].cell[i]);
    for(j=0;j<Nitems;j++) 
      free((char *)(*prcp).veg_var[i][j]);
    free((char *)(*prcp).veg_var[i]);
  }
  for(j=0;j<Nitems;j++) {
    free((char *)prcp[0].energy[j]);
  }
  free((char *)prcp[0].energy);
  for(i=0;i<Nitems;i++)
    free((char *)prcp[0].snow[i]);
  free((char *)prcp[0].snow);
  free((char *)prcp[0].mu);

}
