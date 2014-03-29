#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void free_all_vars(all_vars_struct *all_vars, 
		   int              Nveg)
/**********************************************************************
	free_all_vars	Keith Cherkauer		March 1998

  This routine frees all memory allocated down the all_vars 
  data structure.  This include all grid cell specific
  variables (soil, vegetation, energy, snow).

  modifications:
  06-24-98 modified to account for redesign of distributed precipitation
           data structures                                          KAC
  2007-Apr-21 Replaced loop over Nveg to loop over Nitems, so that lake-
	      specific veg tiles could be freed.			TJB
  2009-Jul-31 Removed extra veg tile for lake/wetland.			TJB
  2013-Jul-29 Added freeing of photosynthesis terms.			TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
**********************************************************************/
{
  extern option_struct options;

  int i, j, k, Nitems;

  Nitems = Nveg + 1;

  for(j=0;j<Nitems;j++) {
    free((char *)all_vars[0].cell[j]);
  }
  free((char *)all_vars[0].cell);
  for(j=0;j<Nitems;j++) {
    if (options.CARBON) {
      for ( k = 0 ; k < options.SNOW_BAND ; k++ ) {
        free((char *)all_vars[0].veg_var[j][k].NscaleFactor);
        free((char *)all_vars[0].veg_var[j][k].aPARLayer);
        free((char *)all_vars[0].veg_var[j][k].CiLayer);
        free((char *)all_vars[0].veg_var[j][k].rsLayer);
      }
    }
    free((char *)(*all_vars).veg_var[j]);
  }
  free((char *)(*all_vars).veg_var);
  for(j=0;j<Nitems;j++) {
    free((char *)all_vars[0].energy[j]);
  }
  free((char *)all_vars[0].energy);
  for(i=0;i<Nitems;i++)
    free((char *)all_vars[0].snow[i]);
  free((char *)all_vars[0].snow);

}
