#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

dist_prcp_struct make_dist_prcp(int  nveg,
				int *Nnodes)
/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.  The value of mu changes with
  the intensity of incoming precipitation, and is set in the routine
  dist_prec.

**********************************************************************/
{
  extern option_struct options;

  dist_prcp_struct temp;
  int              i;
  int              Nitems;

#if LAKE_MODEL
  if ( options.LAKES ) Nitems = nveg + 2;
  else Nitems = nveg + 1;
#else // LAKE_MODEL
  Nitems = nveg + 1;
#endif // LAKE_MODEL

  temp.mu     = (double *)calloc(Nitems,sizeof(double));
  for ( i = 0; i < Nitems; i++ ) temp.mu[i] = 1;
  temp.snow   = make_snow_data(Nitems);
  temp.energy = make_energy_bal(Nitems,Nnodes);
  for ( i = 0; i < 2; i++ ) {
    temp.veg_var[i]  = make_veg_var(Nitems-1);
    temp.cell[i]     = make_cell_data(Nitems,options.Nlayer);
  }

  return (temp);

}
