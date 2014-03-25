#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: calc_Nscale_factors.c,v 5.7 2004/07/07 01:46:14 tbohn Exp $";

void calc_Nscale_factors(char        NscaleFlag,
                         double     *CanopLayerBnd,
                         double      LAItotal,
                         double      lat,
                         double      lng,
                         double      time_zone_lng,
                         dmy_struct  dmy,
                         double     *NscaleFactor)
/**********************************************************************
	calc_Nscale_factors.c	Ted Bohn		Feb 23, 2007

  This subroutine calculates nitrogen scaling factors for all canopy
  layers, following eqns 106 and 107 in Knorr 1997.

  Note: this should only be applied to veg types that have a canopy,
  e.g. trees and shrubs, but not grass or tundra vegetation.

**********************************************************************/
{  
  extern option_struct options;

  dmy_struct dmy_tmp;
  double coszen_noon;
  double k12;
  int cidx;       // canopy layer index

  /* Compute solar zenith angle at local noon */
  dmy_tmp.year = dmy.year;
  dmy_tmp.month = dmy.month;
  dmy_tmp.day = dmy.day;
  dmy_tmp.day_in_year = dmy.day_in_year;
  dmy_tmp.hour = 12;
  coszen_noon = compute_coszen(lat,lng,time_zone_lng,dmy_tmp);
  if (coszen_noon < ZenithMinPar) coszen_noon = ZenithMinPar;

  /* Extinction factor; eqn 119c in Knorr 1997 */
  k12 = 0.5 / coszen_noon;

  /* Condition: LAI > LaiLimit; eqns 107 and 108 in Knorr 1997 */
  for (cidx = 0; cidx < options.Ncanopy; cidx++) {
    if (NscaleFlag && LAItotal > LaiLimit && cidx > 0) {
      NscaleFactor[cidx] = exp( -k12 * CanopLayerBnd[cidx-1] * LAItotal);
      if (NscaleFactor[cidx] < 1e-10) NscaleFactor[cidx] = 1e-10;
    }
    else {
      NscaleFactor[cidx] = 1;
    }
  }

}
