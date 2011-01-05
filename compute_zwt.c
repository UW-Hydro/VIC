#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void compute_zwt(soil_con_struct  *soil_con,
                 cell_data_struct *cell)
/****************************************************************************

  compute_zwt			Ted Bohn		2010-Dec-1
                                                                           
  Function to compute spatial average water table position (zwt).  Water
  table position is measured in cm and is negative below the soil surface.

  modifications:

****************************************************************************/

{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int    i;
  int    lindex;
  double tmp_depth;

  /** Compute each layer's zwt using soil moisture v zwt curve **/
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    i = MAX_ZWTVMOIST-1;
    while (i>=1 && cell->layer[lindex].moist > soil_con->zwtvmoist_moist[lindex][i]) {
      i--;
    }
    if (i == MAX_ZWTVMOIST-1) {
      if (cell->layer[lindex].moist < soil_con->zwtvmoist_moist[lindex][i])
        cell->layer[lindex].zwt = 999; // 999 indicates water table not present in this layer
      else if (cell->layer[lindex].moist == soil_con->zwtvmoist_moist[lindex][i])
        cell->layer[lindex].zwt = soil_con->zwtvmoist_zwt[lindex][i];
    }
    else {
      cell->layer[lindex].zwt = soil_con->zwtvmoist_zwt[lindex][i+1] + (soil_con->zwtvmoist_zwt[lindex][i]-soil_con->zwtvmoist_zwt[lindex][i+1]) * (cell->layer[lindex].moist-soil_con->zwtvmoist_moist[lindex][i+1])/(soil_con->zwtvmoist_moist[lindex][i]-soil_con->zwtvmoist_moist[lindex][i+1]);
    }
  }

  /** Compute total soil column's zwt; this will be the zwt of the lowest layer that isn't completely saturated **/
  tmp_depth = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    tmp_depth += soil_con->depth[lindex];
  }
  lindex = options.Nlayer-1;
  while (lindex>=0 && soil_con->max_moist[lindex]-cell->layer[lindex].moist<=SMALL) {
    tmp_depth -= soil_con->depth[lindex];
    lindex--;
  }
  if (cell->layer[lindex].zwt != 999)
    cell->zwt = cell->layer[lindex].zwt;
  else
    cell->zwt = -tmp_depth*100;

}
