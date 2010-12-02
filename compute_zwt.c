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

  int    lindex;
  double tmp_depth;

  /** Compute Zwt assuming moisture above water table has average value of Wcr **/
  tmp_depth = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    tmp_depth += soil_con->depth[lindex];
  }
  lindex = options.Nlayer-1;
  while (lindex >= 0 && soil_con->max_moist[lindex]-cell->layer[lindex].moist <= SMALL) {
    tmp_depth -= soil_con->depth[lindex];
    lindex--;
  }
  if (lindex < 0) {
    cell->zwt = 0;
  }
  else {
    tmp_depth -= (cell->layer[lindex].moist-soil_con->Wcr[lindex])/(soil_con->max_moist[lindex]-soil_con->Wcr[lindex])*soil_con->depth[lindex];
    cell->zwt = -tmp_depth*100; // convert to cm
  }

}
