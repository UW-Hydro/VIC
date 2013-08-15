#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double compute_zwt(soil_con_struct  *soil_con,
                   int               lindex,
                   double            moist)
/****************************************************************************

  compute_zwt			Ted Bohn		2010-Dec-1
                                                                           
  Function to compute spatial average water table position (zwt).  Water
  table position is measured in cm and is negative below the soil surface.

  modifications:
  2011-Mar-01 Simplified this function and added wrap_compute_zwt() to
	      call it.							TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
****************************************************************************/

{
  extern option_struct options;

  int    i;
  double zwt;

  /** Compute zwt using soil moisture v zwt curve **/
  i = MAX_ZWTVMOIST-1;
  while (i>=1 && moist > soil_con->zwtvmoist_moist[lindex][i]) {
    i--;
  }
  if (i == MAX_ZWTVMOIST-1) {
    if (moist < soil_con->zwtvmoist_moist[lindex][i])
      zwt = 999; // 999 indicates water table not present in this layer
    else if (moist == soil_con->zwtvmoist_moist[lindex][i])
      zwt = soil_con->zwtvmoist_zwt[lindex][i]; // Just barely enough water for a water table
  }
  else {
    zwt = soil_con->zwtvmoist_zwt[lindex][i+1] + (soil_con->zwtvmoist_zwt[lindex][i]-soil_con->zwtvmoist_zwt[lindex][i+1]) * (moist-soil_con->zwtvmoist_moist[lindex][i+1])/(soil_con->zwtvmoist_moist[lindex][i]-soil_con->zwtvmoist_moist[lindex][i+1]); // interpolate to find water table level
  }

  return(zwt);

}


void wrap_compute_zwt(soil_con_struct  *soil_con,
                      cell_data_struct *cell)
/****************************************************************************

  wrap_compute_zwt			Ted Bohn		2011-Mar-1
                                                                           
  Function to compute spatial average water table position (zwt) for
  individual layers as well as various total-column versions of zwt.  Water
  table position is measured in cm and is negative below the soil surface.

  modifications:
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Feb-07 Removed OUT_ZWT2 and OUT_ZWTL; renamed OUT_ZWT3 to
	      OUT_ZWT_LUMPED.						TJB
  2013-Jul-25 Added looping over water table (zwt) distribution.	TJB
****************************************************************************/

{
  extern option_struct options;

  int    i;
  int    lindex;
  int    zwtidx;
  double total_depth;
  double tmp_depth;
  double tmp_moist;

  /** Compute total soil column depth **/
  total_depth = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    total_depth += soil_con->depth[lindex];
  }

  /** Compute each layer's zwt using soil moisture v zwt curve **/
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    cell->layer[lindex].zwt = compute_zwt(soil_con, lindex, cell->layer[lindex].moist);
  }
  if (cell->layer[options.Nlayer-1].zwt == 999) cell->layer[options.Nlayer-1].zwt = -total_depth*100; // in cm

  /** Compute total soil column's zwt; this will be the zwt of the lowest layer that isn't completely saturated **/
  lindex = options.Nlayer-1;
  tmp_depth = total_depth;
  while (lindex>=0 && soil_con->max_moist[lindex]-cell->layer[lindex].moist<=SMALL) {
    tmp_depth -= soil_con->depth[lindex];
    lindex--;
  }
  if (lindex < 0) cell->zwt = 0;
  else if (lindex < options.Nlayer-1) {
    if (cell->layer[lindex].zwt != 999)
      cell->zwt = cell->layer[lindex].zwt;
    else
      cell->zwt = -tmp_depth*100;
  }
  else
    cell->zwt = cell->layer[lindex].zwt;

  /** Compute total soil column's zwt_lumped; this will be the zwt of all N layers lumped together. **/
  tmp_moist = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    tmp_moist += cell->layer[lindex].moist;
  }
  cell->zwt_lumped = compute_zwt(soil_con, options.Nlayer+1, tmp_moist);
  if (cell->zwt_lumped == 999) cell->zwt_lumped = -total_depth*100; // in cm;

  if (options.DIST_ZWT) {

    /** Compute distributed water table depth **/
    for (zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {

      /** Compute each layer's zwt using soil moisture v zwt curve **/
      for (lindex=0; lindex<options.Nlayer; lindex++) {
        cell->layer[lindex].zwt_dist_zwt[zwtidx] = compute_zwt(soil_con, lindex, cell->layer[lindex].moist_dist_zwt[zwtidx]);
      }
      if (cell->layer[options.Nlayer-1].zwt_dist_zwt[zwtidx] == 999) cell->layer[options.Nlayer-1].zwt_dist_zwt[zwtidx] = -total_depth*100; // in cm

      /** Compute total soil column's zwt; this will be the zwt of the lowest layer that isn't completely saturated **/
      lindex = options.Nlayer-1;
      tmp_depth = total_depth;
      while (lindex>=0 && soil_con->max_moist[lindex]-cell->layer[lindex].moist_dist_zwt[zwtidx]<=SMALL) {
        tmp_depth -= soil_con->depth[lindex];
        lindex--;
      }
      if (lindex < 0) cell->zwt_dist_zwt[zwtidx] = 0;
      else if (lindex < options.Nlayer-1) {
        if (cell->layer[lindex].zwt_dist_zwt[zwtidx] != 999)
          cell->zwt_dist_zwt[zwtidx] = cell->layer[lindex].zwt_dist_zwt[zwtidx];
        else
          cell->zwt_dist_zwt[zwtidx] = -tmp_depth*100;
      }
      else
        cell->zwt_dist_zwt[zwtidx] = cell->layer[lindex].zwt_dist_zwt[zwtidx];

      /** Compute total soil column's zwt_lumped; this will be the zwt of all N layers lumped together. **/
      tmp_moist = 0;
      for (lindex=0; lindex<options.Nlayer; lindex++) {
        tmp_moist += cell->layer[lindex].moist_dist_zwt[zwtidx];
      }
      cell->zwt_lumped_dist_zwt[zwtidx] = compute_zwt(soil_con, options.Nlayer+1, tmp_moist);
      if (cell->zwt_lumped_dist_zwt[zwtidx] == 999) cell->zwt_lumped_dist_zwt[zwtidx] = -total_depth*100; // in cm;

    }

  }

}
