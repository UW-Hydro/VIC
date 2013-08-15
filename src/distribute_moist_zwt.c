#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id: $";
int distribute_moist_zwt(cell_data_struct  *cell,
                         soil_con_struct   *soil_con)
/**********************************************************************
	distribute_moist_zwt	Ted Bohn		2013-Jul-25

  This function loops over zwt distribution and apportions soil moisture
  according to local water table depth.

**********************************************************************/
{

  extern option_struct options;
  int                ErrorFlag;
  double             tmp_total_max_moist;
  int                lindex;
  double             resid_moist[MAX_LAYERS];
  double             moist_distrib[MAX_NZWT];
  double             moist_distrib_mean;
  double             area_sum;
  int                zwtidx;
  double             tmp_moist;
  double             tmp_resid_moist;
  double             tmp_max_moist;
  double             mean_deficit;
  double             tmp_moist_dist_zwt[MAX_NZWT];
  double             slop;
  double             tmp_area;
  int                iter, max_iter = 10;
  double             slop_save;
  double             slop_per_area;
  double             deficit;
  double             column_slop;
  double             new_moist;
  int                frost_area;

  ErrorFlag = 0;

  tmp_total_max_moist = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    resid_moist[lindex] = soil_con->resid_moist[lindex]*soil_con->depth[lindex]*1000;
    tmp_total_max_moist += soil_con->max_moist[lindex];
  }

  // Compute total fractional area of exposed wetland (as fraction of tile area)
  moist_distrib[0] = 0.5*soil_con->ZwtDeltaMoist[0];
  moist_distrib_mean = moist_distrib[0]*soil_con->ZwtAreaFract[0];
  area_sum = soil_con->ZwtAreaFract[0];
  zwtidx = 1;
  while (zwtidx < options.Nzwt) {
    moist_distrib[zwtidx] = moist_distrib[zwtidx-1]+0.5*(soil_con->ZwtDeltaMoist[zwtidx-1]+soil_con->ZwtDeltaMoist[zwtidx]);
    moist_distrib_mean += moist_distrib[zwtidx]*soil_con->ZwtAreaFract[zwtidx];
    area_sum += soil_con->ZwtAreaFract[zwtidx];
    zwtidx++;
  }
  moist_distrib_mean /= area_sum;

  // Distribute total moisture

  // Distribute total moisture without regard to bounds
  tmp_moist = 0;
  tmp_resid_moist = 0;
  tmp_max_moist = 0;
  for (lindex=0; lindex<options.Nlayer; lindex++) {
    tmp_moist += cell->layer[lindex].moist;
    tmp_resid_moist += resid_moist[lindex];
    tmp_max_moist += soil_con->max_moist[lindex];
  }
  mean_deficit = tmp_max_moist - tmp_moist;

  // Shape of moisture curve is constant
  for (zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {
    tmp_moist_dist_zwt[zwtidx] = tmp_moist + moist_distrib[zwtidx]-moist_distrib_mean;
  }

  // Take care of moisture that went out of bounds
  slop = 0;
  tmp_area = 0;
  for (zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {
    if (tmp_moist_dist_zwt[zwtidx] > tmp_max_moist) {
      slop += (tmp_moist_dist_zwt[zwtidx] - tmp_max_moist)*soil_con->ZwtAreaFract[zwtidx];
      tmp_moist_dist_zwt[zwtidx] = tmp_max_moist;
    }
    else if (tmp_moist_dist_zwt[zwtidx] < tmp_resid_moist) {
      slop += (tmp_moist_dist_zwt[zwtidx] - tmp_resid_moist)*soil_con->ZwtAreaFract[zwtidx];
      tmp_moist_dist_zwt[zwtidx] = tmp_resid_moist;
    }
    else {
      tmp_area += soil_con->ZwtAreaFract[zwtidx];
    }
  }
  iter = 0;
  slop_save = slop;
  while (iter < max_iter && fabs(slop) > 1e-20 && tmp_area > 0) {
    slop_save = slop;
    slop_per_area = slop/tmp_area;
    for (zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {
      if (fabs(slop) > 1e-20 && slop*slop_save > 0) {
        if (slop > 0) {
          if (slop_per_area < tmp_max_moist-tmp_moist_dist_zwt[zwtidx]) {
            tmp_moist_dist_zwt[zwtidx] += slop_per_area;
            slop -= slop_per_area*soil_con->ZwtAreaFract[zwtidx];
          }
          else if (tmp_max_moist > tmp_moist_dist_zwt[zwtidx]) {
            slop -= (tmp_max_moist-tmp_moist_dist_zwt[zwtidx])*soil_con->ZwtAreaFract[zwtidx];
            tmp_moist_dist_zwt[zwtidx] = tmp_max_moist;
            tmp_area -= soil_con->ZwtAreaFract[zwtidx];
          }
        }
        else {
          if (tmp_moist_dist_zwt[zwtidx] + slop_per_area > tmp_resid_moist) {
            tmp_moist_dist_zwt[zwtidx] += slop_per_area;
            slop -= slop_per_area*soil_con->ZwtAreaFract[zwtidx];
          }
          else if (tmp_moist_dist_zwt[zwtidx] > tmp_resid_moist) {
            slop -= (tmp_resid_moist-tmp_moist_dist_zwt[zwtidx])*soil_con->ZwtAreaFract[zwtidx];
            tmp_moist_dist_zwt[zwtidx] = tmp_resid_moist;
            tmp_area -= soil_con->ZwtAreaFract[zwtidx];
          }
        }
      }
    }
    iter++;
  }

  // At this point, we know the horizontal distribution of deficit has the same mean as the original deficit -> deficit is conserved
  // Now distribute deficit and moisture vertically, conserving deficit horizontally
  // Attempt to distribute the deficit among layers in proportion to each layer's share of total deficit
  // But when a layer cannot accomodate a local deficit, we must pass its excess deficit on to the other layers
  // This may change the vertical distribution of moisture but should conserve moisture
  slop = 0;
  for (zwtidx=options.Nzwt-1; zwtidx>=0; zwtidx--) {
    deficit = tmp_max_moist - tmp_moist_dist_zwt[zwtidx];
    column_slop = 0;
    for (lindex=0; lindex<options.Nlayer; lindex++) {
      if (mean_deficit > 0)
        new_moist = soil_con->max_moist[lindex]-deficit*(soil_con->max_moist[lindex]-cell->layer[lindex].moist)/mean_deficit;
      else
        new_moist = soil_con->max_moist[lindex];
      new_moist += column_slop;
      column_slop = 0;
      if (new_moist > soil_con->max_moist[lindex]) {
        column_slop += (new_moist-soil_con->max_moist[lindex]);
        new_moist = soil_con->max_moist[lindex];
      }
      if (new_moist < resid_moist[lindex]) {
        column_slop -= (resid_moist[lindex]-new_moist);
        new_moist = resid_moist[lindex];
      }
      cell->layer[lindex].moist_dist_zwt[zwtidx] = new_moist;
    }
    slop += column_slop*soil_con->ZwtAreaFract[zwtidx];
  }
  if (fabs(slop)>1e-5) {
    fprintf(stderr,"ERROR: deficit distribution among layers messed up, error = %f\n",slop);
  }

  /** Distribute ice in constant proportion to moisture **/
  for (zwtidx=0; zwtidx<options.Nzwt; zwtidx++) {
    for (lindex=0; lindex<options.Nlayer; lindex++) {
      if (cell->layer[lindex].moist > 0) {
#if SPATIAL_FROST
        for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          cell->layer[lindex].ice_dist_zwt[zwtidx][frost_area] = cell->layer[lindex].ice[frost_area]*cell->layer[lindex].moist_dist_zwt[zwtidx]/cell->layer[lindex].moist;
#else
        cell->layer[lindex].ice_dist_zwt[zwtidx] = cell->layer[lindex].ice*cell->layer[lindex].moist_dist_zwt[zwtidx]/cell->layer[lindex].moist;
#endif
      }
      else {
#if SPATIAL_FROST
        for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          cell->layer[lindex].ice_dist_zwt[zwtidx][frost_area] = 0;
#else
        cell->layer[lindex].ice_dist_zwt[zwtidx] = 0;
#endif
      }
    }
  }

  return (ErrorFlag);

}


