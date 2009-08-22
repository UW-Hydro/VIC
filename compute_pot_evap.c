#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void compute_pot_evap(int veg_class, 
		      dmy_struct *dmy, 
		      int rec, 
		      int dt, 
		      double shortwave,
		      double net_longwave,
		      double tair, 
		      double vpd,
		      double elevation,
		      double **aero_resist,
		      double *pot_evap)
/****************************************************************************
                                                                           
  compute_pot_evap: computes potential evaporation for several different
                    reference land cover types (which are defined in
                    vicNl_def.h and global.h).

  modifications:
  2009-Aug-21 Fixed bug in assignment of albedo for natural vegetation.	TJB

****************************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern char ref_veg_ref_crop[];

  int NVegLibTypes;
  int i;
  double albedo;
  double net_short;
  double net_rad;
  double rs;
  double rarc;
  double RGL;
  double lai;
  double gsm_inv;
  char ref_crop;
  double rc;
  double ra;

  /************************************************
  Estimate and store potential evap estimates using penman equation
  ************************************************/

  NVegLibTypes = veg_lib[0].NVegLibTypes;
  for (i=0; i<N_PET_TYPES; i++) {
    if (i < N_PET_TYPES_NON_NAT) {
      rs = veg_lib[NVegLibTypes+i].rmin;
      rarc = veg_lib[NVegLibTypes+i].rarc;
      RGL = veg_lib[NVegLibTypes+i].RGL;
      lai = veg_lib[NVegLibTypes+i].LAI[dmy[rec].month-1];
      albedo = veg_lib[NVegLibTypes+i].albedo[dmy[rec].month-1];
    }
    else {
      rs = veg_lib[veg_class].rmin;
      if (i == PET_VEGNOCR) rs = 0;
      rarc = veg_lib[veg_class].rarc;
      RGL = veg_lib[veg_class].RGL;
      lai = veg_lib[veg_class].LAI[dmy[rec].month-1];
      albedo = veg_lib[veg_class].albedo[dmy[rec].month-1];
    }
    gsm_inv = 1.0;
    ref_crop = ref_veg_ref_crop[i];
    rc = calc_rc(rs, net_short, RGL, tair, vpd, lai, gsm_inv, ref_crop);
    if (i < N_PET_TYPES_NON_NAT || !veg_lib[veg_class].overstory)
      ra = aero_resist[i][0];
    else
      ra = aero_resist[i][1];
    net_short = (1.0 - albedo) * shortwave;
    net_rad = net_short + net_longwave;
    pot_evap[i] = penman(tair, elevation, net_rad, vpd, ra, rc, rarc) * dt/24.0;
  }

}
