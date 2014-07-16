#include <vic_def.h>
#include <vic_run.h>

/* One element for each non-natural PET type */
char   ref_veg_over[]        = { 0, 0, 0, 0 };
double ref_veg_rarc[]        = { 0.0, 0.0, 25, 25 };
double ref_veg_rmin[]        = { 0.0, 0.0, 100, 100 };
double ref_veg_lai[]         = { 1.0, 1.0, 2.88, 4.45 };
double ref_veg_albedo[]      = { BARE_SOIL_ALBEDO, H2O_SURF_ALBEDO, 0.23, 0.23 };
double ref_veg_rough[]       = { 0.001, 0.001, 0.0148, 0.0615 };
double ref_veg_displ[]       = { 0.0054, 0.0054, 0.08, 0.3333 };
double ref_veg_wind_h[]      = { 10.0, 10.0, 10.0, 10.0 };
double ref_veg_RGL[]         = { 0.0, 0.0, 100, 100 };
double ref_veg_rad_atten[]   = { 0.0, 0.0, 0.0, 0.0 };
double ref_veg_wind_atten[]  = { 0.0, 0.0, 0.0, 0.0 };
double ref_veg_trunk_ratio[] = { 0.0, 0.0, 0.0, 0.0 };
/* One element for each PET type (non-natural or natural) */
char ref_veg_ref_crop[] = { FALSE, FALSE, TRUE, TRUE, FALSE, FALSE };

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
  extern veg_lib_struct *vic_run_veg_lib;

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

  NVegLibTypes = vic_run_veg_lib[0].NVegLibTypes;
  for (i=0; i<N_PET_TYPES; i++) {
    if (i < N_PET_TYPES_NON_NAT) {
      rs = vic_run_veg_lib[NVegLibTypes+i].rmin;
      rarc = vic_run_veg_lib[NVegLibTypes+i].rarc;
      RGL = vic_run_veg_lib[NVegLibTypes+i].RGL;
      lai = vic_run_veg_lib[NVegLibTypes+i].LAI[dmy[rec].month-1];
      albedo = vic_run_veg_lib[NVegLibTypes+i].albedo[dmy[rec].month-1];
    }
    else {
      rs = vic_run_veg_lib[veg_class].rmin;
      if (i == PET_VEGNOCR) rs = 0;
      rarc = vic_run_veg_lib[veg_class].rarc;
      RGL = vic_run_veg_lib[veg_class].RGL;
      lai = vic_run_veg_lib[veg_class].LAI[dmy[rec].month-1];
      albedo = vic_run_veg_lib[veg_class].albedo[dmy[rec].month-1];
    }
    gsm_inv = 1.0;
    ref_crop = ref_veg_ref_crop[i];
    rc = calc_rc(rs, net_short, RGL, tair, vpd, lai, gsm_inv, ref_crop);
    if (i < N_PET_TYPES_NON_NAT || !vic_run_veg_lib[veg_class].overstory)
      ra = aero_resist[i][0];
    else
      ra = aero_resist[i][1];
    net_short = (1.0 - albedo) * shortwave;
    net_rad = net_short + net_longwave;
    pot_evap[i] = penman(tair, elevation, net_rad, vpd, ra, rc, rarc) * dt/24.0;
  }

}