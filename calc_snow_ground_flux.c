#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
double calc_snow_ground_flux(double dp,
                             double moist,
                             double ice0,
			     double Tsnow_surf,
                             double *T1,
                             energy_bal_struct *energy,
                             snow_data_struct *snow,
                             soil_con_struct soil_con,
                             global_param_struct gp) {
/**********************************************************************
  calc_snow_ground_flux.c	Keith Cherkauer		March 31, 1997

  This subroutine computes the surface temperature under the snowpack
  using the thermal flux through the snowpack.

**********************************************************************/

  double T2;
  double Ts_old;
  double T1_old;
  double kappa1;
  double kappa2;
  double Cs1;
  double Cs2;
  double delta_t;
  double snow_density;
  double snow_depth;
  double surf_temp;
  double D1;
  double D2;
  double max_moist;
  double bubble;
  double expt;

  double C1, C2, C3;
  double ice;
  double kappa_snow;
  double snow_flux;
  double error;

  /**************************************************
    Initialize Snow Flux Variables
  **************************************************/

  T2 = energy->T[gp.Ulayer+gp.Llayer+1];
  Ts_old = energy->T[0];
  T1_old = energy->T[1];
  kappa1 = energy->kappa[0];
  kappa2 = energy->kappa[1];
  Cs1 = energy->Cs[0];
  Cs2 = energy->Cs[1];
  delta_t = (double)(gp.dt)*3600.;
  snow_density = snow->density;
  snow_depth = snow->depth;
  D1 = soil_con.depth[0];
  D2 = soil_con.depth[1];
  max_moist = soil_con.max_moist[0]/(soil_con.depth[0]*1000.);
  bubble = soil_con.bubble;
  expt = soil_con.expt[0];

  /**************************************************
    Find Surface Temperature using Root Brent Method
  **************************************************/
  surf_temp = root_brent(energy->T[0]-25., 0., func_snow_ground_flux,
			 T2, Ts_old, T1_old, kappa1, kappa2, Cs1, Cs2,
			 delta_t, snow_density, snow_depth, Tsnow_surf,
			 D1, D2, dp,moist, ice0, max_moist, bubble, expt,
			 &energy->grnd_flux, &energy->deltaH, 
			 &energy->snow_flux, &energy->Trad[1], T1);
 
  /**************************************************
    Recalculate Energy Fluxes Based on Final Temperature
  **************************************************/
  error = solve_snow_ground_flux(surf_temp, T2, Ts_old, T1_old, 
				 kappa1, kappa2, Cs1, Cs2,
				 delta_t, snow_density, snow_depth, 
				 Tsnow_surf, D1, D2, dp, moist, ice0, 
				 max_moist, bubble, expt,
				 &energy->grnd_flux, &energy->deltaH, 
				 &energy->snow_flux, &energy->Trad[1],
				 T1);
 
  energy->error += error;

  return (surf_temp);

}

double solve_snow_ground_flux(double Tsurf, ...) {

  va_list ap;

  double  error;

  va_start(ap,Tsurf);
  error = func_snow_ground_flux(Tsurf, ap);
  va_end(ap);

  return error;

}
