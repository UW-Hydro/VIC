#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
double calc_snow_ground_flux(double dp,
                             double moist,
                             double ice0,
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
  surf_temp = snow->surf_temp;
  D1 = soil_con.depth[0];
  D2 = soil_con.depth[1];
  max_moist = soil_con.max_moist[0]/(soil_con.depth[0]*1000.);
  bubble = soil_con.bubble;
  expt = soil_con.expt[0];

  /**************************************************
    Find Surface Temperature using Root Brent Method
  **************************************************/
  surf_temp = root_brent(energy->T[0]-25.,
              0.,func_snow_ground_flux,
              T2,Ts_old,T1_old,kappa1,kappa2,Cs1,Cs2,
              delta_t,snow_density,snow_depth,surf_temp,
              D1,D2,dp,moist,ice0,max_moist,bubble,expt,
              &energy->grnd_flux,&energy->deltaH,T1);
 
  /**************************************************
    Recalculate Energy Fluxes Based on Final Temperature
  **************************************************/
  surf_temp = 0.5 * (surf_temp + Ts_old);

  kappa_snow = 2.9302 * pow(snow_density*0.001, 2.0);

  snow_flux = kappa_snow * (surf_temp - surf_temp) / snow_depth;

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);
  *T1 = (kappa1/2./D1/D2*surf_temp + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  if((surf_temp+ *T1)/2.<0.) {
    ice = moist - maximum_unfrozen_water((surf_temp+ *T1)/2.,max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;

  energy->deltaH = Cs1 * (Ts_old - surf_temp) * D1 / delta_t;
  energy->deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  energy->grnd_flux = kappa1/D1*(*T1 - surf_temp);

  error = energy->deltaH + energy->grnd_flux - snow_flux;

  return (surf_temp);

}
