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
 
  if(surf_temp<=-9999) 
    error_calc_snow_ground_flux(surf_temp, T2, Ts_old, T1_old, 
				kappa1, kappa2, Cs1, Cs2,
				delta_t, snow_density, snow_depth, 
				Tsnow_surf, D1, D2, dp, moist, ice0, 
				max_moist, bubble, expt,
				&energy->grnd_flux, &energy->deltaH, 
				&energy->snow_flux, &energy->Trad[1],
				T1);

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

double error_calc_snow_ground_flux(double Tsurf, ...) {

  va_list ap;

  double  error;

  va_start(ap,Tsurf);
  error = error_print_snow_ground_flux(Tsurf, ap);
  va_end(ap);

  return error;

}

double error_print_snow_ground_flux(double Ts, va_list ap) {  

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
  double dp;
  double moist;	
  double ice0;
  double max_moist;
  double bubble;
  double expt;
  double *grnd_flux;
  double *deltaH;
  double *snow_flux;
  double *TMean;
  double *T1;

  /** Initialize Variables **/
  T2           = (double) va_arg(ap, double);
  Ts_old       = (double) va_arg(ap, double);
  T1_old       = (double) va_arg(ap, double);
  kappa1       = (double) va_arg(ap, double);
  kappa2       = (double) va_arg(ap, double);
  Cs1          = (double) va_arg(ap, double);
  Cs2          = (double) va_arg(ap, double);
  delta_t      = (double) va_arg(ap, double);
  snow_density = (double) va_arg(ap, double);
  snow_depth   =(double) va_arg(ap, double);
  surf_temp    = (double) va_arg(ap, double);
  D1           = (double) va_arg(ap, double);
  D2           = (double) va_arg(ap, double);
  dp           = (double) va_arg(ap, double);
  moist        = (double) va_arg(ap, double);
  ice0         = (double) va_arg(ap, double);
  max_moist    = (double) va_arg(ap, double);
  bubble       = (double) va_arg(ap, double);
  expt         = (double) va_arg(ap, double);
  grnd_flux    = (double *) va_arg(ap, double *);
  deltaH       = (double *) va_arg(ap, double *);
  snow_flux    = (double *) va_arg(ap, double *);
  TMean        = (double *) va_arg(ap, double *);
  T1           = (double *) va_arg(ap, double *);

  /* Print Variables */
  fprintf(stderr,"T2 = %lf\n",T2);
  fprintf(stderr,"Ts_old = %lf\n",Ts_old);
  fprintf(stderr,"T1_old = %lf\n",T1_old);
  fprintf(stderr,"kappa1 = %lf\n",kappa1);
  fprintf(stderr,"kappa2 = %lf\n",kappa2);
  fprintf(stderr,"Cs1 = %lf\n",Cs1);
  fprintf(stderr,"Cs2 = %lf\n",Cs2);
  fprintf(stderr,"delta_t = %lf\n",delta_t);
  fprintf(stderr,"snow_density = %lf\n",snow_density);
  fprintf(stderr,"snow_depth = %lf\n",snow_depth);
  fprintf(stderr,"surf_temp = %lf\n",surf_temp);
  fprintf(stderr,"D1 = %lf\n",D1);
  fprintf(stderr,"D2 = %lf\n",D2);
  fprintf(stderr,"dp = %lf\n",dp);
  fprintf(stderr,"moist = %lf\n",moist);
  fprintf(stderr,"ice0 = %lf\n",ice0);
  fprintf(stderr,"max_moist = %lf\n",max_moist);
  fprintf(stderr,"bubble = %lf\n",bubble);
  fprintf(stderr,"expt = %lf\n",expt);
  fprintf(stderr,"grnd_flux = %lf\n",grnd_flux[0]);
  fprintf(stderr,"deltaH = %lf\n",deltaH[0]);
  fprintf(stderr,"snow_flux = %lf\n",snow_flux[0]);
  fprintf(stderr,"TMean = %lf\n",TMean[0]);
  fprintf(stderr,"T1 = %lf\n",T1[0]);

  vicerror("Finished Dumping calc_snow_ground_flux variables");

}
