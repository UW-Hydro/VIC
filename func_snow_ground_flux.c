#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double func_snow_ground_flux(double Ts, va_list ap) {
/**********************************************************************
	snow_ground_flux	Keith Cherkauer		January 24, 1997

  This subroutine computes the surface temeperature under a snow pack
  by balancing the ground heat flux, with the flux coming from the 
  snow pack.

  Reference:
	Gel'fan, A. N., "Comparison of Two Methods of Calculating
	Soil Freezing Depth," Meteorologiya i Gidrologiya, No. 2,
	pp 98-104, 1989.

  UNITS: mks

**********************************************************************/
  extern option_struct options;
  extern debug_struct debug;

  double T2;
  double Ts_old;
  double T1_old;
  double kappa1;        /** thermal conductivity of 1st layer */
  double kappa2;        /** thermal conductivity of 2nd layer */
  double Cs1;           /** volumetric heat capacity of 1st layer **/
  double Cs2;           /** volumetric heat capacity of 2nd layer **/
  double delta_t;       /** Time Step in Seconds **/
  double snow_density;	/** density of snow pack **/
  double snow_depth;	/** depth of snow pack **/
  double surf_temp;	/** surface temperature of snow pack **/
  double D1;		/** thickness of top layer **/
  double D2;		/** thickness of second layer **/
  double dp;		/** thermal damping depth of soil column **/
  double moist;		/** moisture content of top layer **/
  double ice0;		/** ice content of top layer **/
  double max_moist;
  double bubble;
  double expt;
  double *grnd_flux;
  double *deltaH;
  double *snow_flux;
  double *TMean;
  double *T1;

  double ice;
  double kappa_snow;	/* thermal conductivity of snow (W/s/K) */
  double C1, C2, C3;
  double error;

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

  *TMean = 0.5 * (Ts + Ts_old);

  kappa_snow = 2.9302e-6 * pow(snow_density, 2.0);
 
  *snow_flux = kappa_snow * (*TMean - surf_temp) / snow_depth;

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);
  *T1 = (kappa1/2./D1/D2*(*TMean) + C1/delta_t*T1_old 
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);
 
  if(options.FROZEN_SOIL && (*TMean+ *T1)/2.<0.) {
    ice = moist - maximum_unfrozen_water((*TMean+ *T1)/2.,max_moist,
					 bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;

  *deltaH = Cs1 * (Ts_old - *TMean) * D1 / delta_t;
  *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  *grnd_flux = kappa1/D1*(*T1 - *TMean);

  error = *deltaH + *grnd_flux - *snow_flux;

  return error;
 
}
