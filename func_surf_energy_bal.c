#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

double func_surf_energy_bal(double Ts, va_list ap)
/**********************************************************************
	func_surf_energy_bal	Keith Cherkauer		January 3, 1996

  This subroutine computes the surface energy balance for bare soil
  and vegetation uncovered by snow.  It computes outgoing longwave,
  sensible heat flux, ground heat flux, and storage of heat in the thin
  upper layer, based on the given surface temperature.

  The Energy Balance Equation used comes from Xu Liang's Paper 
  "Insights of the Ground Heat Flux in Land Surface Parameterization
  Schemes."

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  /** Thermal Properties **/
  double T2;		/** average soil temperature (C) **/
  double Ts_old;	/** last temperature (C) **/
  double T1_old;	/** last layer 1 soil temperature (C) **/
  double Tair;		/** Air Temperature **/
  double ra;		/** aerodynamic reisistance (s/m) **/
  double atmos_density;	/** atmospheric density (kg/m^3) **/
  double shortwave;
  double longwave;
  double albedo;
  double emissivity;
  double kappa1;	/** thermal conductivity of 1st layer */
  double kappa2;	/** thermal conductivity of 2nd layer */
  double Cs1; 		/** volumetric heat capacity of 1st layer **/
  double Cs2; 		/** volumetric heat capacity of 2nd layer **/
  double D1;		/** thickness of 1st layer (m) **/
  double D2;		/** thickness of 2nd layer (m) **/
  double dp;		/** depth to constant temperature (m) */
  double delta_t;	/** Time Step in Seconds **/
  double Le;		/** Latent heat of vapoization **/
  double Ls;		/** Latent heat of sublimation **/
  double Evap;		/** Total evap in m/s **/
  double Vapor;		/** Total vapor flux from snow in m/s **/
  double moist;		/** layer moisture content in m/m **/
  double ice0;		/** layer ice content in m/m **/
  double max_moist;	/** layer maximum moisture content in m/m **/
  double bubble;	/** bubbling pressure in cm **/
  double expt;
  double surf_atten;
  double wind;
  double displacement;
  double roughness;
  double ref_height;
  double dH_height;
  double melt_energy;
  double *grnd_flux;
  double *T1;
  double *latent_heat;
  double *sensible_heat;
  double *deltaH;
  double *store_error;
  char   veg_present;

  double Tmean;
  double error;
  double ice;
  double C1, C2, C3;
  double heat_capacity;

  T2 = (double) va_arg(ap, double);
  Ts_old = (double) va_arg(ap, double);
  T1_old = (double) va_arg(ap, double);
  Tair = (double) va_arg(ap, double);
  ra = (double) va_arg(ap, double);
  atmos_density = (double) va_arg(ap, double);
  shortwave = (double) va_arg(ap, double);
  longwave = (double) va_arg(ap, double);
  albedo = (double) va_arg(ap, double);
  emissivity = (double) va_arg(ap, double);
  kappa1 = (double) va_arg(ap, double);
  kappa2 = (double) va_arg(ap, double);
  Cs1 = (double) va_arg(ap, double);
  Cs2 = (double) va_arg(ap, double);
  D1 = (double) va_arg(ap, double);
  D2 = (double) va_arg(ap, double);
  dp = (double) va_arg(ap, double);
  delta_t = (double) va_arg(ap, double);
  Le = (double) va_arg(ap, double);
  Ls = (double) va_arg(ap, double);
  Evap = (double) va_arg(ap, double);
  Vapor = (double) va_arg(ap, double);
  moist = (double) va_arg(ap, double);
  ice0 = (double) va_arg(ap, double);
  max_moist = (double) va_arg(ap, double);
  bubble = (double) va_arg(ap, double);
  expt = (double) va_arg(ap, double);
  surf_atten = (double) va_arg(ap, double);
  wind = (double) va_arg(ap, double);
  displacement = (double) va_arg(ap, double);
  roughness = (double) va_arg(ap, double);
  ref_height = (double) va_arg(ap, double);
  dH_height = (double) va_arg(ap, double);
  melt_energy = (double) va_arg(ap, double);
  grnd_flux = (double *) va_arg(ap, double *);
  T1 = (double *) va_arg(ap, double *);
  latent_heat = (double *) va_arg(ap, double *);
  sensible_heat = (double *) va_arg(ap, double *);
  deltaH = (double *) va_arg(ap, double *);
  store_error = (double *) va_arg(ap, double *);
  veg_present = (char) va_arg(ap, char);

  Tmean = 0.5 * (Ts + Ts_old);

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);

  *T1 = (kappa1/2./D1/D2*Tmean + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  if((Tmean+ *T1)/2.<0. && options.FROZEN_SOIL) {
    ice = moist - maximum_unfrozen_water((Tmean+ *T1)/2.,
          max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
 
  *grnd_flux = surf_atten * ((kappa1/D1*(*T1-Tmean)));

  *latent_heat = -RHO_W*Le*Evap;
  *latent_heat += -atmos_density*Ls*Vapor;

  *sensible_heat = atmos_density*Cp*(Tair-Tmean)/ra;

  *deltaH = Cs1 * (Ts_old - Tmean) * D1 / delta_t;
  *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  error = (1.-albedo)*shortwave 
        + emissivity*(longwave-STEFAN_B*pow(Tmean + KELVIN,4.))
        + *sensible_heat + *latent_heat + *grnd_flux + *deltaH + melt_energy;

  *store_error = error;

  return error;

}
