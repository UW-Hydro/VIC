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

  Modified 4-14-98 to compute evapotranspiration within this routine
  in the hopes of reducing the number of iteration needed to find a
  solution surface temperature.                                 KAC

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
  double Vapor;		/** Total vapor flux from snow in m/s **/
  double moist;		/** layer moisture content in m/m **/
  double ice0;		/** layer ice content in m/m **/
  double max_moist;	/** layer maximum moisture content in m/m **/
  double bubble;	/** bubbling pressure in cm **/
  double             expt;
  double             surf_atten;
  double             wind;
  double             displacement;
  double             roughness;
  double             ref_height;
  double             elevation;
  double             b_infilt;
  double             max_infil;
  double             dt;
  double             vpd;
  double             rainfall;
  double             Wdew;
  double             snow_energy;
  double            *grnd_flux;
  double            *T1;
  double            *latent_heat;
  double            *sensible_heat;
  double            *deltaH;
  double            *store_error;
  double            *TMean;
  double            *rad;
  double            *depth;
  double            *Wcr;
  double            *Wpwp;
  layer_data_struct *layer;
  veg_var_struct    *veg_var;
  char               VEG;
  char               CALC_EVAP;
  int                veg_class;
  int                month;

  double error;
  double ice;
  double C1, C2, C3;
  double heat_capacity;
  double Evap;		/** Total evap in m/s **/

  T2            = (double) va_arg(ap, double);
  Ts_old        = (double) va_arg(ap, double);
  T1_old        = (double) va_arg(ap, double);
  Tair          = (double) va_arg(ap, double);
  ra            = (double) va_arg(ap, double);
  atmos_density = (double) va_arg(ap, double);
  shortwave     = (double) va_arg(ap, double);
  longwave      = (double) va_arg(ap, double);
  albedo        = (double) va_arg(ap, double);
  emissivity    = (double) va_arg(ap, double);
  kappa1        = (double) va_arg(ap, double);
  kappa2        = (double) va_arg(ap, double);
  Cs1           = (double) va_arg(ap, double);
  Cs2           = (double) va_arg(ap, double);
  D1            = (double) va_arg(ap, double);
  D2            = (double) va_arg(ap, double);
  dp            = (double) va_arg(ap, double);
  delta_t       = (double) va_arg(ap, double);
  Le            = (double) va_arg(ap, double);
  Ls            = (double) va_arg(ap, double);
  Vapor         = (double) va_arg(ap, double);
  moist         = (double) va_arg(ap, double);
  ice0          = (double) va_arg(ap, double);
  max_moist     = (double) va_arg(ap, double);
  bubble        = (double) va_arg(ap, double);
  expt          = (double) va_arg(ap, double);
  surf_atten    = (double) va_arg(ap, double);
  wind          = (double) va_arg(ap, double);
  displacement  = (double) va_arg(ap, double);
  roughness     = (double) va_arg(ap, double);
  ref_height    = (double) va_arg(ap, double);
  elevation     = (double) va_arg(ap, double);
  b_infilt      = (double) va_arg(ap, double);
  max_infil     = (double) va_arg(ap, double);
  dt            = (double) va_arg(ap, double);
  vpd           = (double) va_arg(ap, double);
  rainfall      = (double) va_arg(ap, double);
  Wdew          = (double) va_arg(ap, double);
  snow_energy   = (double) va_arg(ap, double);
  grnd_flux     = (double *) va_arg(ap, double *);
  T1            = (double *) va_arg(ap, double *);
  latent_heat   = (double *) va_arg(ap, double *);
  sensible_heat = (double *) va_arg(ap, double *);
  deltaH        = (double *) va_arg(ap, double *);
  store_error   = (double *) va_arg(ap, double *);
  TMean         = (double *) va_arg(ap, double *);
  rad           = (double *) va_arg(ap, double *);
  depth         = (double *) va_arg(ap, double *);
  Wcr           = (double *) va_arg(ap, double *);
  Wpwp          = (double *) va_arg(ap, double *);
  layer         = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var       = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG           = (char) va_arg(ap, char);
  CALC_EVAP     = (char) va_arg(ap, char);
  veg_class     = (int)  va_arg(ap, int);
  month         = (int)  va_arg(ap, int);

  /**********************************************
    Compute Surface Temperature at Half Time Step
    **********************************************/
  *TMean = 0.5 * (Ts + Ts_old);

  /*********************************************************************
    Estimate the Soil Temperature at the Interface of the Top Two Layers
    *********************************************************************/
  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);

  *T1 = (kappa1/2./D1/D2*(*TMean) + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  /*****************************************************
    Compute the Ground Heat Flux from the Top Soil Layer
    *****************************************************/
  *grnd_flux = surf_atten * (kappa1/D1*(*T1 - (*TMean)));

  /***************************
    Compute Evapotranspiration
    ***************************/
  *rad = (1.0 - albedo) * shortwave + longwave 
       - STEFAN_B * pow(*TMean+KELVIN,4.0) + *grnd_flux;
  if(VEG && CALC_EVAP)
    Evap = canopy_evap(layer,veg_var,TRUE,
		       veg_class,month,Wdew,Tair,dt,rad[0],vpd,
		       (1.0 - albedo) * shortwave,Tair,ra,
		       rainfall,displacement,roughness,ref_height,elevation,
		       depth,Wcr,Wpwp);
  else if(CALC_EVAP)
    Evap = arno_evap(layer, rad[0], Tair, vpd, (1.0 - albedo) * shortwave, 
		     D1, max_moist*depth[0]*1000., elevation, b_infilt, 
		     max_infil, 
		     Tair, displacement, roughness, ref_height, ra, dt);
  else Evap = 0.;
  
  /******************************************************
    Compute the Current Ice Content of the Top Soil Layer
    ******************************************************/
  if((*TMean+ *T1)/2.<0. && options.FROZEN_SOIL) {
    ice = moist - maximum_unfrozen_water((*TMean+ *T1)/2.,
          max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
 
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
    **********************************************************************/
  *latent_heat = -RHO_W*Le*Evap;
  *latent_heat += -atmos_density*Ls*Vapor;

  /************************************************
    Compute the Sensible Heat Flux from the Surface
    ************************************************/
  *sensible_heat = atmos_density*Cp*(Tair - (*TMean))/ra;

  /*******************************************
    Compute Heat Storage in the Top Soil Layer
    *******************************************/
  *deltaH = Cs1 * (Ts_old - *TMean) * D1 / delta_t;
  *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  /*************************************
    Compute Surface Energy Balance Error
    *************************************/
  error = (1.-albedo)*shortwave 
        + emissivity*(longwave-STEFAN_B*pow(*TMean + KELVIN,4.))
        + *sensible_heat + *latent_heat + *grnd_flux + *deltaH + snow_energy;

  *store_error = error;

  return error;

}
