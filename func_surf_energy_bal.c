#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

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

  Modifications:
  04-14-98 modified to compute evapotranspiration within this routine
           in the hopes of reducing the number of iteration 
	  needed to find a solution surface temperature.       KAC
  07-13-98 modified to include elevation bands for vegetation 
           and snow                                             KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct  debug;

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
  double             snow_energy;
  double             mu;
  double            *rainfall;
  double            *Wdew;
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
  double            *T_node;
  double            *Tnew_node;
  double            *dz_node;
  double            *kappa_node;
  double            *Cs_node;
  double            *moist_node;
  double            *expt_node;
  double            *max_moist_node;
  double            *ice_node;
  double            *alpha;
  double            *beta;
  double            *gamma;
  float             *root;
  layer_data_struct *layer_wet;
  layer_data_struct *layer_dry;
  veg_var_struct    *veg_var_wet;
  veg_var_struct    *veg_var_dry;
  int                VEG;
  int                CALC_EVAP;
  int                veg_class;
  int                month;
  int                Nnodes;
  char              *FIRST_SOLN;

  double             error;
  double             ice;
  double             Evap;		/** Total evap in m/s **/
  double             tmp_rainfall[2];

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
  snow_energy   = (double) va_arg(ap, double);
  mu            = (double) va_arg(ap, double);
  rainfall      = (double *) va_arg(ap, double *);
  Wdew          = (double *) va_arg(ap, double *);
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
  T_node        = (double *) va_arg(ap, double *);
  Tnew_node     = (double *) va_arg(ap, double *);
  dz_node       = (double *) va_arg(ap, double *);
  kappa_node    = (double *) va_arg(ap, double *);
  Cs_node       = (double *) va_arg(ap, double *);
  moist_node    = (double *) va_arg(ap, double *);
  expt_node     = (double *) va_arg(ap, double *);
  max_moist_node= (double *) va_arg(ap, double *);
  ice_node      = (double *) va_arg(ap, double *);
  alpha         = (double *) va_arg(ap, double *);
  beta          = (double *) va_arg(ap, double *);
  gamma         = (double *) va_arg(ap, double *);
  root          = (float  *) va_arg(ap, float  *);
  layer_wet     = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry     = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet   = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry   = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG           = (int) va_arg(ap, int);
  CALC_EVAP     = (int) va_arg(ap, int);
  veg_class     = (int) va_arg(ap, int);
  month         = (int) va_arg(ap, int);
  Nnodes        = (int) va_arg(ap, int);
  FIRST_SOLN    = (char *)va_arg(ap, char *);

  /**********************************************
    Compute Surface Temperature at Half Time Step
    **********************************************/
  *TMean = Ts;

  if(!options.FROZEN_SOIL) {
    /*************************************************
      Use Xu's Equations to Calculate Thermal Fluxes
    *************************************************/
    *T1 = estimate_T1(TMean[0], T1_old, T2, D1, D2, kappa1, kappa2, Cs1, 
		      Cs2, dp, delta_t);
    
  }
  else {
    /*************************************************************
      Explicitly Solve Thermal Fluxes For all Soil Thermal Nodes 
    *************************************************************/
    T_node[0] = *TMean;
    solve_T_profile(Tnew_node,T_node,dz_node,kappa_node,Cs_node,
		    moist_node,delta_t,max_moist_node,
		    bubble,expt_node,ice_node,alpha,
		    beta,gamma,Nnodes,FIRST_SOLN,FALSE);
    *T1 = Tnew_node[1];
  }

  /*****************************************************
    Compute the Ground Heat Flux from the Top Soil Layer
  *****************************************************/
  *grnd_flux = surf_atten * (kappa1/D1*((*T1) - (*TMean)));

  /******************************************************
    Compute the Current Ice Content of the Top Soil Layer
  ******************************************************/
  if(options.FROZEN_SOIL && (*TMean+ *T1)/2.<0.) {
    ice = moist - maximum_unfrozen_water((*TMean+ *T1)/2.,
					 max_moist,bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;
 
  *deltaH = Cs1 * ((Ts_old + T1_old)/2. - (*TMean + *T1)/2.) * D1 / delta_t;
  *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  /***************************
    Compute Evapotranspiration
  ***************************/
  *rad = (1.0 - albedo) * shortwave + longwave 
       - STEFAN_B * pow(*TMean+KELVIN,4.0) + *grnd_flux + *deltaH;
  if(VEG && CALC_EVAP) {
    tmp_rainfall[WET] = rainfall[WET];
    tmp_rainfall[DRY] = rainfall[DRY];
    Evap = canopy_evap(layer_wet,layer_dry,veg_var_wet,
		       veg_var_dry,TRUE,
		       veg_class,month,mu,Wdew,Tair,dt,rad[0],
		       vpd,(1.0 - albedo) * shortwave,Tair,ra,
		       displacement,roughness,ref_height,elevation,
		       tmp_rainfall,depth,Wcr,Wpwp,root);
  }
  else if(CALC_EVAP)
    Evap = arno_evap(layer_wet, layer_dry, rad[0], Tair, vpd, 
		     (1.0 - albedo) * shortwave, 
		     D1, max_moist*depth[0]*1000., elevation, b_infilt, 
		     Tair, displacement, roughness, ref_height, ra, dt, mu);
  else Evap = 0.;
  
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
    **********************************************************************/
  *latent_heat = -RHO_W*Le*Evap;
  *latent_heat += -atmos_density*Ls*Vapor;

  /************************************************
    Compute the Sensible Heat Flux from the Surface
    ************************************************/
  *sensible_heat = atmos_density*Cp*(Tair - (*TMean))/ra;

  /*************************************
    Compute Surface Energy Balance Error
    *************************************/
  error = (1.-albedo)*shortwave 
        + emissivity*(longwave-STEFAN_B*pow(*TMean + KELVIN,4.))
        + *sensible_heat + *latent_heat + *grnd_flux + *deltaH + snow_energy;

  *store_error = error;

  return error;

}

double estimate_T1(double Ts, 
		   double T1_old,
		   double T2,
		   double D1, 
		   double D2, 
		   double kappa1, 
		   double kappa2, 
		   double Cs1, 
		   double Cs2, 
		   double dp,
		   double delta_t) {
/**********************************************************************
  estimate_T1                Keith Cherkauer          July 15, 1998

  uses Xu Liangs 3-layer energy balance formulation to estimate the 
  temperature between the first and second layers.  Formerly calculated
  independently in each of the surface energy balance equation routines.

**********************************************************************/

  double C1;
  double C2;
  double C3;
  double T1;

  C1 = Cs2 * dp / D2 * ( 1. - exp(-D2/dp));
  C2 = - ( 1. - exp(D1/dp) ) * exp(-D2/dp);
  C3 = kappa1/D1 - kappa2/D1 + kappa2/D1*exp(-D1/dp);

  T1 = (kappa1/2./D1/D2*(Ts) + C1/delta_t*T1_old
     + (2.*C2-1.+exp(-D1/dp))*kappa2/2./D1/D2*T2)
     / (C1/delta_t + kappa2/D1/D2*C2 + C3/2./D2);

  return(T1);

}
