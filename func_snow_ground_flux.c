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
  int                Nnodes;

  double ice;
  double kappa_snow;	/* thermal conductivity of snow (W/s/K) */
  double C1, C2, C3;
  double error;

  /** Initialize Variables **/
  T2             = (double) va_arg(ap, double);
  Ts_old         = (double) va_arg(ap, double);
  T1_old         = (double) va_arg(ap, double);
  kappa1         = (double) va_arg(ap, double);
  kappa2         = (double) va_arg(ap, double);
  Cs1            = (double) va_arg(ap, double);
  Cs2            = (double) va_arg(ap, double);
  delta_t        = (double) va_arg(ap, double);
  snow_density   = (double) va_arg(ap, double);
  snow_depth     = (double) va_arg(ap, double);
  surf_temp      = (double) va_arg(ap, double);
  D1             = (double) va_arg(ap, double);
  D2             = (double) va_arg(ap, double);
  dp             = (double) va_arg(ap, double);
  moist          = (double) va_arg(ap, double);
  ice0           = (double) va_arg(ap, double);
  max_moist      = (double) va_arg(ap, double);
  bubble         = (double) va_arg(ap, double);
  expt           = (double) va_arg(ap, double);
  grnd_flux      = (double *) va_arg(ap, double *);
  deltaH         = (double *) va_arg(ap, double *);
  snow_flux      = (double *) va_arg(ap, double *);
  TMean          = (double *) va_arg(ap, double *);
  T1             = (double *) va_arg(ap, double *);
  T_node         = (double *) va_arg(ap, double *);
  Tnew_node      = (double *) va_arg(ap, double *);
  dz_node        = (double *) va_arg(ap, double *);
  kappa_node     = (double *) va_arg(ap, double *);
  Cs_node        = (double *) va_arg(ap, double *);
  moist_node     = (double *) va_arg(ap, double *);
  expt_node      = (double *) va_arg(ap, double *);
  max_moist_node = (double *) va_arg(ap, double *);
  ice_node       = (double *) va_arg(ap, double *);
  alpha          = (double *) va_arg(ap, double *);
  beta           = (double *) va_arg(ap, double *);
  gamma          = (double *) va_arg(ap, double *);
  Nnodes         = (int)      va_arg(ap, int);

  *TMean = Ts;

  kappa_snow = 2.9302e-6 * pow(snow_density, 2.0);
 
  *snow_flux = kappa_snow * (*TMean - surf_temp) / snow_depth;

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
		    bubble,expt_node,ice_node,alpha,beta,gamma,Nnodes);
    *T1 = Tnew_node[1];
  }

  /************************************************************
    Compute the Ground Heat Flux Through the Upper Soil Layer
  ************************************************************/
  *grnd_flux = kappa1/D1*((*T1) - (*TMean));

  /**************************************************************
    Compute the Change in Heat Storage in the Upper Soil Layer
  **************************************************************/
  if(options.FROZEN_SOIL && (*TMean+ *T1)/2.<0.) {
    ice = moist - maximum_unfrozen_water((*TMean+ *T1)/2.,max_moist,
					 bubble,expt);
    if(ice<0.) ice=0.;
    if(ice>max_moist) ice=max_moist;
  }
  else ice=0.;

  *deltaH = Cs1 * ((Ts_old + T1_old)/2. - (*TMean + *T1)/2.) * D1 / delta_t;
  *deltaH -= ice_density*Lf*(ice0-ice)*D1/delta_t;

  /*******************************
    Compute Energy Balance Error
  *******************************/
  error = *deltaH + *grnd_flux - *snow_flux;

  return error;
 
}
