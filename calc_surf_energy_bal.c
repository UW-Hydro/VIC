#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

double calc_surf_energy_bal(int                rec,
			    int                iveg,
			    int                nlayer,
			    int                Nveg,
			    int                dt,
			    int                Nnodes,
			    int                veg_class,
			    int                band,
			    int                hour,
			    double             ice0,
			    double             moist,
			    double             dp,
			    double             surf_atten,
			    double             T0,
			    double             shortwave,
			    double             longwave,
			    double             air_temp,
			    double             Le,
			    double             Ls,
			    double             mu,
			    double             displacement,
			    double             roughness,
			    double             ref_height,
			    double             snow_energy,
			    double             vapor_flux,
			    double             bare_albedo,
			    double            *aero_resist,
			    double            *wind,
			    double            *rainfall,
			    double            *ppt,
			    float             *root,
			    atmos_data_struct *atmos,
			    veg_var_struct    *veg_var_wet,
			    veg_var_struct    *veg_var_dry,
			    energy_bal_struct *energy,
			    snow_data_struct  *snow,
			    layer_data_struct *layer_wet,
			    layer_data_struct *layer_dry,
			    soil_con_struct   *soil_con,
			    dmy_struct        *dmy)
/**************************************************************
  calc_surf_energy_bal.c  Greg O'Donnell and Keith Cherkauer  Sept 9 1997
  
  This function calculates the surface temperature, in the
  case of no snow cover.  Evaporation is computed using the
  previous ground heat flux, and then used to comput latent heat
  in the energy balance routine.  Surface temperature is found
  using the Root Brent method (Numerical Recipies).
  

  modifications:
    02-29-00  Included variables needed to compute energy flux
              through the snow pack.  The ground surface energy
              balance will now be a mixture of snow covered
	      and bare ground, controlled by the snow cover 
	      fraction set in solve_snow.c                 KAC
    04-Jun-04 Placed "ERROR" at beginning of screen dump in
	      error_print_surf_energy_bal.		TJB

***************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;

  char     FIRST_SOLN[1];
  double   T2;
  double   Ts_old;
  double   T1_old;
  double   Tair;
  double   ra;
  double   atmos_density;
  double   albedo;
  double   emissivity;
  double   kappa1;
  double   kappa2;
  double   Cs1;
  double   Cs2;
  double   D1;
  double   D2;
  double   delta_t;
  double   Vapor;
  double   max_moist;
  double   bubble;
  double   expt;
  double   T1;
  double   snow_depth;
  double   snow_density;
  double   Tsnow_surf;
  double   snow_cover_fraction;
  double   snow_flux;

  int      VEG;
  double   Tsurf;
  double   U;
  double   error;
  double   Wdew[2];
  double  *T_node;
  double  Tnew_node[MAX_NODES];
  double  *dz_node;
  double  *kappa_node;
  double  *Cs_node;
  double  *moist_node;
  double  *bubble_node;
  double  *expt_node;
  double  *max_moist_node;
  double  *ice_node;
  double  *alpha;
  double  *beta;
  double  *gamma;
  layer_data_struct layer[MAX_LAYERS];

  double   T_lower, T_upper;

  /**************************************************
    Correct Aerodynamic Resistance for Stability
  **************************************************/
  ra = aero_resist[0];
  U = wind[0];
  if (U > 0.0)
    ra /= StabilityCorrection(ref_height, displacement, T0,
          air_temp, U, roughness);
  else
    ra = HUGE_RESIST;
  
  /**************************************************
    Compute Evaporation and Transpiration 
  **************************************************/
  if(iveg!=Nveg) {
    if(veg_lib[veg_class].LAI[dmy->month-1] > 0.0) VEG = TRUE;
    else VEG = FALSE;
  }
  else VEG = FALSE;
  Vapor = vapor_flux / (double)dt / 3600.;

  /**************************************************
    Set All Variables For Use
  **************************************************/

  T2                  = energy->T[Nnodes-1];
  Ts_old              = energy->T[0];
  T1_old              = energy->T[1];
  Tair                = air_temp;
  atmos_density       = atmos->density[hour];
  albedo              = bare_albedo;
  emissivity          = 1.;
  kappa1              = energy->kappa[0];
  kappa2              = energy->kappa[1];
  Cs1                 = energy->Cs[0];
  Cs2                 = energy->Cs[1];
  D1                  = soil_con->depth[0];
/*   D2                  = soil_con->depth[1]; */
  D2                  = soil_con->depth[0];
  delta_t             = (double)dt*3600.;
  max_moist           = soil_con->max_moist[0]/(soil_con->depth[0]*1000.);
  bubble              = soil_con->bubble[0];
  expt                = soil_con->expt[0];
  snow_depth          = snow->depth;
  snow_density        = snow->density;
  Tsnow_surf          = snow->surf_temp;
  snow_cover_fraction = snow->coverage;
  Wdew[WET]           = veg_var_wet->Wdew;
  if(options.DIST_PRCP) Wdew[DRY] = veg_var_dry->Wdew;
  FIRST_SOLN[0] = TRUE;

  /*************************************************************
    Prepare soil node variables for finite difference solution
  *************************************************************/

  if(!options.QUICK_FLUX) {

    bubble_node    = soil_con->bubble_node; 
    expt_node      = soil_con->expt_node; 
    max_moist_node = soil_con->max_moist_node;  
    alpha          = soil_con->alpha; 
    beta           = soil_con->beta; 
    gamma          = soil_con->gamma; 
    moist_node     = energy->moist;
    kappa_node     = energy->kappa_node;
    Cs_node        = energy->Cs_node;
    T_node         = energy->T;
    dz_node        = soil_con->dz_node;
    ice_node       = energy->ice;

  }
  else {

    bubble_node    = NULL; 
    expt_node      = NULL; 
    max_moist_node = NULL;  
    alpha          = NULL; 
    beta           = NULL; 
    gamma          = NULL; 
    moist_node     = NULL;
    kappa_node     = NULL;
    Cs_node        = NULL;
    T_node         = NULL;
    dz_node        = NULL;
    ice_node       = NULL;

  }

  /**************************************************
    Find Surface Temperature Using Root Brent Method
  **************************************************/
  if(options.FULL_ENERGY) {

    /** Added for temporary backwards compatability **/
    if(snow->swq > 0) {
      T_lower = 0.;
      T_upper = energy->T[0]-SURF_DT;
    }
    else {
      T_lower = 0.5*(energy->T[0]+air_temp)-SURF_DT;
      T_upper = 0.5*(energy->T[0]+air_temp)+SURF_DT;
    }

#if QUICK_FS
    Tsurf = root_brent(T_upper, T_lower, 
		       func_surf_energy_bal, T2, Ts_old, T1_old, Tair, 
		       ra, atmos_density, shortwave, longwave, albedo, 
		       emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, dp, 
		       delta_t, Le, Ls, Vapor, moist, ice0, max_moist, 
		       bubble, expt, snow_depth, snow_density, Tsnow_surf, 
		       snow_cover_fraction, surf_atten, U, displacement, 
		       roughness, ref_height, (double)soil_con->elevation, 
		       soil_con->b_infilt, soil_con->max_infil, 
		       (double)dt, atmos->vpd[hour], snow_energy, mu, 
		       rainfall, Wdew, &energy->grnd_flux, &T1, 
		       &energy->latent, &energy->sensible, 
		       &energy->deltaH, &energy->snow_flux,
		       &energy->error, soil_con->depth, 
		       soil_con->Wcr, soil_con->Wpwp, 
		       soil_con->resid_moist, T_node, Tnew_node, dz_node, 
		       kappa_node, Cs_node, moist_node, bubble_node, 
		       expt_node, max_moist_node, ice_node, alpha, beta, 
		       gamma, soil_con->ufwc_table_layer[0], 
		       soil_con->ufwc_table_node, root, layer_wet, 
		       layer_dry, veg_var_wet, veg_var_dry, VEG, 
		       veg_class, dmy->month, Nnodes, 
		       FIRST_SOLN, snow->snow, soil_con->FS_ACTIVE);
#else
    Tsurf = root_brent(T_upper, T_lower, 
		       func_surf_energy_bal, T2, Ts_old, T1_old, Tair, 
		       ra, atmos_density, shortwave, longwave, albedo, 
		       emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, dp, 
		       delta_t, Le, Ls, Vapor, moist, ice0, max_moist, 
		       bubble, expt, snow_depth, snow_density, Tsnow_surf, 
		       snow_cover_fraction, surf_atten, U, displacement, 
		       roughness, ref_height, (double)soil_con->elevation, 
		       soil_con->b_infilt, soil_con->max_infil, 
		       (double)dt, atmos->vpd[hour], snow_energy, mu, 
		       rainfall, Wdew, &energy->grnd_flux, &T1, 
		       &energy->latent, &energy->sensible, &energy->deltaH, 
		       &energy->snow_flux, &energy->error, soil_con->depth, 
		       soil_con->Wcr, soil_con->Wpwp, 
		       soil_con->resid_moist, T_node, Tnew_node, dz_node, 
		       kappa_node, Cs_node, moist_node, bubble_node, 
		       expt_node, max_moist_node, ice_node, alpha, beta, 
		       gamma, root, layer_wet, layer_dry, veg_var_wet, 
		       veg_var_dry, VEG, veg_class, 
		       dmy->month, Nnodes, FIRST_SOLN, snow->snow, 
		       soil_con->FS_ACTIVE);
#endif

    if(Tsurf <= -9998) {  
#if QUICK_FS
      error = error_calc_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, 
					 ra, atmos_density, shortwave, 
					 longwave, albedo, emissivity, kappa1, 
					 kappa2, Cs1, Cs2, D1, D2, dp, 
					 delta_t, Le, Ls, Vapor, moist, ice0, 
					 max_moist, bubble, expt, snow_depth, 
					 snow_density, Tsnow_surf, 
					 snow_cover_fraction, surf_atten, 
					 U, displacement, roughness, 
					 ref_height, 
					 (double)soil_con->elevation, 
					 soil_con->b_infilt, 
					 soil_con->max_infil, (double)dt, 
					 atmos->vpd[hour], snow_energy, mu, 
					 rainfall, Wdew, &energy->grnd_flux, 
					 &T1, &energy->latent, 
					 &energy->sensible, &energy->deltaH, 
					 &energy->snow_flux,
					 &energy->error, soil_con->depth, 
					 soil_con->Wcr, soil_con->Wpwp, 
					 soil_con->resid_moist, T_node, 
					 Tnew_node, dz_node, kappa_node, 
					 Cs_node, moist_node, bubble_node, 
					 expt_node, max_moist_node, 
					 ice_node, alpha, beta, gamma, 
					 soil_con->ufwc_table_layer[0], 
					 soil_con->ufwc_table_node, 
					 root, layer_wet, layer_dry, 
					 veg_var_wet, veg_var_dry, VEG, 
					 veg_class, 
					 dmy->month, iveg, Nnodes, FIRST_SOLN,
					 snow->snow, soil_con->FS_ACTIVE);
#else
      error = error_calc_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, 
					 ra, atmos_density, shortwave, 
					 longwave, albedo, emissivity, kappa1, 
					 kappa2, Cs1, Cs2, D1, D2, dp, 
					 delta_t, Le, Ls, Vapor, moist, ice0, 
					 max_moist, bubble, expt, snow_depth, 
					 snow_density, Tsnow_surf, 
					 snow_cover_fraction, surf_atten, 
					 U, displacement, roughness, 
					 ref_height, 
					 (double)soil_con->elevation, 
					 soil_con->b_infilt, 
					 soil_con->max_infil, (double)dt, 
					 atmos->vpd[hour], snow_energy, mu, 
					 rainfall, Wdew, &energy->grnd_flux, 
					 &T1, &energy->latent, 
					 &energy->sensible, &energy->deltaH, 
					 &energy->snow_flux,
					 &energy->error, soil_con->depth, 
					 soil_con->Wcr, soil_con->Wpwp, 
					 soil_con->resid_moist, T_node, 
					 Tnew_node, dz_node, kappa_node, 
					 Cs_node, moist_node, bubble_node, 
					 expt_node, max_moist_node, ice_node, 
					 alpha, beta, gamma, root, layer_wet, 
					 layer_dry, veg_var_wet, veg_var_dry, 
					 VEG, veg_class, 
					 dmy->month, iveg, Nnodes, FIRST_SOLN, 
					 snow->snow, soil_con->FS_ACTIVE);
#endif
    }
  }
  else {

    /** Frozen soil model run with no surface energy balance **/
    Tsurf = Tair;

  }

  /**************************************************
    Recalculate Energy Balance Terms for Final Surface Temperature
  **************************************************/
#if QUICK_FS
  error = solve_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, ra, 
				atmos_density, shortwave, longwave, albedo, 
				emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, 
				dp, delta_t, Le, Ls, Vapor, moist, ice0, 
				max_moist, bubble, expt, snow_depth, 
				snow_density, Tsnow_surf, 
				snow_cover_fraction, surf_atten, U, 
				displacement, roughness, ref_height, 
				(double)soil_con->elevation, 
				soil_con->b_infilt, soil_con->max_infil, 
				(double)dt, atmos->vpd[hour], snow_energy, mu, 
				rainfall, Wdew, &energy->grnd_flux, 
				&T1, &energy->latent, &energy->sensible, 
				&energy->deltaH, &energy->snow_flux, 
				&energy->error, soil_con->depth, 
				soil_con->Wcr, soil_con->Wpwp, 
				soil_con->resid_moist, T_node, 
				Tnew_node, dz_node, kappa_node, Cs_node, 
				moist_node, bubble_node, expt_node, 
				max_moist_node, ice_node, alpha, beta, gamma,  
				soil_con->ufwc_table_layer[0],  
				soil_con->ufwc_table_node, root, layer_wet, 
				layer_dry, veg_var_wet, veg_var_dry, VEG, 
				veg_class, dmy->month, Nnodes, 
				FIRST_SOLN, snow->snow, soil_con->FS_ACTIVE);
#else
  error = solve_surf_energy_bal(Tsurf, T2, Ts_old, T1_old, Tair, ra, 
				atmos_density, shortwave, longwave, albedo, 
				emissivity, kappa1, kappa2, Cs1, Cs2, D1, D2, 
				dp, delta_t, Le, Ls, Vapor, moist, ice0, 
				max_moist, bubble, expt, snow_depth, 
				snow_density, Tsnow_surf, 
				snow_cover_fraction, surf_atten, U, 
				displacement, roughness, ref_height, 
				(double)soil_con->elevation, 
				soil_con->b_infilt, soil_con->max_infil, 
				(double)dt, atmos->vpd[hour], snow_energy, mu, 
				rainfall, Wdew, &energy->grnd_flux, 
				&T1, &energy->latent, &energy->sensible, 
				&energy->deltaH, &energy->snow_flux, 
				&energy->error, soil_con->depth, 
				soil_con->Wcr, soil_con->Wpwp, 
				soil_con->resid_moist, T_node, 
				Tnew_node, dz_node, kappa_node, Cs_node, 
				moist_node, bubble_node, expt_node, 
				max_moist_node, ice_node, alpha, beta, gamma, 
				root, layer_wet, layer_dry, veg_var_wet, 
				veg_var_dry, VEG, 
				veg_class, dmy->month, Nnodes, FIRST_SOLN, 
				snow->snow, soil_con->FS_ACTIVE);
#endif
  
  energy->error = error;

  /***************************************************
    Recalculate Soil Moisture and Thermal Properties
  ***************************************************/
  if(options.GRND_FLUX) {
    if(options.QUICK_FLUX) {
      
      energy->T[0] = Tsurf;
      energy->T[1] = T1;
      
    }
    else {
      
      finish_frozen_soil_calcs(energy, layer_wet, layer_dry, layer, soil_con, 
			       Nnodes, iveg, mu, Tnew_node, kappa_node, 
			       Cs_node, moist_node);
      
    }
    
  }
  else {

    energy->T[0] = Tsurf;

  }

  /** Store precipitation that reaches the surface */
  if(!snow->snow) {
    if(iveg!=Nveg) {
      if(veg_lib[veg_class].LAI[dmy->month-1] <= 0.0) { 
	veg_var_wet->throughfall = rainfall[WET];
	if(options.DIST_PRCP) veg_var_dry->throughfall = rainfall[DRY];
	ppt[WET] = veg_var_wet->throughfall;
	if(options.DIST_PRCP) ppt[DRY] = veg_var_dry->throughfall;
      }
      else {
	ppt[WET] = veg_var_wet->throughfall;
	if(options.DIST_PRCP) ppt[DRY] = veg_var_dry->throughfall;
      }
    }
    else {
      ppt[WET] = rainfall[WET];
      if(options.DIST_PRCP) ppt[DRY] = rainfall[DRY];
    }
  }
   
  /** Store net longwave radiation **/
  if(hour < 24)
    energy->longwave = longwave 
      - STEFAN_B * (Tsurf+KELVIN) * (Tsurf+KELVIN) 
		    * (Tsurf+KELVIN) * (Tsurf+KELVIN);
  else energy->longwave = longwave;

  /** Store surface albedo **/
  energy->albedo = bare_albedo;

  /** Store net shortwave radiation **/
  energy->shortwave = (1. - energy->albedo) * shortwave;

  /** Return soil surface temperature **/
  return (Tsurf);
    
}

double solve_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap, Tsurf);
  error = func_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double error_calc_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap, Tsurf);
  error = error_print_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double error_print_surf_energy_bal(double Ts, va_list ap) {

  extern option_struct options;

  /** Thermal Properties **/
  double             T2; 
  double             Ts_old; 
  double             T1_old; 
  double             Tair; 
  double             ra; 
  double             atmos_density; 
  double             shortwave;
  double             longwave;
  double             albedo;
  double             emissivity;
  double             kappa1; 
  double             kappa2; 
  double             Cs1; 
  double             Cs2; 
  double             D1; 
  double             D2; 
  double             dp; 
  double             delta_t;
  double             Le; 
  double             Ls; 
  double             Vapor; 
  double             moist; 
  double             ice0; 
  double             max_moist; 
  double             bubble; 
  double             expt;
  double             snow_depth;
  double             snow_density;
  double             Tsnow_surf;
  double             snow_cover_fraction;
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
  double            *snow_flux;
  double            *store_error;
  double            *depth;
  double            *Wcr;
  double            *Wpwp;
  double            *resid_moist;
  double            *T_node;
  double            *Tnew_node;
  double            *dz_node;
  double            *kappa_node;
  double            *Cs_node;
  double            *moist_node;
  double            *bubble_node;
  double            *expt_node;
  double            *max_moist_node;
  double            *ice_node;
  double            *alpha;
  double            *beta;
  double            *gamma;
#if QUICK_FS
  double           **ufwc_table_layer;
  double          ***ufwc_table_node;
#endif
  float             *root;
  layer_data_struct *layer_wet;
  layer_data_struct *layer_dry;
  veg_var_struct    *veg_var_wet;
  veg_var_struct    *veg_var_dry;
  int                VEG;
  int                veg_class;
  int                month;
  int                iveg;
  int                Nnodes;
  char              *FIRST_SOLN;
  int                SNOWING;
  int                FS_ACTIVE;

  int                i;

  /* Initialize Variables */
  T2                  = (double) va_arg(ap, double);
  Ts_old              = (double) va_arg(ap, double);
  T1_old              = (double) va_arg(ap, double);
  Tair                = (double) va_arg(ap, double);
  ra                  = (double) va_arg(ap, double);
  atmos_density       = (double) va_arg(ap, double);
  shortwave           = (double) va_arg(ap, double);
  longwave            = (double) va_arg(ap, double);
  albedo              = (double) va_arg(ap, double);
  emissivity          = (double) va_arg(ap, double);
  kappa1              = (double) va_arg(ap, double);
  kappa2              = (double) va_arg(ap, double);
  Cs1                 = (double) va_arg(ap, double);
  Cs2                 = (double) va_arg(ap, double);
  D1                  = (double) va_arg(ap, double);
  D2                  = (double) va_arg(ap, double);
  dp                  = (double) va_arg(ap, double);
  delta_t             = (double) va_arg(ap, double);
  Le                  = (double) va_arg(ap, double);
  Ls                  = (double) va_arg(ap, double);
  Vapor               = (double) va_arg(ap, double);
  moist               = (double) va_arg(ap, double);
  ice0                = (double) va_arg(ap, double);
  max_moist           = (double) va_arg(ap, double);
  bubble              = (double) va_arg(ap, double);
  expt                = (double) va_arg(ap, double);
  snow_depth          = (double) va_arg(ap, double);
  snow_density        = (double) va_arg(ap, double);
  Tsnow_surf          = (double) va_arg(ap, double);
  snow_cover_fraction = (double) va_arg(ap, double);
  surf_atten          = (double) va_arg(ap, double);
  wind                = (double) va_arg(ap, double);
  displacement        = (double) va_arg(ap, double);
  roughness           = (double) va_arg(ap, double);
  ref_height          = (double) va_arg(ap, double);
  elevation           = (double) va_arg(ap, double);
  b_infilt            = (double) va_arg(ap, double);
  max_infil           = (double) va_arg(ap, double);
  dt                  = (double) va_arg(ap, double);
  vpd                 = (double) va_arg(ap, double);
  snow_energy         = (double) va_arg(ap, double);
  mu                  = (double) va_arg(ap, double);
  rainfall            = (double *) va_arg(ap, double *);
  Wdew                = (double *) va_arg(ap, double *);
  grnd_flux           = (double *) va_arg(ap, double *);
  T1                  = (double *) va_arg(ap, double *);
  latent_heat         = (double *) va_arg(ap, double *);
  sensible_heat       = (double *) va_arg(ap, double *);
  deltaH              = (double *) va_arg(ap, double *);
  snow_flux           = (double *) va_arg(ap, double *);
  store_error         = (double *) va_arg(ap, double *);
  depth               = (double *) va_arg(ap, double *);
  Wcr                 = (double *) va_arg(ap, double *);
  Wpwp                = (double *) va_arg(ap, double *);
  resid_moist         = (double *) va_arg(ap, double *);
  T_node              = (double *) va_arg(ap, double *);
  Tnew_node           = (double *) va_arg(ap, double *);
  dz_node             = (double *) va_arg(ap, double *);
  kappa_node          = (double *) va_arg(ap, double *);
  Cs_node             = (double *) va_arg(ap, double *);
  moist_node          = (double *) va_arg(ap, double *);
  bubble_node         = (double *) va_arg(ap, double *);
  expt_node           = (double *) va_arg(ap, double *);
  max_moist_node      = (double *) va_arg(ap, double *);
  ice_node            = (double *) va_arg(ap, double *);
  alpha               = (double *) va_arg(ap, double *);
  beta                = (double *) va_arg(ap, double *);
  gamma               = (double *) va_arg(ap, double *);
#if QUICK_FS
  ufwc_table_layer    = (double **) va_arg(ap, double **);
  ufwc_table_node     = (double ***) va_arg(ap, double ***);
#endif
  root                = (float  *) va_arg(ap, float  *);
  layer_wet           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry           = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry         = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG                 = (int) va_arg(ap, int);
  veg_class           = (int) va_arg(ap, int);
  month               = (int) va_arg(ap, int);
  iveg                = (int) va_arg(ap, int);
  Nnodes              = (int) va_arg(ap, int);
  FIRST_SOLN          = (char *) va_arg(ap, char *);
  SNOWING             = (int) va_arg(ap, int);
  FS_ACTIVE           = (int) va_arg(ap, int);

  /* Print Variables */
  fprintf(stderr, "ERROR: calc_surf_energy_bal failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  fprintf(stderr,"T2 = %f\n",T2);
  fprintf(stderr,"Ts_old = %f\n",Ts_old);
  fprintf(stderr,"T1_old = %f\n",T1_old);
  fprintf(stderr,"Tair = %f\n",Tair);
  fprintf(stderr,"ra = %f\n",ra);
  fprintf(stderr,"atmos_density = %f\n",atmos_density);
  fprintf(stderr,"shortwave = %f\n",shortwave);
  fprintf(stderr,"longwave = %f\n",longwave);
  fprintf(stderr,"albedo = %f\n",albedo);
  fprintf(stderr,"emissivity = %f\n",emissivity);
  fprintf(stderr,"kappa1 = %f\n",kappa1);
  fprintf(stderr,"kappa2 = %f\n",kappa2);
  fprintf(stderr,"Cs1 = %f\n",Cs1);
  fprintf(stderr,"Cs2 = %f\n",Cs2);
  fprintf(stderr,"D1 = %f\n",D1);
  fprintf(stderr,"D2 = %f\n",D2);
  fprintf(stderr,"dp = %f\n",dp);
  fprintf(stderr,"delta_t = %f\n",delta_t);
  fprintf(stderr,"Le = %f\n",Le);
  fprintf(stderr,"Ls = %f\n",Ls);
  fprintf(stderr,"Vapor = %f\n",Vapor);
  fprintf(stderr,"moist = %f\n",moist);
  fprintf(stderr,"ice0 = %f\n",ice0);
  fprintf(stderr,"max_moist = %f\n",max_moist);
  fprintf(stderr,"bubble = %f\n",bubble);
  fprintf(stderr,"expt = %f\n",expt);
  fprintf(stderr,"surf_atten = %f\n",surf_atten);
  fprintf(stderr,"wind = %f\n",wind);
  fprintf(stderr,"displacement = %f\n",displacement);
  fprintf(stderr,"roughness = %f\n",roughness);
  fprintf(stderr,"ref_height = %f\n",ref_height);
  fprintf(stderr,"elevation = %f\n",elevation);
  fprintf(stderr,"b_infilt = %f\n",b_infilt);
  fprintf(stderr,"max_infil = %f\n",max_infil);
  fprintf(stderr,"dt = %f\n",dt);
  fprintf(stderr,"vpd = %f\n",vpd);
  fprintf(stderr,"snow_energy = %f\n",snow_energy);
  fprintf(stderr,"mu = %f\n",mu);
  fprintf(stderr,"rainfall = %f\n",rainfall[0]);
  fprintf(stderr,"Wdew = %f\n",Wdew[0]);
  fprintf(stderr,"grnd_flux = %f\n",grnd_flux[0]);
  fprintf(stderr,"T1 = %f\n",T1[0]);
  fprintf(stderr,"latent_heat = %f\n",latent_heat[0]);
  fprintf(stderr,"sensible_heat = %f\n",sensible_heat[0]);
  fprintf(stderr,"deltaH = %f\n",deltaH[0]);
  fprintf(stderr,"store_error = %f\n",store_error[0]);
  fprintf(stderr,"depth = %f\n",depth[0]);
  fprintf(stderr,"Wcr = %f\n",Wcr[0]);
  fprintf(stderr,"Wpwp = %f\n",Wpwp[0]);
  fprintf(stderr,"Residual Moisture = %f\n",resid_moist[0]);
  fprintf(stderr,"VEG = %i\n",VEG);
  fprintf(stderr,"veg_class = %i\n",veg_class);
  fprintf(stderr,"month = %i\n",month);
  write_layer(layer_wet,iveg,options.Nlayer,depth);
  if(options.DIST_PRCP) 
    write_layer(layer_dry,iveg,options.Nlayer,depth);
  write_vegvar(&(veg_var_wet[0]),iveg);
  if(options.DIST_PRCP) 
    write_vegvar(&(veg_var_dry[0]),iveg);

  if(!options.QUICK_FLUX) {
    fprintf(stderr,"Node\tT\tTnew\tdz\tkappa\tCs\tmoist\tbubble\texpt\tmax_moist\tice\n");
    for(i=0;i<Nnodes;i++) 
      fprintf(stderr,"%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	      i, T_node[i], Tnew_node[i], dz_node[i], kappa_node[i], 
	      Cs_node[i], moist_node[i], bubble_node[i], expt_node[i], 
	      max_moist_node[i], ice_node[i]);
  }

  vicerror("Finished writing calc_surf_energy_bal variables.\nTry increasing SURF_DT to get model to complete cell.\nThen check output for instabilities.\n");

  return(0.0);
    
}

