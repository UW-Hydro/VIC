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
  01-20-00 modified to work with the updated radiation estimation
           routines, as well as the simplified frozen soil moisture
           storage                                              KAC
  01-08-01 sensible heat is now set to 0 W/m^2 when the ground is
           fully covered by snow.                                KAC
  04-12-01 fixed error where sensible heat flux for partial bare
           ground was not multiplied by the snow cover fraction. KAC
  04-29-02 moved calculation of sensible heat so that it is computed
           even in water balance mode.  This assures that it is set
           to 0 in water balance mode and does not yield the 
           cumulative sum of sensible heat from the snowpack.    KAC
  11-18-02 modified to compute the effects of blowing snow on the
           surface energy balance.                               LCB
  10-May-04 Added check that both FS_ACTIVE and FROZEN_SOIL are true
	    before computing *fusion.  This fixes a bug caused when
	    FS_ACTIVE was false and FROZEN_SOIL was true, in which
	    non-zero ice content was calculated at the beginning of
	    the time step but always assumed zero at the point of
	    *fusion computation, leading to *fusion always being
	    calculated as if ice had melted, leading to lower soil
	    temperatures and snow packs that never melted.	TJB

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct  debug;
#endif

  /* define routine input variables */

  /* general model terms */
  int iveg;
  int month;
  int VEG;
  int veg_class;

  double delta_t;

  /* soil layer terms */
  double Cs1;
  double Cs2;
  double D1;
  double D2;
  double T1_old;
  double T2;
  double Ts_old;
  double b_infilt;
  double bubble;
  double dp;
  double expt;
  double ice0;
  double kappa1;
  double kappa2;
  double max_infil;
  double max_moist;
  double moist;

  double *Wcr;
  double *Wpwp;
  double *depth;
  double *resid_moist;

  float *root;

  /* meteorological forcing terms */
  int UnderStory;
  int overstory;

  double NetShortBare;  // net SW that reaches bare ground
  double NetShortGrnd;  // net SW that penetrates snowpack
  double NetShortSnow;  // net SW that reaches snow surface
  double Tair;          // temperature of canopy air or atmosphere
  double atmos_density;
  double atmos_pressure;
  double elevation;
  double emissivity;
  double LongBareIn; // incoming LW to snow-free surface
  double LongSnowIn; // incoming LW to snow surface - if INCLUDE_SNOW
  double mu;
  double surf_atten;
  double vp;
  double vpd;

  double *Wdew;
  double *displacement;
  double *ra;
  double *rainfall;
  double *ref_height;
  double *roughness;
  double *wind;
 
  /* latent heat terms */
  double  Le;

  /* snowpack terms */
  double Advection;
  double OldTSurf;
  double TPack;
  double Tsnow_surf;
  double kappa_snow; // snow conductance / depth
  double melt_energy; // energy consumed in reducing the snowpack coverage
  double snow_coverage; // snowpack coverage fraction
  double snow_density;
/*   double snow_depth; */
  double snow_swq;
  double snow_water;
  int LastSnow;

  double *deltaCC;
  double *refreeze_energy;
  double *VaporMassFlux;
  double *BlowingMassFlux;
  double *SurfaceMassFlux;

  /* soil node terms */
  int     Nnodes;

  double *Cs_node;
  double *T_node;
  double *Tnew_node;
  double *alpha;
  double *beta;
  double *bubble_node;
  double *dz_node;
  double *expt_node;
  double *gamma;
  double *ice_node;
  double *kappa_node;
  double *max_moist_node;
  double *moist_node;

  /* spatial frost terms */
#if SPATIAL_FROST    
  double *frost_fract;
#endif

  /* quick solution frozen soils terms */
#if QUICK_FS
  double **ufwc_table_layer;
  double ***ufwc_table_node;
#endif

  /* model structures */
  layer_data_struct *layer_wet;
  layer_data_struct *layer_dry;
  veg_var_struct *veg_var_wet;
  veg_var_struct *veg_var_dry;

  /* control flags */
  int INCLUDE_SNOW;
  int FS_ACTIVE;
  int NOFLUX;
  int SNOWING;

  int *FIRST_SOLN;

  /* returned energy balance terms */
  double *NetLongBare; // net LW from snow-free ground
  double *NetLongSnow; // net longwave from snow surface - if INCLUDE_SNOW
  double *T1;
  double *deltaH;
  double *fusion;
  double *grnd_flux;
  double *latent_heat;
  double *latent_heat_sub;
  double *sensible_heat;
  double *snow_flux;
  double *store_error;
  double dt;
  double SnowDepth;
  float lag_one;
  float sigma_slope;
  float fetch;
  int Nveg;

  /* Define internal routine variables */
  double Evap;		/** Total evap in m/s **/
  double LongBareOut; // outgoing LW from snow-free ground
  double NetBareRad;
  double TMean;
  double Tmp;
  double error;
  double ice;
/*   double             kappa_snow; */
  double out_long;
  double ra_under;
  double temp_latent_heat;
  double temp_latent_heat_sub;

  /************************************
    Read variables from variable list 
  ************************************/

  /* general model terms */
  iveg                    = (int) va_arg(ap, int);
  month                   = (int) va_arg(ap, int);
  VEG                     = (int) va_arg(ap, int);
  veg_class               = (int) va_arg(ap, int);

  delta_t                 = (double) va_arg(ap, double);

  /* soil layer terms */
  Cs1                     = (double) va_arg(ap, double);
  Cs2                     = (double) va_arg(ap, double);
  D1                      = (double) va_arg(ap, double);
  D2                      = (double) va_arg(ap, double);
  T1_old                  = (double) va_arg(ap, double);
  T2                      = (double) va_arg(ap, double);
  Ts_old                  = (double) va_arg(ap, double);
  b_infilt                = (double) va_arg(ap, double);
  bubble                  = (double) va_arg(ap, double);
  dp                      = (double) va_arg(ap, double);
  expt                    = (double) va_arg(ap, double);
  ice0                    = (double) va_arg(ap, double);
  kappa1                  = (double) va_arg(ap, double);
  kappa2                  = (double) va_arg(ap, double);
  max_infil               = (double) va_arg(ap, double);
  max_moist               = (double) va_arg(ap, double);
  moist                   = (double) va_arg(ap, double);

  Wcr                     = (double *) va_arg(ap, double *);
  Wpwp                    = (double *) va_arg(ap, double *);
  depth                   = (double *) va_arg(ap, double *);
  resid_moist             = (double *) va_arg(ap, double *);

  root                    = (float  *) va_arg(ap, float  *);

  /* meteorological forcing terms */
  UnderStory              = (int) va_arg(ap, int);
  overstory               = (int) va_arg(ap, int);

  NetShortBare            = (double) va_arg(ap, double);
  NetShortGrnd            = (double) va_arg(ap, double);
  NetShortSnow            = (double) va_arg(ap, double);
  Tair                    = (double) va_arg(ap, double);
  atmos_density           = (double) va_arg(ap, double);
  atmos_pressure          = (double) va_arg(ap, double);
  elevation               = (double) va_arg(ap, double);
  emissivity              = (double) va_arg(ap, double);
  LongBareIn              = (double) va_arg(ap, double);
  LongSnowIn              = (double) va_arg(ap, double);
  mu                      = (double) va_arg(ap, double);
  surf_atten              = (double) va_arg(ap, double);
  vp                      = (double) va_arg(ap, double);
  vpd                     = (double) va_arg(ap, double);

  Wdew                    = (double *) va_arg(ap, double *);
  displacement            = (double *) va_arg(ap, double *);
  ra                      = (double *) va_arg(ap, double *);
  rainfall                = (double *) va_arg(ap, double *);
  ref_height              = (double *) va_arg(ap, double *);
  roughness               = (double *) va_arg(ap, double *);
  wind                    = (double *) va_arg(ap, double *);

  /* latent heat terms */
  Le                      = (double) va_arg(ap, double);

  /* snowpack terms */
  Advection               = (double) va_arg(ap, double);
  OldTSurf                = (double) va_arg(ap, double);
  TPack                   = (double) va_arg(ap, double);
  Tsnow_surf              = (double) va_arg(ap, double);
  kappa_snow              = (double) va_arg(ap, double);
  melt_energy             = (double) va_arg(ap, double);
  snow_coverage           = (double) va_arg(ap, double);
  snow_density            = (double) va_arg(ap, double);
  snow_swq                = (double) va_arg(ap, double);
  snow_water              = (double) va_arg(ap, double);
  LastSnow                = (int) va_arg(ap, int);
    
  deltaCC                 = (double *) va_arg(ap, double *);
  refreeze_energy         = (double *) va_arg(ap, double *);
  VaporMassFlux           = (double *) va_arg(ap, double *);
  BlowingMassFlux         = (double *) va_arg(ap, double *);
  SurfaceMassFlux         = (double *) va_arg(ap, double *);

  /* soil node terms */
  Nnodes                  = (int) va_arg(ap, int);

  Cs_node                 = (double *) va_arg(ap, double *);
  T_node                  = (double *) va_arg(ap, double *);
  Tnew_node               = (double *) va_arg(ap, double *);
  alpha                   = (double *) va_arg(ap, double *);
  beta                    = (double *) va_arg(ap, double *);
  bubble_node             = (double *) va_arg(ap, double *);
  dz_node                 = (double *) va_arg(ap, double *);
  expt_node               = (double *) va_arg(ap, double *);
  gamma                   = (double *) va_arg(ap, double *);
  ice_node                = (double *) va_arg(ap, double *);
  kappa_node              = (double *) va_arg(ap, double *);
  max_moist_node          = (double *) va_arg(ap, double *);
  moist_node              = (double *) va_arg(ap, double *);

#if SPATIAL_FROST    
  frost_fract             = (double *) va_arg(ap, double *);
#endif // SPATIAL_FROST

#if QUICK_FS
  ufwc_table_layer        = (double **) va_arg(ap, double **);
  ufwc_table_node         = (double ***) va_arg(ap, double ***);
#endif // QUICK_FS

  /* model structures */
  layer_wet               = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry               = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet             = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry             = (veg_var_struct *) va_arg(ap, veg_var_struct *);

  /* control flags */
  INCLUDE_SNOW            = (int) va_arg(ap, int);
  FS_ACTIVE               = (int) va_arg(ap, int);
  NOFLUX                  = (int) va_arg(ap, int);
  SNOWING                 = (int) va_arg(ap, int);

  FIRST_SOLN              = (int *) va_arg(ap, int *);

  /* returned energy balance terms */
  NetLongBare             = (double *) va_arg(ap, double *);
  NetLongSnow             = (double *) va_arg(ap, double *);
  T1                      = (double *) va_arg(ap, double *);
  deltaH                  = (double *) va_arg(ap, double *);
  fusion                  = (double *) va_arg(ap, double *);
  grnd_flux               = (double *) va_arg(ap, double *);
  latent_heat             = (double *) va_arg(ap, double *);
  latent_heat_sub         = (double *) va_arg(ap, double *);
  sensible_heat           = (double *) va_arg(ap, double *);
  snow_flux               = (double *) va_arg(ap, double *);
  store_error             = (double *) va_arg(ap, double *);
  dt   = (double) va_arg(ap, double);
  SnowDepth  = (double) va_arg(ap, double);
  lag_one = (float) va_arg(ap, float);
  sigma_slope = (float) va_arg(ap, float);
  fetch = (float) va_arg(ap, float);
  Nveg = (int) va_arg(ap, int);

  /***************
    MAIN ROUTINE
  ***************/

  TMean = Ts;
  Tmp = TMean + KELVIN;

  if(options.GRND_FLUX) {
  
    /**********************************************
      Compute Surface Temperature at Half Time Step
    **********************************************/

    if ( snow_coverage > 0 && !INCLUDE_SNOW ) {

      /****************************************
        Compute energy flux through snow pack
      ****************************************/

      *snow_flux = ( kappa_snow * (Tsnow_surf - TMean) );

    }
    else if ( INCLUDE_SNOW ) {
      *snow_flux = 0;
      Tsnow_surf = TMean;
    }
    else *snow_flux = 0;

    /***************************************************************
      Estimate soil temperatures for ground heat flux calculations
    ***************************************************************/

    if ( options.QUICK_FLUX ) {
      /**************************************************************
        Use Liang et al. 1999 Equations to Calculate Ground Heat 
	Flux
      **************************************************************/
      *T1 = estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs1, 
			Cs2, dp, delta_t);
    
      /*****************************************************
	Compute the Ground Heat Flux from the Top Soil Layer
      *****************************************************/
      *grnd_flux = (snow_coverage + (1. - snow_coverage) 
		    * surf_atten) * (kappa1 / D1 * ((*T1) - TMean) 
				     + (kappa2 / D1 * ( 1. - exp( -D1 / dp )) 
					* (T2 - (*T1)))) / 2.;
      
    }
    else {
    /*************************************************************
      Use Finite Difference Method to Explicitly Solve Ground Heat
      Flux at Soil Thermal Nodes (Cherkauer and Lettenmaier, 1999)
    *************************************************************/
      T_node[0] = TMean;
#if QUICK_FS
      solve_T_profile(Tnew_node, T_node, dz_node, kappa_node, Cs_node, 
		      moist_node, delta_t, max_moist_node, bubble_node, 
		      expt_node, ice_node, alpha, beta, gamma, 
		      ufwc_table_node, Nnodes, FIRST_SOLN, FALSE, FS_ACTIVE, 
		      NOFLUX);
#else
      solve_T_profile(Tnew_node, T_node, dz_node, kappa_node, Cs_node, 
		      moist_node, delta_t, max_moist_node, bubble_node, 
		      expt_node, ice_node, alpha, beta, gamma, Nnodes, 
		      FIRST_SOLN, FALSE, FS_ACTIVE, NOFLUX);
#endif
      *T1 = Tnew_node[1];

      /*****************************************************
        Compute the Ground Heat Flux from the Top Soil Layer
      *****************************************************/
      *grnd_flux = (snow_coverage + (1. - snow_coverage) 
		    * surf_atten) * (kappa1 / D1 * ((*T1) - TMean) 
				     + (kappa2 / D2 * (Tnew_node[2] 
						       - (*T1)))) / 2.;
      
    }

    /******************************************************
      Compute the Current Ice Content of the Top Soil Layer
    ******************************************************/
    if((FS_ACTIVE && options.FROZEN_SOIL) && (TMean+ *T1)/2.<0.) {
      ice = moist - maximum_unfrozen_water((TMean+ *T1)/2.,
					   max_moist,bubble,expt);
      if(ice<0.) ice=0.;
    }
    else ice=0.;

    /* compute the change in heat storage */
    *deltaH = (Cs1 * ((Ts_old + T1_old) 
				       - (TMean + *T1)) * D1 / delta_t / 2.);
    
    /* compute the change in heat due to solid - liquid phase changes */
    if (FS_ACTIVE && options.FROZEN_SOIL)
      *fusion = (-ice_density * Lf * (ice0 - ice) * D1 / delta_t);
    
    /* if thin snowpack, compute the change in energy stored in the pack */
    if ( INCLUDE_SNOW ) {
      if ( TMean > 0 )
	*deltaCC = CH_ICE * (snow_swq - snow_water) * (0 - OldTSurf) / delta_t;
      else
	*deltaCC = CH_ICE * (snow_swq - snow_water) * (TMean - OldTSurf) 
	  / delta_t;
      *refreeze_energy = (snow_water * Lf * snow_density) / delta_t;
      *deltaCC *= snow_coverage; // adjust for snow cover fraction
      *refreeze_energy *= snow_coverage; // adjust for snow cover fraction
    }
    
    /** Compute net surface radiation of snow-free area for evaporation 
	estimates **/
    LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp;
    if ( INCLUDE_SNOW ) { // compute net LW at snow surface
      (*NetLongSnow) = (LongSnowIn - snow_coverage 
			* LongBareOut);
    }
    (*NetLongBare)   = (LongBareIn - (1. - snow_coverage) 
			* LongBareOut); // net LW snow-free area
    NetBareRad = (NetShortBare + (*NetLongBare) + *grnd_flux + *deltaH 
		  + *fusion);

  } /* End computation for ground heat flux */
  else { /* ground heat flux not estimated */

    if ( delta_t < SEC_PER_DAY ) {
      
      /** Compute net surface radiation of snow-free area for evaporation 
	  estimates **/
      LongBareOut = (1. - snow_coverage) * STEFAN_B * Tmp * Tmp * Tmp * Tmp;
      if ( INCLUDE_SNOW ) { // compute net LW at snow surface
	(*NetLongSnow) = (LongSnowIn - snow_coverage 
			  * LongBareOut);
      }
      (*NetLongBare)   = (LongBareIn - (1. - snow_coverage) 
			  * LongBareOut); // net LW snow-free area
      NetBareRad = NetShortBare + (*NetLongBare);

    }
    else {

      /** Daily water balance model provides average net shortwave and 
	  longwave radiation **/
      LongBareOut = 
	(1. - snow_coverage) * STEFAN_B * Tmp * Tmp * Tmp * Tmp;
      if ( INCLUDE_SNOW ) {
	(*NetLongSnow) = snow_coverage * LongBareIn / (1. - snow_coverage);
      }
      (*NetLongBare) = LongBareIn;
      NetBareRad = (1. - snow_coverage) * (NetShortBare + (*NetLongBare));

    }
  }

  /** Compute atmospheric stability correction **/
  /** CHECK THAT THIS WORKS FOR ALL SUMMER SITUATIONS **/
  if ( wind[UnderStory] > 0.0 && overstory && SNOWING )
    ra_under = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], 0.f, TMean, Tair, 
			    wind[UnderStory], roughness[UnderStory]);
  else if ( wind[UnderStory] > 0.0 )
    ra_under = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], displacement[UnderStory], 
			    TMean, Tair, wind[UnderStory], 
			    roughness[UnderStory]);
  else
    ra_under = HUGE_RESIST;

  /*************************************************
    Compute Evapotranspiration if not snow covered

    Should evapotranspiration be active when the 
    ground is only paritially covered with snow????
  *************************************************/
  if ( VEG && !SNOWING ) 
    Evap = canopy_evap(layer_wet, layer_dry, veg_var_wet, veg_var_dry, TRUE, 
		       veg_class, month, mu, Wdew, delta_t, NetBareRad, vpd, 
		       NetShortBare, Tair, ra[1], 
		       displacement[1], roughness[1], ref_height[1], 
		       elevation, rainfall, depth, Wcr, Wpwp, 
#if SPATIAL_FROST
		       frost_fract,
#endif // SPATIAL_FROST
		       root);
  else if(!SNOWING) {
    Evap = arno_evap(layer_wet, layer_dry, NetBareRad, Tair, vpd, 
		     NetShortBare, D1, max_moist * depth[0] * 1000., 
		     elevation, b_infilt, displacement[0], 
		     roughness[0], ref_height[0], ra_under, delta_t, mu, 
#if SPATIAL_FROST
		     resid_moist[0], frost_fract);
#else
                     resid_moist[0]);
#endif // SPATIAL_FROST
  }
  else Evap = 0.;
  
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
  **********************************************************************/
  *latent_heat  = -RHO_W * Le * Evap;
  *latent_heat_sub = 0.;

  /** Compute the latent heat flux from a thin snowpack if present **/
  if (INCLUDE_SNOW) {
    latent_heat_from_snow(atmos_density, ice_density, vp, Le, atmos_pressure, 
			  ra_under, TMean, vpd, &temp_latent_heat, 
			  &temp_latent_heat_sub, VaporMassFlux, BlowingMassFlux, SurfaceMassFlux,
			  dt,Tair, LastSnow, snow_water, wind[2],roughness, ref_height[2],
			  SnowDepth, overstory, lag_one, sigma_slope, fetch, Nveg, iveg);
    *latent_heat += temp_latent_heat * snow_coverage;
    *latent_heat_sub = temp_latent_heat_sub * snow_coverage;
  }
  else *latent_heat *= (1. - snow_coverage);

  /************************************************
    Compute the Sensible Heat Flux from the Surface
  ************************************************/
  if ( snow_coverage < 1 || INCLUDE_SNOW ) {
    *sensible_heat = atmos_density * Cp * (Tair - (TMean)) / ra_under;
    if ( !INCLUDE_SNOW ) (*sensible_heat) *= (1. - snow_coverage);
  }
  else *sensible_heat = 0.;

  if(options.GRND_FLUX) {
  
    /*************************************
      Compute Surface Energy Balance Error
    *************************************/
    error = (NetBareRad // net radiation on snow-free area
	     + NetShortGrnd + NetShortSnow // net surface SW 
	     + emissivity * (*NetLongSnow)) // net surface LW 
	     + *sensible_heat // surface sensible heat
	     + (*latent_heat + *latent_heat_sub) // surface latent heats
	     /* heat flux through snowpack - for snow covered fraction */
	     + *snow_flux * snow_coverage
	     /* energy used in reducing snow coverage area */
	     + melt_energy 
	     /* snow energy terms - values are 0 unless INCLUDE_SNOW */
	     + Advection - *deltaCC;
    
    if ( INCLUDE_SNOW ) {
      if (Tsnow_surf == 0.0 && error > -(*refreeze_energy)) {
	*refreeze_energy = -error;
	error = 0.0;
      }
      else {
	error += *refreeze_energy;
      }
    }

    *store_error = error;
  }
  else error = MISSING;

  return error;

}

