#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

#if CLOSE_ENERGY
#define MAX_ITER 250 /* Max number of iterations for total energy balance */
#else
#define MAX_ITER 0   /* No iterations */
#endif // CLOSE_ENERGY
#define GRND_TOL 0.001
#define OVER_TOL 0.001

void surface_fluxes(char                 overstory,
		    double               BareAlbedo,
		    double               height,
		    double               ice0,
		    double               moist,
		    double               mu,
		    double               surf_atten,
		    double              *Melt,
		    double              *Le,
		    double              *aero_resist,
		    double              *baseflow_dry,
		    double              *baseflow_wet,
		    double              *displacement,
		    double              *gauge_correction,
		    double              *inflow_dry,
		    double              *inflow_wet,
		    double              *out_prec,
		    double              *ref_height,
		    double              *roughness,
		    double              *runoff_dry,
		    double              *runoff_wet,
		    double              *snow_inflow,
		    double              *wind,
		    float               *root,
		    int                  Nbands,
		    int                  Ndist,
		    int                  Nlayers,
		    int                  Nveg,
		    int                  band,
		    int                  dp,
		    int                  iveg,
		    int                  rec,
		    int                  veg_class,
		    atmos_data_struct   *atmos,
		    dmy_struct          *dmy,
		    energy_bal_struct   *energy,
		    global_param_struct *gp,
		    layer_data_struct   *layer_dry,
		    layer_data_struct   *layer_wet,
		    snow_data_struct    *snow,
		    soil_con_struct     *soil_con,
		    veg_var_struct      *veg_var_dry,
		    veg_var_struct      *veg_var_wet,
		    float              lag_one,
		    float              sigma_slope,
		    float              fetch)
/**********************************************************************
	surface_fluxes	Keith Cherkauer		February 29, 2000

  Formerly a part of full_energy.c this routine computes all surface
  fluxes, and solves the snow accumulation and ablation algorithm.
  Solutions are for the current snow band and vegetation type (these
  are defined in full_energy before the routine is called).

  modifications:
  10-06-00 modified to handle partial snow cover                KAC
  10-31-00 modified to iterate a solution for the exchange of
           energy between the snowpack and the ground surface.  KAC
  11-18-02 modified to add the effects of blowing snow.         LCB

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif // LINK_DEBUG
 
  double                 total_store_moist[3];
  double                 step_store_moist[3];

  int                    BISECT_OVER;
  int                    BISECT_UNDER;
  int                    INCLUDE_SNOW = FALSE;
  int                    UNSTABLE_CNT;
  int                    UNSTABLE_SNOW = FALSE;
  int                    N_steps;
  int                    UnderStory;
  int                    dist;
  int                    endhour;
  int                    hour;
  int                    lidx;
  int                    over_iter;
  int                    step_dt;
  int                    under_iter;
  double                 Evap;
  double Ls;
  double                 LongUnderIn; // inmoing LW to ground surface
  double                 LongUnderOut; // outgoing LW from ground surface
  double                 NetLongSnow; // net LW over snowpack
  double                 NetShortSnow; // net SW over understory
  double                 NetShortGrnd; // net SW over snow-free surface
  double                 OldTSurf; // previous snow surface temperature
  double                 ShortUnderIn; // incoming SW to understory
  double                 Tair; // air temperature
  double                 Tcanopy; // canopy air temperature
  double                 Tgrnd; // soil surface temperature
  double                 Tsurf; // ground surface temperature
  double                 VPDcanopy; // vapor pressure deficit in canopy/atmos
  double                 VPcanopy; // vapor pressure in canopy/atmos
  double                 coverage; // mid-step snow cover fraction
  double                 delta_coverage; // change in snow cover fraction
  double                 delta_snow_heat; // change in snowpack heat storage
  double                 last_LongUnderOut;
  double                 last_Tcanopy;
  double                 last_Tgrnd;
  double                 last_Tsurf;
  double                 last_latent_ground_heat;
  double                 last_snow_coverage; // previous snow covered area
  double                 last_snow_flux;
  double                 last_tol_under; // previous surface iteration tol
  double                 last_tol_over; // previous overstory iteration tol
  double                 latent_ground_heat; // latent heat from understory
  double                 ppt[2]; // precipitation/melt reaching soil surface
  double                 rainfall[2]; // rainfall
  double                 snowfall[2]; // snowfall
  double                 snow_flux; // heat flux through snowpack
  double                 snow_grnd_flux; // ground heat flux into snowpack
  double                 tol_under;
  double                 tol_over;

  double                 step_Evap;
  double                 step_Wdew;
  double                 step_melt;
  double                 step_melt_energy;  /* energy used to reduce snow coverage */
  double                 step_out_prec;
  double                 step_ppt[2];
  double                 step_prec[2];
  double                 store_AlbedoOver;
  double                 store_AlbedoUnder;
  double                 store_AtmosLatent;
  double                 store_AtmosLatentSub;
  double                 store_AtmosSensible;
  double                 store_LongOverIn;
  double                 store_LongUnderIn;
  double                 store_LongUnderOut;
  double                 store_NetLongAtmos;
  double                 store_NetLongOver;
  double                 store_NetLongUnder;
  double                 store_NetShortAtmos;
  double                 store_NetShortGrnd;
  double                 store_NetShortOver;
  double                 store_NetShortUnder;
  double                 store_ShortOverIn;
  double                 store_ShortUnderIn;
  double                 store_advected_sensible;
  double                 store_canopyevap[2];
  double                 store_canopy_vapor_flux;
  double                 store_layerevap[2][MAX_LAYERS];
  double                 store_melt;
  double                 store_melt_energy;
  double                 store_ppt[2];
  double                 store_sensible;
  double                 store_latent;
  double                 store_latent_sub;
  double                 store_throughfall[2];
  double                 store_vapor_flux;
  double                 store_grnd_flux;
  double                 store_deltaH;
  double                 store_fusion;
  double                 store_advection; 
  double                 store_deltaCC; 
  double                 store_snow_flux; 
  double                 store_refreeze_energy; 
  double                 store_canopy_advection;
  double                 store_canopy_latent;
  double                 store_canopy_latent_sub;
  double                 store_canopy_sensible;
  double                 store_canopy_refreeze;
  layer_data_struct      step_layer[2][MAX_LAYERS];
  energy_bal_struct      snow_energy;
  energy_bal_struct      bare_energy;
  snow_data_struct       step_snow;
  veg_var_struct         snow_veg_var[2];
  veg_var_struct         bare_veg_var[2];

  // handle bisection of understory solution
  double store_tol_under;
  double A_tol_under;
  double B_tol_under;
  double A_snow_flux;
  double B_snow_flux;

  // handle bisection of overstory solution
  double store_tol_over;
  double A_tol_over;
  double B_tol_over;
  double A_Tcanopy;
  double B_Tcanopy;

  /***********************************************************************
    Set temporary variables - preserves original values until iterations
    are completed
  ***********************************************************************/

  energy->advection       = 0; 
  energy->deltaCC         = 0; 
  if ( snow->swq > 0 ) {
      snow_flux = energy->snow_flux;
  }
  else 
    snow_flux = -(energy->grnd_flux + energy->deltaH + energy->fusion);
  coverage                = snow->coverage;
  energy->refreeze_energy = 0; 
  bare_energy.LongUnderOut = energy->LongUnderOut;
  snow_veg_var[WET]       = (*veg_var_wet);
  snow_veg_var[DRY]       = (*veg_var_dry);
  bare_veg_var[WET]       = (*veg_var_wet);
  bare_veg_var[DRY]       = (*veg_var_dry);
  for ( lidx = 0; lidx < Nlayers; lidx++ ) {
    step_layer[WET][lidx] = layer_wet[lidx];
    step_layer[DRY][lidx] = layer_dry[lidx];
  }

  /********************************
    Set-up sub-time step controls
    (May eventually want to set this up so that it is also true 
    if frozen soils are present)
  ********************************/

  if(snow->swq > 0 || snow->snow_canopy > 0 || atmos->snowflag[NR]) {
    hour      = 0;
    endhour   = hour + NF;
    step_dt   = 1;
    //endhour   = hour + NF * options.SNOW_STEP;
    //step_dt   = options.SNOW_STEP;
  }
  else {
    hour      = NR;
    endhour   = hour + gp->dt;
    step_dt   = gp->dt;
  }

  /*******************************************
    Initialize sub-model time step variables
  *******************************************/

  for ( dist = 0; dist < Ndist; dist++ ) {
    store_throughfall[dist] = 0.;
    store_canopyevap[dist]  = 0.;
    for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
      store_layerevap[dist][lidx] = 0.;
    }
  }
  store_canopy_vapor_flux = 0;
  store_vapor_flux        = 0;
  store_ppt[WET]          = 0;
  store_ppt[DRY]          = 0;
  store_melt              = 0;
  step_prec[DRY]          = 0;
  (*snow_inflow)          = 0;
  store_AtmosLatent       = 0;
  store_AtmosLatentSub    = 0;
  store_AtmosSensible     = 0;
  store_LongOverIn        = 0;
  store_LongUnderIn       = 0;
  store_LongUnderOut      = 0;
  store_NetLongAtmos      = 0;
  store_NetLongOver       = 0;
  store_NetLongUnder      = 0;
  store_NetShortAtmos     = 0;
  store_NetShortGrnd      = 0;
  store_NetShortOver      = 0;
  store_NetShortUnder     = 0;
  store_ShortOverIn       = 0;
  store_ShortUnderIn      = 0;
  store_sensible          = 0;
  store_advected_sensible = 0;
  store_latent            = 0;
  store_latent_sub        = 0;
  store_grnd_flux         = 0;
  store_deltaH            = 0;
  store_fusion            = 0;
  store_advection         = 0; 
  store_deltaCC           = 0; 
  store_snow_flux         = 0; 
  store_refreeze_energy   = 0; 
  store_melt_energy       = 0;
  store_canopy_advection  = 0;
  store_canopy_latent     = 0;
  store_canopy_latent_sub = 0;
  store_canopy_sensible   = 0;
  store_canopy_refreeze   = 0;
  store_AlbedoOver        = 0;
  store_AlbedoUnder       = 0;
  N_steps                 = 0;
  last_snow_coverage      = snow->coverage;
  step_Wdew               = veg_var_wet->Wdew;
      
  /*************************
    Compute surface fluxes 
  *************************/

  do {

    /** Solve energy balance for all sub-model time steps **/


    /* set air temperature and precipitation for this snow band */
    Tair = atmos->air_temp[hour] + soil_con->Tfactor[band];
    step_prec[WET] = atmos->prec[hour] / mu * soil_con->Pfactor[band];
    
    // initialize ground surface temperaure
    if ( options.GRND_FLUX ) Tgrnd = energy->T[0];
    else Tgrnd = Tair;

    // initialize canopy terms
    Tcanopy = Tair;
    VPcanopy = atmos->vp[hour];
    VPDcanopy = atmos->vpd[hour];

    over_iter  = 0;
    tol_over   = 999;

    last_Tcanopy      = 999;
    last_LongUnderOut = 999;
    last_snow_flux    = 999;

    // Compute mass flux of blowing snow
    if( !overstory && options.BLOWING && snow->swq > 0.) {
      Ls = (677. - 0.07 * snow->surf_temp) * JOULESPCAL * GRAMSPKG;
      snow->blowing_flux = CalcBlowingSnow((double) step_dt, Tair, 
					   snow->last_snow, snow->surf_water, 
					   wind[2], Ls, atmos->density[hour], 
					   atmos->pressure[hour], 
					   atmos->vp[hour], roughness, 
					   ref_height[2], snow->depth, 
					   lag_one, sigma_slope, 
					   snow->surf_temp, iveg, Nveg, fetch, 
					   veg_lib[iveg].displacement[dmy[rec].month], 
					   veg_lib[iveg].roughness[dmy[rec].month]);
      snow->blowing_flux *= 3600./RHO_W;
    }
    else
      snow->blowing_flux = 0.0;

    // initialize bisection startup
    BISECT_OVER  = FALSE;
    A_tol_over   = 999;
    B_tol_over   = 999;

    do {

      /** Iterate for overstory solution **/

      over_iter++;
      last_tol_over  = tol_over;

      under_iter = 0;
      tol_under  = 999;
      UnderStory = 999;

      UNSTABLE_CNT = 0;
      
      // bisect understory
      BISECT_UNDER = FALSE;
      A_tol_under  = 999;
      B_tol_under  = 999;
      store_tol_under = 999;
      
      A_tol_over   = 999;
      B_tol_over   = 999;
      store_tol_over = 999;

      do {

	/** Iterate for understory solution - itererates to find snow flux **/

	under_iter++;
	last_tol_under = tol_under;

	if ( last_Tcanopy != 999 ) Tcanopy = (last_Tcanopy + Tcanopy) / 2.;
	last_Tcanopy       = Tcanopy;
	A_tol_over         = store_tol_over;
	A_Tcanopy          = Tcanopy;

	last_LongUnderOut  = LongUnderOut;
	LongUnderOut       = bare_energy.LongUnderOut;
	
	// update understory energy balance terms for iteration
	if ( last_snow_flux != 999 ) { 
	  if ( ( fabs(store_tol_under) > fabs(A_tol_under) 
		 && A_tol_under != 999 
		 && fabs(store_tol_under - A_tol_under) > 1. ) 
	       || tol_under < 0 ) { // stepped the correct way
	    UNSTABLE_CNT++;
	    if ( UNSTABLE_CNT > 3 || tol_under < 0 )
	      UNSTABLE_SNOW = TRUE;
	  }
	  else if ( !INCLUDE_SNOW ) { // stepped the wrong way
	    snow_flux = (last_snow_flux + bare_energy.snow_flux) / 2.;
	  } 
	}
	last_snow_flux = snow_flux;
	A_tol_under    = store_tol_under;
	A_snow_flux    = snow_flux;

	snow_grnd_flux = -snow_flux;
	
	// Initialize structures for new iteration
	step_snow      = (*snow);
	snow_energy    = (*energy);
	bare_energy    = (*energy);
	snow_veg_var[WET].Wdew = step_Wdew;
	bare_veg_var[WET].Wdew = step_Wdew;
	snow_veg_var->canopyevap = 0;
	bare_veg_var->canopyevap = 0;
	for ( dist = 0; dist < Ndist; dist ++ ) 
	  for ( lidx = 0; lidx < Nlayers; lidx ++ ) 
	    step_layer[dist][lidx].evap = 0;
	step_snow.canopy_vapor_flux = 0;
	step_snow.vapor_flux = 0;


	/** Solve snow accumulation, ablation and interception **/
	step_melt = solve_snow(overstory, BareAlbedo, LongUnderOut, 
			       gp->MIN_RAIN_TEMP, gp->MAX_SNOW_TEMP, 
			       Tcanopy, Tgrnd, Tair, atmos->density[hour], 
			       dp, ice0, atmos->longwave[hour], moist, mu, 
			       step_prec[WET], atmos->pressure[hour], 
			       atmos->shortwave[hour], snow_grnd_flux, 
			       VPcanopy, VPDcanopy, gp->wind_h, 
			       &energy->AlbedoUnder, &step_Evap, Le, 
			       &LongUnderIn, &NetLongSnow, &NetShortGrnd, 
			       &NetShortSnow, &ShortUnderIn, &OldTSurf, 
			       aero_resist, &coverage, &delta_coverage, 
			       &delta_snow_heat, displacement, 
			       gauge_correction, &step_melt_energy, 
			       &step_out_prec, step_ppt, rainfall, ref_height, 
			       roughness, snow_inflow, snowfall, &surf_atten, 
			       wind, root, UNSTABLE_SNOW, options.Nnode, 
			       Nveg, band, dmy[rec].hour, iveg, 
			       dmy[rec].day_in_year, step_dt, 
			       dmy[rec].month, dmy[rec].day, dmy[rec].year, 
			       rec, veg_class, &UnderStory, &(snow_energy), 
			       step_layer[DRY], step_layer[WET], &(step_snow), 
			       soil_con, &(snow_veg_var[DRY]), 
			       &(snow_veg_var[WET]),
			       lag_one, sigma_slope, fetch);
      
// snow_energy.sensible + snow_energy.latent + snow_energy.latent_sub + NetShortSnow + NetLongSnow + ( snow_grnd_flux + snow_energy.advection - snow_energy.deltaCC + snow_energy.refreeze_energy + snow_energy.advected_sensible ) * snow->coverage

	/* Check that the snow surface temperature was estimated, if not
	   prepare to include thin snowpack in the estimation of the
	   snow-free surface energy balance */
	if ( ( step_snow.surf_temp == 999 || UNSTABLE_SNOW ) 
	     && step_snow.swq > 0 ) {
	  INCLUDE_SNOW = UnderStory + 1;
	  bare_energy.advection = snow_energy.advection;
	  step_snow.surf_temp = snow->surf_temp;
	  step_melt_energy = 0;
	}
	else {
	  INCLUDE_SNOW = FALSE;
	}

	/**************************************************
          Solve Energy Balance Components at Soil Surface
        **************************************************/
	      
	Tsurf = calc_surf_energy_bal((*Le), LongUnderIn, NetLongSnow, 
				     NetShortGrnd, NetShortSnow, OldTSurf, 
				     ShortUnderIn, step_snow.albedo, 
				     snow_energy.latent, 
				     snow_energy.latent_sub, 
				     snow_energy.sensible, 
				     Tcanopy, VPDcanopy, 
				     VPcanopy, snow_energy.advection, 
				     snow->coldcontent, delta_coverage, dp, 
				     ice0, step_melt_energy, moist, mu, 
				     step_snow.coverage, 
				     (snow->depth + step_snow.depth) / 2., 
				     BareAlbedo, surf_atten, 
				     step_snow.vapor_flux, aero_resist, 
				     displacement, &step_melt, step_ppt, 
				     rainfall, ref_height, roughness, 
				     snowfall, wind, root, INCLUDE_SNOW, 
				     UnderStory, options.Nnode, Nveg, band, 
				     step_dt, hour, iveg, options.Nlayer, 
				     (int)overstory, rec, veg_class, atmos, 
				     &(dmy[rec]), &bare_energy, 
				     step_layer[DRY], step_layer[WET], 
				     &(step_snow), soil_con, 
				     &bare_veg_var[DRY], &bare_veg_var[WET], 
				     lag_one, sigma_slope, fetch); 
	if ( INCLUDE_SNOW ) {
	  /* store melt from thin snowpack */
	  step_melt *= 1000.;
	  step_ppt[WET] += step_melt;
	}
	
	/*****************************************
          Compute energy balance with atmosphere
        *****************************************/
	if ( step_snow.snow && overstory && MAX_ITER > 0 ) {
	  // do this if overstory is active and energy balance is closed
	  Tcanopy = calc_atmos_energy_bal(snow_energy.canopy_sensible,
					  bare_energy.sensible, 
					  snow_energy.canopy_latent, 
					  bare_energy.latent, 
					  snow_energy.canopy_latent_sub, 
					  bare_energy.latent_sub, 
					  (*Le),
					  snow_energy.NetLongOver, 
					  bare_energy.NetLongUnder, 
					  snow_energy.NetShortOver, 
					  bare_energy.NetShortUnder, 
					  aero_resist[1], Tair, 
					  atmos->density[hour], 
					  atmos->vp[hour], atmos->vpd[hour], 
					  &bare_energy.AtmosError, 
					  &bare_energy.AtmosLatent,
					  &bare_energy.AtmosLatentSub,
					  &bare_energy.NetLongAtmos, 
					  &bare_energy.NetShortAtmos, 
					  &bare_energy.AtmosSensible, 
					  &VPcanopy, &VPDcanopy);
	  /* iterate to find Tcanopy which will solve the atmospheric energy
	     balance.  Since I do not know vp in the canopy, use the
	     sum of latent heats from the ground and foliage, and iterate
	     on the temperature used for the sensible heat flux from the
	     canopy air to the mixing level */
	}
	else {
	  // else put surface fluxes into atmospheric flux storage so that 
	  // the model will continue to function
	  bare_energy.AtmosLatent    = bare_energy.latent;
	  bare_energy.AtmosLatentSub = bare_energy.latent_sub;
	  bare_energy.AtmosSensible  = bare_energy.sensible;
	  bare_energy.NetLongAtmos   = bare_energy.NetLongUnder;
	  bare_energy.NetShortAtmos  = bare_energy.NetShortUnder;
	}
	bare_energy.Tcanopy = Tcanopy;

	/*****************************************
          Compute iteration tolerance statistics 
        *****************************************/

	// compute understory tolerance
	if ( INCLUDE_SNOW || ( step_snow.swq == 0 && delta_coverage == 0 ) ) {
	  store_tol_under = 0;
	  tol_under       = 0;
	}
	else {
	  store_tol_under = snow_flux - bare_energy.snow_flux;
	  tol_under       = fabs(store_tol_under);
	}
	if ( fabs( tol_under - last_tol_under ) < GRND_TOL && tol_under > 1. )
	  tol_under = -999;

	// compute overstory tolerance
	if ( overstory && step_snow.snow ) {
	  store_tol_over = Tcanopy - last_Tcanopy;
	  tol_over       = fabs( store_tol_over );
	}
	else {
	  store_tol_over = 0;
	  tol_over       = 0;
	}
		
      } while ( ( fabs( tol_under - last_tol_under ) > GRND_TOL )
		&& ( tol_under != 0 ) && (under_iter < MAX_ITER) );

    } while ( ( fabs( tol_over - last_tol_over ) > OVER_TOL 
		&& overstory ) && ( tol_over != 0 ) 
	      && (over_iter < MAX_ITER) );
 
    /**************************************
      Store sub-model time step variables 
    **************************************/

    for(dist = 0; dist < Ndist; dist++) {

      if(iveg != Nveg) {
	if(step_snow.snow) {
	  store_throughfall[dist] += snow_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += snow_veg_var[dist].canopyevap;
	  bare_veg_var[dist].Wdew  = snow_veg_var[dist].Wdew;
	}
	else {
	  store_throughfall[dist] += bare_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += bare_veg_var[dist].canopyevap;
	  snow_veg_var[dist].Wdew  = bare_veg_var[dist].Wdew;
	}
      }

      for(lidx = 0; lidx < options.Nlayer; lidx++)
	store_layerevap[dist][lidx] += step_layer[dist][lidx].evap;

      store_ppt[dist] += step_ppt[dist];
    }

    if(iveg != Nveg) 
      store_canopy_vapor_flux += step_snow.canopy_vapor_flux;
    store_vapor_flux += step_snow.vapor_flux;
      
    store_melt  += step_melt;
    out_prec[0] += step_out_prec * mu;

    if ( INCLUDE_SNOW ) {
      /* copy needed flux terms to the snowpack */
      snow_energy.latent             = bare_energy.latent;
      snow_energy.latent_sub         = bare_energy.latent_sub;
      snow_energy.sensible           = bare_energy.sensible;
      snow_energy.advection          = bare_energy.advection; 
      snow_energy.deltaCC            = bare_energy.deltaCC; 
      snow_energy.snow_flux          = bare_energy.snow_flux; 
      snow_energy.refreeze_energy    = bare_energy.refreeze_energy; 
      snow_energy.advected_sensible  = bare_energy.advected_sensible; 
    }

    store_AtmosLatent       += bare_energy.AtmosLatent;
    store_AtmosLatentSub    += bare_energy.AtmosLatentSub;
    store_AtmosSensible     += bare_energy.AtmosSensible;
    store_LongOverIn        += snow_energy.LongOverIn; 
    store_LongUnderIn       += LongUnderIn; 
    store_LongUnderOut      += bare_energy.LongUnderOut; 
    store_NetShortAtmos     += bare_energy.NetShortAtmos; 
    store_NetShortOver      += snow_energy.NetShortOver; 
    store_NetShortUnder     += bare_energy.NetShortUnder; 
    store_NetLongAtmos      += bare_energy.NetLongAtmos; 
    store_NetLongOver       += snow_energy.NetLongOver; 
    store_NetLongUnder      += bare_energy.NetLongUnder; 
    store_ShortOverIn       += snow_energy.ShortOverIn; 
    store_ShortUnderIn      += bare_energy.ShortUnderIn; 
    store_latent            += bare_energy.latent; 
    store_latent_sub        += bare_energy.latent_sub; 
    store_sensible          += bare_energy.sensible; 
    store_grnd_flux         += bare_energy.grnd_flux; 
    store_deltaH            += bare_energy.deltaH; 
    store_fusion            += bare_energy.fusion; 
    store_canopy_advection  += snow_energy.canopy_advection; 
    store_canopy_latent     += snow_energy.canopy_latent; 
    store_canopy_latent_sub += snow_energy.canopy_latent_sub; 
    store_canopy_sensible   += snow_energy.canopy_sensible; 
    store_canopy_refreeze   += snow_energy.canopy_refreeze; 
    store_AlbedoOver        += snow_energy.AlbedoOver; 
    store_AlbedoUnder       += bare_energy.AlbedoUnder;
    if ( step_snow.swq == 0 && INCLUDE_SNOW ) {
      if ( last_snow_coverage == 0 && step_prec > 0 ) last_snow_coverage = 1;
      store_advection       += snow_energy.advection * last_snow_coverage; 
      store_deltaCC         += snow_energy.deltaCC * last_snow_coverage; 
      store_snow_flux       += bare_energy.snow_flux * last_snow_coverage; 
      store_refreeze_energy += snow_energy.refreeze_energy * last_snow_coverage; 
    }
    else if ( step_snow.snow || INCLUDE_SNOW ) {
      store_advection       += snow_energy.advection * (step_snow.coverage + delta_coverage); 
      store_deltaCC         += snow_energy.deltaCC * (step_snow.coverage + delta_coverage); 
      store_snow_flux       += bare_energy.snow_flux * (step_snow.coverage + delta_coverage); 
      store_refreeze_energy += snow_energy.refreeze_energy * (step_snow.coverage + delta_coverage); 
    }
    store_melt_energy        += step_melt_energy;
    store_advected_sensible  += (step_snow.coverage + delta_coverage) 
      * snow_energy.advected_sensible;

    /* increment time step */
    N_steps ++;
    hour += step_dt;

  } while (hour < endhour);

  /************************************************
    Store snow variables for sub-model time steps 
  ************************************************/

  (*snow) = step_snow;
  snow->vapor_flux = store_vapor_flux;
  snow->canopy_vapor_flux = store_canopy_vapor_flux;
  (*Melt) = store_melt;
  snow->melt = store_melt;
  for(dist = 0; dist < 2; dist++) ppt[dist] = store_ppt[dist];

  /******************************************************
    Store energy flux averages for sub-model time steps 
  ******************************************************/

  (*energy) = bare_energy;
  energy->AtmosLatent       = store_AtmosLatent / (double)N_steps; 
  energy->AtmosLatentSub    = store_AtmosLatentSub / (double)N_steps; 
  energy->AtmosSensible     = store_AtmosSensible / (double)N_steps; 
  energy->NetShortAtmos     = store_NetShortAtmos / (double)N_steps; 
  energy->NetShortOver      = store_NetShortOver / (double)N_steps; 
  energy->NetShortUnder     = store_NetShortUnder / (double)N_steps; 
  energy->NetLongAtmos      = store_NetLongAtmos / (double)N_steps; 
  energy->NetLongOver       = store_NetLongOver / (double)N_steps; 
  energy->NetLongUnder      = store_NetLongUnder / (double)N_steps; 
  energy->ShortOverIn       = store_ShortOverIn / (double)N_steps; 
  energy->ShortUnderIn      = store_ShortUnderIn / (double)N_steps; 
  energy->LongOverIn        = store_LongOverIn / (double)N_steps; 
  energy->LongUnderIn       = store_LongUnderIn / (double)N_steps; 
  energy->latent            = store_latent / (double)N_steps; 
  energy->latent_sub        = store_latent_sub / (double)N_steps; 
  energy->sensible          = store_sensible / (double)N_steps; 
  energy->grnd_flux         = store_grnd_flux / (double)N_steps; 
  energy->deltaH            = store_deltaH / (double)N_steps; 
  energy->fusion            = store_fusion / (double)N_steps; 
  energy->canopy_advection  = store_canopy_advection / (double)N_steps; 
  energy->canopy_latent     = store_canopy_latent / (double)N_steps; 
  energy->canopy_latent_sub = store_canopy_latent_sub / (double)N_steps; 
  energy->canopy_sensible   = store_canopy_sensible / (double)N_steps; 
  energy->canopy_refreeze   = store_canopy_refreeze / (double)N_steps; 
  energy->AlbedoOver        = store_AlbedoOver / (double)N_steps; 
  energy->AlbedoUnder       = store_AlbedoUnder / (double)N_steps;
  if (snow->snow || INCLUDE_SNOW) {
    energy->advection       = store_advection / (double)N_steps; 
    energy->deltaCC         = store_deltaCC / (double)N_steps; 
    energy->snow_flux       = store_snow_flux / (double)N_steps; 
    energy->refreeze_energy = store_refreeze_energy / (double)N_steps; 
  }
  energy->melt_energy       = store_melt_energy / (double)N_steps;
  energy->advected_sensible = store_advected_sensible / (double)N_steps; 
  energy->Tfoliage          = snow_energy.Tfoliage;

// energy->AtmosSensible + energy->AtmosLatent + energy->AtmosLatentSub + energy->NetShortAtmos + energy->NetLongAtmos + energy->grnd_flux + energy->deltaH + energy->fusion + energy->advection - energy->deltaCC + energy->refreeze_energy + energy->advected_sensible

  /**********************************************************
    Store vegetation variable sums for sub-model time steps 
  **********************************************************/

  if(iveg != Nveg) {
    veg_var_wet->throughfall = store_throughfall[WET];
    veg_var_dry->throughfall = store_throughfall[DRY];
    veg_var_wet->canopyevap  = store_canopyevap[WET];
    veg_var_dry->canopyevap  = store_canopyevap[DRY];
    if(snow->snow) {
      veg_var_wet->Wdew        = snow_veg_var[WET].Wdew;
      veg_var_dry->Wdew        = snow_veg_var[DRY].Wdew;
    }
    else {
      veg_var_wet->Wdew        = bare_veg_var[WET].Wdew;
      veg_var_dry->Wdew        = bare_veg_var[DRY].Wdew;
    }
  }

  /**********************************************************
    Store soil layer variables for sub-model time steps 
  **********************************************************/

  for(lidx=0;lidx<Nlayers;lidx++) {
    layer_wet[lidx]      = step_layer[WET][lidx];
    layer_dry[lidx]      = step_layer[DRY][lidx];
    layer_wet[lidx].evap = store_layerevap[WET][lidx];
    layer_dry[lidx].evap = store_layerevap[DRY][lidx];
  }

  /********************************************************
    Compute Runoff, Baseflow, and Soil Moisture Transport
  ********************************************************/

  (*inflow_wet) = ppt[WET];
  (*inflow_dry) = ppt[DRY];

  runoff(layer_wet, layer_dry, energy, soil_con, runoff_wet, runoff_dry, 
	 baseflow_wet, baseflow_dry, ppt, 
#if SPATIAL_FROST
	 soil_con->frost_fract,
#endif // SPATIAL_FROST
	 mu, gp->dt, options.Nnode, band, rec, iveg);

}

#undef MAX_ITER
#undef GRND_TOL
#undef OVER_TOL
