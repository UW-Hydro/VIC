#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

double calc_surf_energy_bal(double             Le,
			    double             LongUnderIn,
			    double             NetLongSnow, // net LW at snow surface
			    double             NetShortGrnd, // net SW transmitted thru snow
			    double             NetShortSnow, // net SW at snow surface
			    double             OldTSurf,
			    double             ShortUnderIn,
			    double             SnowAlbedo,
			    double             SnowLatent,
			    double             SnowLatentSub,
			    double             SnowSensible,
			    double             Tair, // T of canopy or air
			    double             VPDcanopy,
			    double             VPcanopy,
			    double             advection,
			    double             coldcontent,
			    double             delta_coverage, // change in coverage fraction
			    double             dp,
			    double             ice0,
			    double             melt_energy,
			    double             moist,
			    double             snow_coverage,
			    double             snow_depth,
			    double             BareAlbedo,
			    double             surf_atten,
			    double             vapor_flux,
			    double            *aero_resist,
			    double            *aero_resist_used,
			    double            *displacement,
			    double            *melt,
			    double            *ppt,
			    double             rainfall,
			    double            *ref_height,
			    double            *roughness,
			    double             snowfall,
			    double            *wind,
			    float             *root,
			    int                INCLUDE_SNOW,
			    int                UnderStory,
			    int                Nnodes,
			    int                Nveg,
			    int                band,
			    int                dt,
			    int                hour,
			    int                iveg,
			    int                nlayer,
			    int                overstory,
			    int                rec,
			    int                veg_class,
			    double            *CanopLayerBnd,
			    double            *dryFrac,
			    atmos_data_struct *atmos,
			    dmy_struct        *dmy,
			    energy_bal_struct *energy,
			    layer_data_struct *layer,
			    snow_data_struct  *snow,
			    soil_con_struct   *soil_con,
			    veg_var_struct    *veg_var,
			    int                nrecs)
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
    6-8-2000  Modified to make use of spatially distributed 
              soil frost                                   KAC
    03-09-01  Added QUICK_SOLVE options for finite difference
              soil thermal solution.  By iterating on only a
              few choice nodes near the soil surface the 
              simulation time can be significantly reduced
              with minimal additional energy balance errors.  KAC
    11-18-02  Modified to include the effects of blowing snow
              on the surface energy balance.                 LCB
    07-30-03  Made sure that local NOFLUX variable is always set
              to the options flag value.                      KAC
    04-Jun-04 Placed "ERROR" at beginning of screen dump in
	      error_print_surf_energy_bal.				TJB
    16-Jul-04 Cast the last 6 variables in the parameter list
	      passed to root_brent, error_calc_surf_energy_bal
	      and solve_surf_energy_bal as double, since for
	      some reason letting them remain ints or floats
	      caused them to become garbage in the child
	      functions.						TJB
    16-Jul-04 Modified the cap on vapor_flux to re-scale
	      blowing_flux and surface_flux proportionally
	      so that vapor_flux still = their sum.			TJB
    05-Aug-04 Removed lag_one, sigma_slope, and fetch from
	      parameter list, since these were only used in
	      the call to root_brent/func_surf_energy_bal(),
	      which no longer needs them.				TJB
    24-Aug-04 Modified the re-scaling of surface_flux to reduce
	      round-off errors.						TJB
    21-Sep-04 Added ErrorString to store error messages from
	      root_brent.						TJB
    28-Sep-04 Added aero_resist_used to store the aerodynamic
	      resistance used in flux calculations.			TJB
  2007-Apr-06 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.	GCT/KAC
  2007-Apr-24 Changed the read-in order of iveg, and VEG in
	      error_print_surf_energy_bal to be consistent with the
	      call order, also added year, day, and hour to the
	      argument list.						JCA
  2007-Apr-24 Features included for IMPLICIT frozen soils option.	JCA 
	      including:
	        passing in nrecs
	        passing nrec, nrecs, and iveg to func_surf_energy_bal
	        passing bulk_density, soil_density, and quartz to
	          func_surf_energy_bal
  2007-Apr-24 Features included for EXP_TRANS option for frozen soils
	      algorithm.						JCA
  2007-Apr-24 Passing in Zsum_node rather than recalculating.		JCA
  2007-Aug-08 Features included for EXCESS_ICE option for frozen soils
	      algorithm.						JCA
	      including:
	        passing in entire soil_con structure.
  2007-Aug-31 Checked root_brent return value against -998 rather
	      than -9998.						JCA
  2007-Sep-01 Checked for return value of ERROR from
	      solve_surf_energy_bal.					JCA
  2007-Nov-09 Modified code to reset NOFLUX boundary to global option
	      value before computing final soil column temperatures.
	      Previously NOFLUX was set to FALSE for initial QUICK_SOLVE
	      estimates, but never reset to reflect actual bottom
	      boundary conditions for final pass through solver.	KAC
  2009-Feb-09 Removed dz_node from variables passed to
	      func_surf_energy_bal.					KAC via TJB
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.  This
	      allows simulation to continue when energy balance fails
	      to converge by using previous T value.			TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Nov-15 Fixed initialization of Tsurf_fbcount.			TJB
  2009-Nov-15 Changed definitions of D1 and D2 to work for arbitrary
	      node spacing.						TJB
  2009-Dec-11 Replaced "assert" statements with "if" statements.	TJB
  2011-May-24 Replaced call to finish_frozen_soil_calcs() with a call
	      to calc_layer_average_thermal_props(); expanded the set of
	      cases for which this function is called to include
	      FROZEN_SOIL FALSE and QUICK_FLUX TRUE, so that a soil T
	      profile can be estimated and output in these cases.	TJB
  2011-May-31 Removed options.GRND_FLUX.  Now soil temperatures and
	      ground flux are always computed.				TJB
  2011-Aug-09 Now initialize soil thermal properties for all modes of
	      operation.						TJB
  2012-Jan-28 Changed the calculations of ground flux etc for the 
	      exponential node distribution case to be over the same
	      control volume (first soil layer) as for the linear node
	      distribution and quick flux cases.  Now we need to pass
	      entire array of node temperatures to
	      func_surf_energy_bal().					TJB
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
  2013-Jul-25 Added photosynthesis.					TJB
  2013-Dec-27 Moved SPATIAL_SNOW to options_struct.			TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-Apr-25 Added non-climatological LAI.				TJB
  2014-May-05 Added non-climatological vegcover fraction.		TJB
***************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;

  int      FIRST_SOLN[2];
  int      VEG;
  int      i;
  int      nidx;
  int      tmpNnodes;

  double   Cs1;
  double   Cs2;
  double   D1;
  double   D2;
  double   LongBareIn;
  double   NetLongBare;
  double   NetShortBare;
  double   T1;
  double   T1_old;
  double   T2;
  double   Ts_old;
  double   Tsnow_surf;
  double   Tsurf; 
  char     Tsurf_fbflag; 
  int      Tsurf_fbcount; 
  double   albedo;
  double   atmos_density;
  double   atmos_pressure;
  double   atmos_shortwave;
  double   atmos_Catm;
  double   bubble;
  double   delta_t;
  double   emissivity;
  double   error;
  double   expt;
  double   kappa1;
  double   kappa2;
  double   kappa_snow;
  double   max_moist;
  double   ra;
  double   refrozen_water;

  double   Wdew;
  double  *T_node;
  double   Tnew_node[MAX_NODES];
  char     Tnew_fbflag[MAX_NODES];
  int      Tnew_fbcount[MAX_NODES];
  double  *Zsum_node;
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

  double   T_lower, T_upper;
  double   LongSnowIn;
  double   TmpNetLongSnow;
  double   TmpNetShortSnow;
  double   old_swq, old_depth;
  char ErrorString[MAXSTRING];

  /**************************************************
    Set All Variables For Use
  **************************************************/
  /* Initialize T_fbflag */
  Tsurf_fbflag = 0;
  Tsurf_fbcount = 0;
  for (nidx=0; nidx<Nnodes; nidx++) {
    Tnew_fbflag[nidx] = 0;
    Tnew_fbcount[nidx] = 0;
  }

  if(iveg!=Nveg) {
    if(veg_var->vegcover > 0.0) VEG = TRUE;
    else VEG = FALSE;
  }
  else VEG = FALSE;

  // Define control volume for ground flux etc calculations to be first soil layer
  T2                  = soil_con->avg_temp; // soil temperature at very deep depth (>> dp; *NOT* at depth D2)
  Ts_old              = energy->T[0]; // previous surface temperature
  /* Compute previous temperature at boundary between first and second layers */
  if (options.QUICK_FLUX || !options.EXP_TRANS) {
    // T[1] is defined to be the temperature at the boundary between first and second layers
    T1_old              = energy->T[1];
  }
  else {
    // need to interpolate to find temperature at the boundary between first and second layers
    i=0;
    while (soil_con->Zsum_node[i] < soil_con->depth[0]) { i++; }
    T1_old              = energy->T[i-1] + (energy->T[i]-energy->T[i-1])*(soil_con->depth[0]-soil_con->Zsum_node[i-1])/(soil_con->Zsum_node[i]-soil_con->Zsum_node[i-1]); // i should be greater than 0, else code would have aborted in initialize_model_state()
  }
  D1                  = soil_con->depth[0]; // top layer thickness
  D2                  = soil_con->depth[0]; // Distance below layer boundary to consider (= top layer thickness)
  kappa1              = energy->kappa[0]; // top layer conductivity
  kappa2              = energy->kappa[1]; // second layer conductivity
  Cs1                 = energy->Cs[0]; // top layer heat capacity
  Cs2                 = energy->Cs[1]; // second layer heat capacity
  atmos_density       = atmos->density[hour]; // atmospheric density
  atmos_pressure      = atmos->pressure[hour]; // atmospheric pressure
  atmos_shortwave     = atmos->shortwave[hour]; // incoming shortwave radiation
  atmos_Catm          = atmos->Catm[hour]; // CO2 mixing ratio
  emissivity          = 1.; // longwave emissivity
  delta_t             = (double)dt * 3600.;
  max_moist           = soil_con->max_moist[0] / (soil_con->depth[0]*1000.);
  bubble              = soil_con->bubble[0];
  expt                = soil_con->expt[0];
  Tsnow_surf          = snow->surf_temp;
  Wdew                = veg_var->Wdew;
  FIRST_SOLN[0] = TRUE;
  FIRST_SOLN[1] = TRUE;
  if ( snow->depth > 0. ) 
    kappa_snow = K_SNOW * (snow->density) * (snow->density) / snow_depth;
  else 
    kappa_snow = 0;

  /** compute incoming and net radiation **/
  NetShortBare  = ( ShortUnderIn * (1. - ( snow_coverage + delta_coverage ) ) 
		    * (1. - BareAlbedo) + ShortUnderIn * ( delta_coverage ) 
		    * ( 1. - SnowAlbedo ) );
  LongBareIn    = (1. - snow_coverage ) * LongUnderIn;
  if ( INCLUDE_SNOW || snow->swq == 0 ) { 
    TmpNetLongSnow  = NetLongSnow;
    TmpNetShortSnow = NetShortSnow;
    LongSnowIn      = snow_coverage * LongUnderIn;
  }
  else {
    TmpNetShortSnow = 0.; 
    TmpNetLongSnow  = 0.; 
    LongSnowIn      = 0.;
  }

  /*************************************************************
    Prepare soil node variables for finite difference solution
  *************************************************************/

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
  Zsum_node      = soil_con->Zsum_node;   
  ice_node       = energy->ice;

  /**************************************************
    Find Surface Temperature Using Root Brent Method
  **************************************************/
  if(options.FULL_ENERGY) {

    /** If snow included in solution, temperature cannot exceed 0C  **/
    if ( INCLUDE_SNOW ) {
      T_lower = energy->T[0]-SURF_DT;
      T_upper = 0.;
    }
    else {
      T_lower = 0.5*(energy->T[0]+Tair)-SURF_DT;
      T_upper = 0.5*(energy->T[0]+Tair)+SURF_DT;
    }

    if ( options.QUICK_SOLVE && !options.QUICK_FLUX ) {
      // Set iterative Nnodes using the depth of the thaw layer
      tmpNnodes = 0;
      for ( nidx = Nnodes-5; nidx >= 0; nidx-- ) 
	if ( T_node[nidx] >= 0 && T_node[nidx+1] < 0 ) tmpNnodes = nidx+1;
      if ( tmpNnodes == 0 ) { 
	if ( T_node[0] <= 0 && T_node[1] >= 0 ) tmpNnodes = Nnodes;
	else tmpNnodes = 3;
      }
      else tmpNnodes += 4;
    }
    else { 
      tmpNnodes = Nnodes;
    }

    Tsurf = root_brent(T_lower, T_upper, ErrorString, func_surf_energy_bal,
		       rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t,
		       Cs1, Cs2, D1, D2,
		       T1_old, T2, Ts_old, energy->T, bubble, dp, 
		       expt, ice0, kappa1, kappa2, 
		       max_moist, moist, root, CanopLayerBnd,
		       UnderStory, overstory, NetShortBare, NetShortGrnd, 
		       TmpNetShortSnow, Tair, atmos_density, 
		       atmos_pressure,
		       emissivity, LongBareIn, LongSnowIn, surf_atten, 
		       VPcanopy, VPDcanopy, 
		       atmos_shortwave, atmos_Catm, dryFrac,
		       &Wdew, displacement, aero_resist, aero_resist_used, 
		       rainfall, ref_height, roughness, wind, Le, 
		       energy->advection, OldTSurf, snow->pack_temp, 
		       Tsnow_surf, kappa_snow, melt_energy, 
		       snow_coverage, 
		       snow->density, snow->swq, snow->surf_water,
		       &energy->deltaCC, &energy->refreeze_energy, 
		       &snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
		       tmpNnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
		       alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
		       ice_node, kappa_node, max_moist_node, moist_node, 
		       soil_con, layer, veg_var, 
		       INCLUDE_SNOW, options.NOFLUX, options.EXP_TRANS, snow->snow, 
		       FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
		       &energy->deltaH, &energy->fusion, &energy->grnd_flux, 
		       &energy->latent, &energy->latent_sub, 
		       &energy->sensible, &energy->snow_flux, &energy->error);
 
    if(Tsurf <= -998 ) {  
      if (options.TFALLBACK) {
        Tsurf = Ts_old;
        Tsurf_fbflag = 1;
        Tsurf_fbcount++;
      }
      else {
        fprintf(stderr, "SURF_DT = %.2f\n", SURF_DT);
        error = error_calc_surf_energy_bal(Tsurf, dmy->year, dmy->month, dmy->day, dmy->hour, VEG, iveg,
					   veg_class, delta_t, Cs1, Cs2, D1, D2, 
					   T1_old, T2, Ts_old, energy->T,
					   soil_con->b_infilt, bubble, dp, 
					   expt, ice0, kappa1, kappa2, 
					   soil_con->max_infil, max_moist, 
					   moist, soil_con->Wcr, soil_con->Wpwp, 
					   soil_con->depth, 
					   soil_con->resid_moist, root, CanopLayerBnd,
					   UnderStory, overstory, NetShortBare, 
					   NetShortGrnd, TmpNetShortSnow, Tair, 
					   atmos_density, atmos_pressure, 
					   (double)soil_con->elevation, 
					   emissivity, LongBareIn, LongSnowIn, 
					   surf_atten, VPcanopy, VPDcanopy, 
					   atmos_shortwave, atmos_Catm, dryFrac,
					   &Wdew, displacement, aero_resist, aero_resist_used, 
					   rainfall, ref_height, roughness, 
					   wind, Le, energy->advection, 
					   OldTSurf, snow->pack_temp, 
					   Tsnow_surf, 
					   kappa_snow, melt_energy, 
					   snow_coverage, snow->density, 
					   snow->swq, snow->surf_water, 
					   &energy->deltaCC, 
					   &energy->refreeze_energy, 
					   &snow->vapor_flux, Nnodes, Cs_node, 
					   T_node, Tnew_node, alpha, beta, 
					   bubble_node, Zsum_node, expt_node, 
					   gamma, ice_node, kappa_node, 
					   max_moist_node, moist_node, 
					   soil_con->frost_fract, 
					   layer, veg_var, 
					   INCLUDE_SNOW, soil_con->FS_ACTIVE, 
					   options.NOFLUX, options.EXP_TRANS,
					   snow->snow, FIRST_SOLN, &NetLongBare, 
					   &TmpNetLongSnow, &T1, &energy->deltaH, 
					   &energy->fusion, &energy->grnd_flux, 
					   &energy->latent, 
					   &energy->latent_sub, 
					   &energy->sensible, 
					   &energy->snow_flux, &energy->error, ErrorString);
        return ( ERROR );
      }
    }

    /**************************************************
      Recalculate Energy Balance Terms for Final Surface Temperature
    **************************************************/

    if ( Ts_old * Tsurf < 0 && options.QUICK_SOLVE ) {
      tmpNnodes = Nnodes;
      FIRST_SOLN[0] = TRUE;
      
      Tsurf = root_brent(T_lower, T_upper, ErrorString, func_surf_energy_bal,
			 rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t,
			 Cs1, Cs2, D1, D2, 
			 T1_old, T2, Ts_old, energy->T, bubble, dp, 
			 expt, ice0, kappa1, kappa2, 
			 max_moist, moist, root, CanopLayerBnd,
			 UnderStory, overstory, NetShortBare, NetShortGrnd, 
			 TmpNetShortSnow, Tair, atmos_density, atmos_pressure, 
			 emissivity, LongBareIn, LongSnowIn, surf_atten, 
			 VPcanopy, VPDcanopy, 
		         atmos_shortwave, atmos_Catm, dryFrac,
			 &Wdew, displacement, aero_resist, aero_resist_used, 
			 rainfall, ref_height, roughness, wind, Le, 
			 energy->advection, OldTSurf, snow->pack_temp, 
			 Tsnow_surf, kappa_snow, melt_energy, 
			 snow_coverage, 
			 snow->density, snow->swq, snow->surf_water,
			 &energy->deltaCC, &energy->refreeze_energy, 
			 &snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
			 tmpNnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
			 alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
			 ice_node, kappa_node, max_moist_node, moist_node, 
			 soil_con, layer, veg_var, 
			 INCLUDE_SNOW, options.NOFLUX, options.EXP_TRANS, snow->snow, 
			 FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
			 &energy->deltaH, &energy->fusion, &energy->grnd_flux, 
			 &energy->latent, &energy->latent_sub, 
			 &energy->sensible, &energy->snow_flux, &energy->error);
      
      if(Tsurf <=  -998 ) {  
        if (options.TFALLBACK) {
          Tsurf = Ts_old;
          Tsurf_fbflag = 1;
          Tsurf_fbcount++;
        }
        else {
	  error = error_calc_surf_energy_bal(Tsurf, dmy->year, dmy->month, dmy->day, dmy->hour, VEG, iveg,
					     veg_class, delta_t, Cs1, Cs2, D1, 
					     D2, T1_old, T2, Ts_old, energy->T,
					     soil_con->b_infilt, bubble, dp, 
					     expt, ice0, kappa1, kappa2, 
					     soil_con->max_infil, max_moist, 
					     moist, soil_con->Wcr, 
					     soil_con->Wpwp, soil_con->depth, 
					     soil_con->resid_moist, root, CanopLayerBnd,
					     UnderStory, overstory, 
					     NetShortBare, NetShortGrnd, 
					     TmpNetShortSnow, Tair, 
					     atmos_density, atmos_pressure, 
					     (double)soil_con->elevation, 
					     emissivity, LongBareIn, LongSnowIn, 
					     surf_atten, VPcanopy, VPDcanopy,
					     atmos_shortwave, atmos_Catm, dryFrac,
					     &Wdew, displacement, 
					     aero_resist, aero_resist_used, rainfall, ref_height, 
					     roughness, wind, Le, 
					     energy->advection, 
					     OldTSurf, snow->pack_temp, 
					     Tsnow_surf, 
					     kappa_snow, melt_energy, 
					     snow_coverage, snow->density, 
					     snow->swq, snow->surf_water, 
					     &energy->deltaCC, 
					     &energy->refreeze_energy, 
					     &snow->vapor_flux, Nnodes, Cs_node, 
					     T_node, Tnew_node, alpha, beta, 
					     bubble_node, Zsum_node, expt_node, 
					     gamma, ice_node, kappa_node, 
					     max_moist_node, moist_node, 
					     soil_con->frost_fract, 
					     layer, veg_var, INCLUDE_SNOW, 
					     soil_con->FS_ACTIVE, options.NOFLUX, options.EXP_TRANS,
					     snow->snow, FIRST_SOLN, &NetLongBare, 
					     &TmpNetLongSnow, &T1, 
					     &energy->deltaH, &energy->fusion, 
					     &energy->grnd_flux, 
					     &energy->latent, 
					     &energy->latent_sub, 
					     &energy->sensible, 
					     &energy->snow_flux, &energy->error, ErrorString);
          return ( ERROR );
        }
      }
    }
  }
  else {

    /** Frozen soil model run with no surface energy balance **/
    Tsurf  = Tair;

  }
  
  if ( options.QUICK_SOLVE && !options.QUICK_FLUX ) 
    // Reset model so that it solves thermal fluxes for full soil column
    FIRST_SOLN[0] = TRUE;
  
  error = solve_surf_energy_bal(Tsurf, 
				rec, nrecs, dmy->month, VEG, veg_class, iveg, delta_t,
			        Cs1, Cs2, D1, D2, 
				T1_old, T2, Ts_old, energy->T, bubble, dp, 
				expt, ice0, kappa1, kappa2, 
				max_moist, moist, root, CanopLayerBnd,
				UnderStory, overstory, NetShortBare, NetShortGrnd, 
				TmpNetShortSnow, Tair, atmos_density, atmos_pressure,
				emissivity, LongBareIn, LongSnowIn, surf_atten, 
				VPcanopy, VPDcanopy, 
				atmos_shortwave, atmos_Catm, dryFrac,
				&Wdew, displacement, aero_resist, aero_resist_used, 
				rainfall, ref_height, roughness, wind, Le, 
				energy->advection, OldTSurf, snow->pack_temp, 
				Tsnow_surf, kappa_snow, melt_energy, 
				snow_coverage, 
				snow->density, snow->swq, snow->surf_water,
				&energy->deltaCC, &energy->refreeze_energy, 
				&snow->vapor_flux, &snow->blowing_flux, &snow->surface_flux,
				Nnodes, Cs_node, T_node, Tnew_node, Tnew_fbflag, Tnew_fbcount,
				alpha, beta, bubble_node, Zsum_node, expt_node, gamma, 
				ice_node, kappa_node, max_moist_node, moist_node, 
				soil_con, layer, veg_var, 
				INCLUDE_SNOW, options.NOFLUX, options.EXP_TRANS, snow->snow, 
				FIRST_SOLN, &NetLongBare, &TmpNetLongSnow, &T1, 
				&energy->deltaH, &energy->fusion, &energy->grnd_flux, 
				&energy->latent, &energy->latent_sub, 
				&energy->sensible, &energy->snow_flux, &energy->error);
  if(error == ERROR)
    return(ERROR);
  else
    energy->error = error;
  
  /***************************************************
    Recalculate Soil Moisture and Thermal Properties
  ***************************************************/
    if(options.QUICK_FLUX) {

      Tnew_node[0] = Tsurf;
      Tnew_node[1] = T1;
      Tnew_node[2] = soil_con->avg_temp + (T1-soil_con->avg_temp)*exp(-(soil_con->Zsum_node[2]-D1)/dp);

    }
    calc_layer_average_thermal_props(energy, layer, soil_con, Nnodes, iveg, Tnew_node);

  /** Store precipitation that reaches the surface */
  if ( !snow->snow && !INCLUDE_SNOW ) {
    if ( iveg != Nveg ) {
      *ppt = veg_var->throughfall;
    }
    else {
      *ppt = rainfall;
    }
  }

  /****************************************
    Store understory energy balance terms 
  ****************************************/

// energy->sensible + energy->latent + energy->latent_sub + NetShortBare + NetLongBare + energy->grnd_flux + energy->deltaH + energy->fusion + energy->snow_flux

  energy->NetShortGrnd = NetShortGrnd;
  if ( INCLUDE_SNOW ) {
    energy->NetLongUnder  = NetLongBare + TmpNetLongSnow;
    energy->NetShortUnder = NetShortBare + TmpNetShortSnow + NetShortGrnd;
    energy->latent        = (energy->latent);
    energy->latent_sub    = (energy->latent_sub);
    energy->sensible      = (energy->sensible);
  }
  else {
    energy->NetLongUnder  = NetLongBare + NetLongSnow;
    energy->NetShortUnder = NetShortBare + NetShortSnow + NetShortGrnd;
/*     energy->latent        = (SnowLatent + (1. - snow_coverage)  */
/* 			     * energy->latent); */
/*     energy->latent_sub    = (SnowLatentSub  */
/* 			     + (1. - snow_coverage) * energy->latent_sub); */
/*     energy->sensible      = (SnowSensible  */
/* 			     + (1. - snow_coverage) * energy->sensible); */
    energy->latent        = (SnowLatent + energy->latent);
    energy->latent_sub    = (SnowLatentSub + energy->latent_sub);
    energy->sensible      = (SnowSensible + energy->sensible);
  }
  energy->LongUnderOut  = LongUnderIn - energy->NetLongUnder;
  energy->AlbedoUnder   = ((1. - ( snow_coverage + delta_coverage ) ) 
			   * BareAlbedo + ( snow_coverage + delta_coverage ) 
			   * SnowAlbedo );
  energy->melt_energy   = melt_energy;
  energy->Tsurf         = (snow->coverage * snow->surf_temp 
			   + (1. - snow->coverage) * Tsurf);

  /*********************************************************************
    adjust snow water balance for vapor mass flux if snowpack included 
  *********************************************************************/

//NEED TO ADJUST SNOW COVERAGE FRACTION - AND STORAGE

  if ( INCLUDE_SNOW ) {

    // don't allow vapor_flux to exceed snow pack
    if (-(snow->vapor_flux) > snow->swq) {
      // if vapor_flux exceeds snow pack, we not only need to
      // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
      snow->blowing_flux *= -( snow->swq / snow->vapor_flux );
//      snow->surface_flux *= -( snow->swq / snow->vapor_flux );
      snow->vapor_flux = -(snow->swq);
      snow->surface_flux = snow->vapor_flux - snow->blowing_flux;
    }

    /* adjust snowpack for vapor flux */
    old_swq           = snow->swq;
    snow->swq        += snow->vapor_flux;
    snow->surf_water += snow->vapor_flux;
    snow->surf_water  = ( snow->surf_water < 0 ) ? 0. : snow->surf_water;

    /* compute snowpack melt or refreeze */
    if (energy->refreeze_energy >= 0.0) {
      refrozen_water = energy->refreeze_energy / ( Lf * RHO_W ) * delta_t; 
      if ( refrozen_water > snow->surf_water) {
        refrozen_water = snow->surf_water;
        energy->refreeze_energy = refrozen_water * Lf * RHO_W / delta_t;
      } 
      snow->surf_water -= refrozen_water;
      if (snow->surf_water < 0.0) snow->surf_water = 0.0;
      (*melt)           = 0.0;

    }
    else {
      
      /* Calculate snow melt */      
      (*melt) = fabs(energy->refreeze_energy) / (Lf * RHO_W) * delta_t;
      snow->swq -= *melt;
      if ( snow->swq < 0 ) { 
	(*melt) += snow->swq;
	snow->swq = 0;
      }
    }

    if ( snow->swq > 0 ) {

      // set snow energy terms
      snow->surf_temp   = ( Tsurf > 0 ) ? 0 : Tsurf;
      snow->coldcontent = CH_ICE * snow->surf_temp * snow->swq;

      // recompute snow depth
      old_depth   = snow->depth;
      snow->depth = 1000. * snow->swq / snow->density; 
      
      /** Check for Thin Snowpack which only Partially Covers Grid Cell
	  exists only if not snowing and snowpack has started to melt **/
      if (options.SPATIAL_SNOW) {
        snow->coverage = calc_snow_coverage(&snow->store_snow, 
					    soil_con->max_snow_distrib_slope, 
					    snow_coverage, snow->swq,
					    old_swq, snow->depth, old_depth, 
					    (*melt) - snow->vapor_flux, 
					    &snow->max_snow_depth, snowfall, 
					    &snow->store_swq, 
					    &snow->snow_distrib_slope,
					    &snow->store_coverage);
      }
      else {
        if ( snow->swq > 0 ) snow->coverage = 1.;
        else snow->coverage = 0.;
      }

      if ( snow->surf_temp > 0 ) 
	energy->snow_flux = ( energy->grnd_flux + energy->deltaH 
			      + energy->fusion );

    }
    else {
      /* snowpack melts completely */
      snow->density    = 0.;
      snow->depth      = 0.;
      snow->surf_water = 0;
      snow->pack_water = 0;
      snow->surf_temp  = 0;
      snow->pack_temp  = 0;
      snow->coverage   = 0;
      if (options.SPATIAL_SNOW)
        snow->store_swq = 0.;
    }
    snow->vapor_flux *= -1;
  }

  /** record T flag values **/
  energy->Tsurf_fbflag = Tsurf_fbflag;
  energy->Tsurf_fbcount += Tsurf_fbcount;
  for (nidx=0; nidx<Nnodes; nidx++) {
    energy->T_fbflag[nidx] = Tnew_fbflag[nidx];
    energy->T_fbcount[nidx] += Tnew_fbcount[nidx];
  }

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
/**********************************************************************
  Modifications:
  2009-Mar-03 Fixed format string for print statement, eliminates
	      compiler WARNING.						KAC via TJB
  2012-Jan-28 Added Told_node array.					TJB
**********************************************************************/

  extern option_struct options;

  /* Define imported variables */

  /* general model terms */
  int year,month,day,hour;
  int iveg;
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
  double *Told_node;
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
  double *CanopLayerBnd;

  /* meteorological forcing terms */
  int UnderStory;
  int overstory;

  double NetShortBare;  // net SW that reaches bare ground
  double NetShortGrnd;  // net SW that penetrates snowpack
  double NetShortSnow;  // net SW that reaches snow surface
  double Tair;
  double atmos_density;
  double atmos_pressure;
  double elevation;
  double emissivity;
  double LongBareIn; 
  double LongSnowIn; 
  double surf_atten;
  double vp;
  double vpd;
  double atmos_shortwave;
  double atmos_Catm;
  double dryFrac;

  double *Wdew;
  double *displacement;
  double *ra;
  double *ra_used;
  double rainfall;
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
  double kappa_snow;
  double melt_energy;
  double snow_coverage;
  double snow_density;
  double snow_swq;
  double snow_water;

  double *deltaCC;
  double *refreeze_energy;
  double *VaporMassFlux;

  /* soil node terms */
  int Nnodes;

  double *Cs_node;
  double *T_node;
  double *Tnew_node;
  double *alpha;
  double *beta;
  double *bubble_node;
  double *Zsum_node;
  double *expt_node;
  double *gamma;
  double *ice_node;
  double *kappa_node;
  double *max_moist_node;
  double *moist_node;

  /* spatial frost terms */
  double *frost_fract;

  /* model structures */
  layer_data_struct *layer;
  veg_var_struct *veg_var;

  /* control flags */
  int INCLUDE_SNOW;
  int FS_ACTIVE;
  int NOFLUX;
  int EXP_TRANS;

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

  char *ErrorString;

  /* Define internal routine variables */
  int                i;

  /***************************
    Read Variables from List
  ***************************/

  /* general model terms */
  year                    = (int) va_arg(ap, int);
  month                   = (int) va_arg(ap, int);
  day                     = (int) va_arg(ap, int);
  hour                    = (int) va_arg(ap, int);
  VEG                     = (int) va_arg(ap, int);
  iveg                    = (int) va_arg(ap, int);
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
  Told_node               = (double *) va_arg(ap, double *);
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
  CanopLayerBnd           = (double *) va_arg(ap, double *);

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
  surf_atten              = (double) va_arg(ap, double);
  vp                      = (double) va_arg(ap, double);
  vpd                     = (double) va_arg(ap, double);
  atmos_shortwave         = (double) va_arg(ap, double);
  atmos_Catm              = (double) va_arg(ap, double);
  dryFrac                 = (double) va_arg(ap, double);

  Wdew                    = (double *) va_arg(ap, double *);
  displacement            = (double *) va_arg(ap, double *);
  ra                      = (double *) va_arg(ap, double *);
  ra_used                 = (double *) va_arg(ap, double *);
  rainfall                = (double) va_arg(ap, double);
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

  deltaCC                 = (double *) va_arg(ap, double *);
  refreeze_energy         = (double *) va_arg(ap, double *);
  VaporMassFlux           = (double *) va_arg(ap, double *);

  /* soil node terms */
  Nnodes                  = (int) va_arg(ap, int);

  Cs_node                 = (double *) va_arg(ap, double *);
  T_node                  = (double *) va_arg(ap, double *);
  Tnew_node               = (double *) va_arg(ap, double *);
  alpha                   = (double *) va_arg(ap, double *);
  beta                    = (double *) va_arg(ap, double *);
  bubble_node             = (double *) va_arg(ap, double *);
  Zsum_node               = (double *) va_arg(ap, double *);
  expt_node               = (double *) va_arg(ap, double *);
  gamma                   = (double *) va_arg(ap, double *);
  ice_node                = (double *) va_arg(ap, double *);
  kappa_node              = (double *) va_arg(ap, double *);
  max_moist_node          = (double *) va_arg(ap, double *);
  moist_node              = (double *) va_arg(ap, double *);
  frost_fract             = (double *) va_arg(ap, double *);

  /* model structures */
  layer               = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var             = (veg_var_struct *) va_arg(ap, veg_var_struct *);

  /* control flags */
  INCLUDE_SNOW            = (int) va_arg(ap, int);
  FS_ACTIVE               = (int) va_arg(ap, int);
  NOFLUX                  = (int) va_arg(ap, int);
  EXP_TRANS               = (int) va_arg(ap, int);
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

  ErrorString             = (char *) va_arg(ap, char *);

  /***************
    Main Routine
  ***************/

  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: calc_surf_energy_bal failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  /* Print Variables */
  /* general model terms */
  fprintf(stderr, "iveg = %i\n", iveg);
  fprintf(stderr, "year = %i\n", year);
  fprintf(stderr, "month = %i\n", month);
  fprintf(stderr, "day = %i\n", day);
  fprintf(stderr, "hour = %i\n", hour);
  fprintf(stderr, "VEG = %i\n", VEG);
  fprintf(stderr, "veg_class = %i\n", veg_class);
  fprintf(stderr, "delta_t = %f\n",  delta_t);

  /* soil layer terms */
  fprintf(stderr, "Cs1 = %f\n",  Cs1);
  fprintf(stderr, "Cs2 = %f\n",  Cs2);
  fprintf(stderr, "D1 = %f\n",  D1);
  fprintf(stderr, "D2 = %f\n",  D2);
  fprintf(stderr, "T1_old = %f\n",  T1_old);
  fprintf(stderr, "T2 = %f\n",  T2);
  fprintf(stderr, "Ts_old = %f\n",  Ts_old);
  fprintf(stderr, "b_infilt = %f\n",  b_infilt);
  fprintf(stderr, "bubble = %f\n",  bubble);
  fprintf(stderr, "dp = %f\n",  dp);
  fprintf(stderr, "expt = %f\n",  expt);
  fprintf(stderr, "ice0 = %f\n",  ice0);
  fprintf(stderr, "kappa1 = %f\n",  kappa1);
  fprintf(stderr, "kappa2 = %f\n",  kappa2);
  fprintf(stderr, "max_infil = %f\n",  max_infil);
  fprintf(stderr, "max_moist = %f\n",  max_moist);
  fprintf(stderr, "moist = %f\n",  moist);

  fprintf(stderr, "*Wcr = %f\n",  *Wcr);
  fprintf(stderr, "*Wpwp = %f\n",  *Wpwp);
  fprintf(stderr, "*depth = %f\n",  *depth);
  fprintf(stderr, "*resid_moist = %f\n",  *resid_moist);

  fprintf(stderr, "*root = %f\n",  *root);

  /* meteorological forcing terms */
  fprintf(stderr, "UnderStory = %i\n", UnderStory);
  fprintf(stderr, "overstory = %i\n", overstory);

  fprintf(stderr, "NetShortBare = %f\n",  NetShortBare); 
  fprintf(stderr, "NetShortGrnd = %f\n",  NetShortGrnd); 
  fprintf(stderr, "NetShortSnow = %f\n",  NetShortSnow); 
  fprintf(stderr, "Tair = %f\n",  Tair);
  fprintf(stderr, "atmos_density = %f\n",  atmos_density);
  fprintf(stderr, "atmos_pressure = %f\n",  atmos_pressure);
  fprintf(stderr, "elevation = %f\n",  elevation);
  fprintf(stderr, "emissivity = %f\n",  emissivity);
  fprintf(stderr, "LongBareIn = %f\n",  LongBareIn); 
  fprintf(stderr, "LongSnowIn = %f\n",  LongSnowIn); 
  fprintf(stderr, "surf_atten = %f\n",  surf_atten);
  fprintf(stderr, "vp = %f\n",  vp);
  fprintf(stderr, "vpd = %f\n",  vpd);
  fprintf(stderr, "atmos_shortwave = %f\n",  atmos_shortwave);
  fprintf(stderr, "atmos_Catm = %f\n",  atmos_Catm);
  fprintf(stderr, "dryFrac = %f\n",  dryFrac);

  fprintf(stderr, "*Wdew = %f\n",  *Wdew);
  fprintf(stderr, "*displacement = %f\n",  *displacement);
  fprintf(stderr, "*ra = %f\n",  *ra);
  fprintf(stderr, "*ra_used = %f\n",  *ra_used);
  fprintf(stderr, "*rainfall = %f\n",  rainfall);
  fprintf(stderr, "*ref_height = %f\n",  *ref_height);
  fprintf(stderr, "*roughness = %f\n",  *roughness);
  fprintf(stderr, "*wind = %f\n",  *wind);
 
  /* latent heat terms */
  fprintf(stderr, "Le = %f\n",   Le);

  /* snowpack terms */
  fprintf(stderr, "Advection = %f\n",  Advection);
  fprintf(stderr, "OldTSurf = %f\n",  OldTSurf);
  fprintf(stderr, "TPack = %f\n",  TPack);
  fprintf(stderr, "Tsnow_surf = %f\n",  Tsnow_surf);
  fprintf(stderr, "kappa_snow = %f\n",  kappa_snow);
  fprintf(stderr, "melt_energy = %f\n",  melt_energy);
  fprintf(stderr, "snow_coverage = %f\n",  snow_coverage);
  fprintf(stderr, "snow_density = %f\n",  snow_density);
  fprintf(stderr, "snow_swq = %f\n",  snow_swq);
  fprintf(stderr, "snow_water = %f\n",  snow_water);

  fprintf(stderr, "*deltaCC = %f\n",  *deltaCC);
  fprintf(stderr, "*refreeze_energy = %f\n",  *refreeze_energy);
  fprintf(stderr, "*VaporMassFlux = %f\n",  *VaporMassFlux);

  /* soil node terms */
  fprintf(stderr, "Nnodes = %i\n", Nnodes);

  /* spatial frost terms */
  fprintf(stderr, "*frost_fract = %f\n",  *frost_fract);

  /* control flags */
  fprintf(stderr, "INCLUDE_SNOW = %i\n", INCLUDE_SNOW);
  fprintf(stderr, "FS_ACTIVE = %i\n", FS_ACTIVE);
  fprintf(stderr, "NOFLUX = %i\n", NOFLUX);
  fprintf(stderr, "EXP_TRANS = %i\n", EXP_TRANS);
  fprintf(stderr, "SNOWING = %i\n", SNOWING);

  fprintf(stderr, "*FIRST_SOLN = %i\n", *FIRST_SOLN);

  /* returned energy balance terms */
  fprintf(stderr, "*NetLongBare = %f\n",  *NetLongBare); 
  fprintf(stderr, "*NetLongSnow = %f\n",  *NetLongSnow); 
  fprintf(stderr, "*T1 = %f\n",  *T1);
  fprintf(stderr, "*deltaH = %f\n",  *deltaH);
  fprintf(stderr, "*fusion = %f\n",  *fusion);
  fprintf(stderr, "*grnd_flux = %f\n",  *grnd_flux);
  fprintf(stderr, "*latent_heat = %f\n",  *latent_heat);
  fprintf(stderr, "*latent_heat_sub = %f\n",  *latent_heat_sub);
  fprintf(stderr, "*sensible_heat = %f\n",  *sensible_heat);
  fprintf(stderr, "*snow_flux = %f\n",  *snow_flux);
  fprintf(stderr, "*store_error = %f\n",  *store_error);

  write_layer(layer, iveg, options.Nlayer, frost_fract, depth);
  write_vegvar(&(veg_var[0]),iveg);

  if(!options.QUICK_FLUX) {
    fprintf(stderr,"Node\tT\tTnew\tZsum\tkappa\tCs\tmoist\tbubble\texpt\tmax_moist\tice\n");
    for(i=0;i<Nnodes;i++) 
      fprintf(stderr,"%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	      i, T_node[i], Tnew_node[i], Zsum_node[i], kappa_node[i], 
	      Cs_node[i], moist_node[i], bubble_node[i], expt_node[i], 
	      max_moist_node[i], ice_node[i]);
  }
  
  fprintf(stderr,"**********\n**********\nFinished writing calc_surf_energy_bal variables.\nTry increasing SURF_DT to get model to complete cell.\nThen check output for instabilities.\n**********\n**********\n");

  return (ERROR);
    
}

