#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double solve_snow(snow_data_struct    *snow,
		  layer_data_struct   *layer_wet,
		  layer_data_struct   *layer_dry,
		  veg_var_struct      *veg_var_wet,
		  veg_var_struct      *veg_var_dry,
		  int                  month,
		  int                  day_in_year,
		  energy_bal_struct   *energy,
		  soil_con_struct     *soil_con,
		  char                 overstory,
		  int                  dt,
		  int                  rec,
		  int                  veg_class,
		  int                  iveg,
		  int                  Nveg,
		  int                  band,
		  int                  hour,
		  int                  Nnodes,
		  double               shortwave,
		  double               longwave,
		  double               air_temp,
		  double               prec,
		  double               density,
		  double               vp,
		  double               vpd,
		  double               pressure,
		  double               mu,
		  double               roughness,
		  double               displacement,
		  double               ref_height,
		  double               surf_atten,
		  double               MAX_SNOW_TEMP,
		  double               MIN_RAIN_TEMP,
		  double               wind_h,
		  double               moist,
		  double               ice0,
		  double               dp,
		  double               bare_albedo,
		  double              *rainfall,
		  double              *out_prec,
		  double              *Le,
		  double              *Ls,
		  double              *aero_resist,
		  double              *aero_resist_used,
		  double              *tmp_wind,
		  double              *net_short,
		  double              *out_short,
		  double              *rad,
		  double              *Evap,
		  double              *tmp_snow_energy,
		  double              *snow_inflow,
		  double              *ppt,
		  double              *gauge_correction,
		  float               *root) {
/*********************************************************************
  solve_snow.c                Keith Cherkauer       July 2, 1998

  This routine was written to handle the various calls and data 
  handling needed to solve the various components of the new VIC
  snow code for both the full_energy and water_balance models.

  Returns snow, veg_var, and energy variables for each elevation
  band.  Variable ppt[] is defined for elevation bands with snow.

  07-13-98 modified to use elevation bands when solving the 
           snow model                                        KAC
  11-30-98 reworked the way the snow/rain fraction is computed
           and added to check to assure that very small amounts 
           of snow do not fall, causing snow sublimation to 
           be calculated.  (found by Greg)                   KAC
  02-29-00 removed ground heat flux computations, will now 
           make those outside of this routine, in the same function
	   that is used to compute the ground heat flux when
	   there is no snow cover.                           KAC
  06-04-03 Counter for number of days since last snowfall was 
           incremented twice in the MELTING update.  This has been
           fixed.                                            KAC
  06-04-03 Added check so that MELTING flag is only TRUE if melt
           occurs in the melt season - currently this is defined
           between March 1 and October 1.  Otherwise the MELTING
           flag can trigger rapid very early season melt     KAC
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    actually used in flux calculations.			TJB

*********************************************************************/

  extern option_struct   options;
  extern veg_lib_struct *veg_lib;

  char                ErrStr[MAXSTRING];
  char                FIRST_SOLN[1];
  double              rainonly;
  double              canopy_temp;
  double              tmp_energy_val = 0.;
  double              old_swq;
  double              tmp_Wdew[2];
  double              melt;
  double              snowfall[2];
  double              snow_coverage;
  double              grnd_temp;
  double              tmp_rain;
  double              surf_long;
  double              store_snowfall;
  int                 curr_snow;

  melt     = 0.;
  ppt[WET] = 0.; 
  ppt[DRY] = 0.; 

  /** Calculate Fraction of Precipitation that falls as Rain **/
  rainonly      = calc_rainonly(air_temp, prec, 
				MAX_SNOW_TEMP, MIN_RAIN_TEMP, mu);
  snowfall[WET] = gauge_correction[SNOW] * (prec - rainonly);
  rainfall[WET] = gauge_correction[RAIN] * rainonly;
  snowfall[DRY] = 0.;
  rainfall[DRY] = 0.;
  if(snowfall[WET] < 1e-5) snowfall[WET] = 0.;
  (*out_prec) = snowfall[WET] + rainfall[WET];
  store_snowfall = snowfall[WET];

  /** Compute latent heats **/
  (*Le) = (2.501 - 0.002361 * air_temp) * 1.0e6;
  (*Ls) = 0.;

  /** Error checks **/
  if((snow->swq > 0 || snowfall[WET] > 0.
      || (snow->snow_canopy>0. && overstory))) {
    if(mu!=1 && options.FULL_ENERGY) {
      sprintf(ErrStr,"Snow model cannot be used if mu (%f) is not equal to 1.\n\tsolve_snow.c: record = %i,\t vegetation type = %i",
	      mu, rec, iveg);
      vicerror(ErrStr);
    }
    else if(mu!=1) {
      fprintf(stderr,"WARNING: Snow is falling, but mu not equal to 1 (%f)\n", mu);
      fprintf(stderr,"\trec = %i, veg = %i, hour = %i\n",rec,iveg,hour);
    }
  }

  if((snow->swq > 0 || snowfall[WET] > 0.
      || (snow->snow_canopy>0. && overstory)) && mu==1) {
    
    /************************************************
      Snow is Present or Falling 
    ************************************************/

    /** Snow is present or falling **/
    snow->snow = TRUE;
	
    /** Initialize variables **/
    (*tmp_snow_energy) = 0.;
    if(!overstory) surf_atten = 1.;
    snow_coverage = snow->coverage;
      
    /** Compute Snow Pack Albedo **/
    if(snow->swq > 0 || snowfall[WET] > 0.) {
      snow->albedo   = snow_albedo( snowfall[WET], snow->swq, 
				    snow->coldcontent, dt, snow->last_snow,
				    snow->MELTING );
      energy->albedo = snow->albedo;
    }
    else {
      snow->albedo   = NEW_SNOW_ALB;
      energy->albedo = bare_albedo;
    }
    
    /** Age Snow Pack **/
    //if( snowfall[WET] > 0 )
      //snow->last_snow = 1;
    //else snow->last_snow++;
    
    /** Compute Radiation Balance over Snow **/ 
    (*out_short) = energy->albedo * shortwave;
    (*net_short) = (1.0 - energy->albedo) * shortwave;
    if(snow->swq > 0) {
      (*rad) = (*net_short) + longwave - STEFAN_B 
	* (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN) 
	* (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN);
    }
    else {
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
	(*rad) = (*net_short) + longwave - STEFAN_B 
	  * (energy->T[0]+KELVIN) * (energy->T[0]+KELVIN) 
	  * (energy->T[0]+KELVIN) * (energy->T[0]+KELVIN);
      }
      else {
	(*rad) = (*net_short) + longwave - STEFAN_B 
	  * (air_temp+KELVIN) * (air_temp+KELVIN) 
	  * (air_temp+KELVIN) * (air_temp+KELVIN);
      }
    }

    if(iveg!=Nveg) {
      
      /****************************************
	Check Vegetation for Intercepted Snow
      ****************************************/
      
      if(overstory) {
	if( snowfall[WET] > 0. || snow->snow_canopy > 0. ) {
	  
	  /** Compute Canopy Interception, if Snow and Canopy **/
	  if(snow->swq > 0) 
	    surf_long = STEFAN_B * (snow->surf_temp + KELVIN) 
	      * (snow->surf_temp + KELVIN) * (snow->surf_temp + KELVIN) 
	      * (snow->surf_temp + KELVIN);
	  else {
	    if(options.FULL_ENERGY || options.FROZEN_SOIL)
	      surf_long = STEFAN_B * (energy->T[0] + KELVIN) 
		* (energy->T[0] + KELVIN) * (energy->T[0] + KELVIN) 
		* (energy->T[0] + KELVIN);
	    else
	      surf_long = STEFAN_B * (air_temp + KELVIN) 
		* (air_temp + KELVIN) * (air_temp + KELVIN) 
		* (air_temp + KELVIN);
	  }
	  snow_intercept((double)dt, 1., veg_lib[veg_class].LAI[month-1], 
			 veg_lib[veg_class].Wdmax[month-1], 
			 aero_resist[1], aero_resist_used,
			 density, vp, (*Le), shortwave, 
			 longwave+surf_long, pressure, air_temp, vpd, 
			 tmp_wind[1], rainfall, snowfall, &veg_var_wet->Wdew, 
			 &snow->snow_canopy, &snow->tmp_int_storage, 
			 &snow->canopy_vapor_flux, &canopy_temp, 
			 &tmp_energy_val, month, rec, hour);

	  /* Store throughfall from canopy */
	  veg_var_wet->throughfall = rainfall[0] + snowfall[0];

	  /* Estimate longwave radiation from canopy */
	  energy->longwave = STEFAN_B * (canopy_temp+KELVIN) 
	    * (canopy_temp+KELVIN) * (canopy_temp+KELVIN) 
	    * (canopy_temp+KELVIN);

	}
	else {
	  /** Compute Canopy Evaporation, if Canopy and No Snow **/
	  tmp_Wdew[WET] = veg_var_wet->Wdew;
	  if(options.DIST_PRCP) tmp_Wdew[DRY] = 0.;
          *aero_resist_used = aero_resist[0];
	  Evap[0] = canopy_evap(layer_wet, layer_dry, veg_var_wet, 
				veg_var_dry, FALSE, veg_class, month, mu, 
				tmp_Wdew, (double)dt, rad[0], vpd, 
				net_short[0], air_temp, *aero_resist_used,
				displacement, roughness, ref_height, 
				(double)soil_con->elevation, rainfall, 
				soil_con->depth, soil_con->Wcr, 
				soil_con->Wpwp, root);

	  /* Store throughfall from canopy */
	  rainfall[WET] = veg_var_wet->throughfall;
	  if(options.DIST_PRCP) 
	    rainfall[DRY] = veg_var_dry->throughfall;

	  /* Estimate longwave radiation from canopy */
	  energy->longwave = (STEFAN_B * (air_temp+KELVIN) * (air_temp+KELVIN) 
			      * (air_temp+KELVIN) * (air_temp+KELVIN));

	}

      }

      else if(snowfall[0] > 0. && veg_var_wet->Wdew > 0.) {

	/** If No Overstory, Empty Vegetation of Stored Water **/

	tmp_rain = calc_rainonly(air_temp, veg_var_wet->Wdew,
				 MAX_SNOW_TEMP, MIN_RAIN_TEMP, mu);
	rainfall[WET]            += tmp_rain;
	snowfall[WET]            += veg_var_wet->Wdew - tmp_rain;
	veg_var_wet->throughfall  = rainfall[WET] + snowfall[WET];
	veg_var_wet->Wdew         = 0.;
	energy->longwave          = longwave;

      }
      else {

	/** Precipitation "Passes Through" Vegetation which 
	    is Under Snow (used only for accounting purposes)**/

	veg_var_wet->throughfall = rainfall[WET] + snowfall[WET];
	veg_var_dry->throughfall = rainfall[DRY] + snowfall[DRY];
	energy->longwave         = longwave;

      }

    }
    else energy->longwave = longwave;
    
    old_swq = snow->swq;
    
    if(snow->swq>0.0 || snowfall[0] > 0) {
      
      /******************************
	Snow Pack Present on Ground
      ******************************/

      // store snowfall reaching the ground for determining the albedo
      store_snowfall            = snowfall[WET];

      /** Age Snowpack **/
      if( snowfall[WET] > 0 ) curr_snow = 1; // new snow - reset pack age
      else curr_snow = snow->last_snow + 1; // age pack by one time step

      if(overstory && snowfall[0] > 0.) {

	/** recompute surface properties if overstory drops snow **/

	snow->albedo    = NEW_SNOW_ALB;
	energy->albedo  = snow->albedo;
	snow->last_snow = 1;
	(*out_short)    = energy->albedo * shortwave;
	(*net_short)    = (1.0 - energy->albedo) * shortwave;
	(*rad)          = (*net_short) + longwave - STEFAN_B 
	  * (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN) 
	  * (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN);
      }
    
      (*snow_inflow) += rainfall[WET] + snowfall[WET];
      if(options.FULL_ENERGY || options.FROZEN_SOIL) grnd_temp = energy->T[0];
      else grnd_temp = air_temp;

      /** Call snow pack accumulation and ablation algorithm **/

      snow_melt(soil_con, rec, iveg, wind_h+soil_con->snow_rough, 
		aero_resist[2], aero_resist_used, Le[0], snow, (double)dt, 0.00, 
		soil_con->snow_rough, surf_atten, rainfall[WET], 
		snowfall[WET], tmp_wind[2], grnd_temp, air_temp, net_short[0], 
		energy->longwave, density, pressure, vpd, vp, &melt, 
		&energy->advection, &energy->deltaCC, &energy->grnd_flux, 
		&energy->latent, &energy->sensible, &energy->error, 
		&energy->refreeze_energy);

      ppt[WET] += melt;
      energy->albedo   = snow->albedo;
      tmp_snow_energy[0] = energy->advection - energy->deltaCC 
	+ energy->refreeze_energy;
      (*Ls) = (677. - 0.07 * snow->surf_temp) * 4.1868 * 1000;
      
      /** Compute Snow Parameters **/
      if(snow->swq > 0.) {

	/** Calculate Snow Density **/
	if ( snow->surf_temp <= 0 )
	  // snowpack present, compress and age density
	  snow->density = snow_density(day_in_year, snowfall[WET], air_temp, 
				       old_swq, snow->depth, snow->coldcontent, 
				       (double)dt, snow->surf_temp);
	else 
	  // no snowpack present, start with new snow density
	  if ( curr_snow == 1 ) 
	    snow->density = new_snow_density(air_temp);

	/** Calculate Snow Depth (H.B.H. 7.2.1) **/
	snow->depth = 1000. * snow->swq / snow->density; 
	
	/** Record if snowpack is melting this time step **/
	if ( snow->coldcontent >= 0 && day_in_year > 60 // ~ March 1
	     && day_in_year < 273 // ~ October 1
	     ) snow->MELTING = TRUE;
	else if ( snow->MELTING && snowfall[WET] > TraceSnow ) 
	  snow->MELTING = FALSE;
	
	/** Check for Thin Snowpack which only Partially Covers Grid Cell **/
	if(snow->swq < MAX_FULL_COVERAGE_SWQ) {
	  snow->coverage = 1.; 
/* 	  snow->coverage = 1. / MAX_FULL_COVERAGE_SWQ   */
/* 	    * snow->swq;  */
/* 	  if(snow_coverage <= snow->coverage) { */
/* 	    energy->albedo = bare_albedo; */
/* 	  } */
/* 	  else { */
/* 	    energy->albedo = (snow_coverage - snow->coverage)  */
/* 	      / (1. - snow->coverage) * snow->albedo; */
/* 	    energy->albedo += (1. - snow_coverage)  */
/* 	      / (1. - snow->coverage) * bare_albedo; */
/* 	  } */
	}
	else {
	  snow->coverage = 1.;
	}

	/** Estimate net longwave at the snow surface */
	energy->longwave -= (STEFAN_B * (snow->surf_temp+KELVIN) 
			     * (snow->surf_temp+KELVIN) 
			     * (snow->surf_temp+KELVIN) 
			     * (snow->surf_temp+KELVIN));

      }
      else {

	/** Reset Snow Pack Variables after Complete Melt **/

	snow->density    = 0.;
	snow->depth      = 0.;
	snow->surf_water = 0;
	snow->pack_water = 0;
	snow->surf_temp  = 0;
	snow->pack_temp  = 0;
	snow->coverage   = 0;
	snow->MELTING    = FALSE;

 	energy->albedo = (snow_coverage - snow->coverage)  
 	  / (1. - snow->coverage) * snow->albedo; 
 	energy->albedo += (1. - snow_coverage)  
 	  / (1. - snow->coverage) * bare_albedo; 
	
      }
    }

    else {

      /** Ground Snow not Present, and Falling Snow Does not Reach Ground **/

      ppt[WET] += rainfall[WET];
      tmp_snow_energy[0] = 0.;
      curr_snow               = 0;
      snow->MELTING           = FALSE;

    }
    
    if ( store_snowfall > TraceSnow || store_snowfall == 0 ) {
      // reset snow albedo ago if new snow is sufficiently deep
      snow->last_snow = curr_snow;
    }
    else {
      snow->last_snow++;
    }

  }
  else {
    
    /*****************************
      No Snow Present or Falling
    *****************************/

    /** Initialize variables **/
    snow->snow         = FALSE;
    tmp_snow_energy[0] = 0;
    energy->albedo     = bare_albedo;
    energy->longwave   = longwave;

    /** Compute Radiation Balance for Bare Surface **/ 
    out_short[0] = energy->albedo * shortwave;
    net_short[0] = (1.0 - energy->albedo) * shortwave;
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      rad[0]       = net_short[0] + longwave 
	- STEFAN_B * (energy->T[0]+KELVIN) * (energy->T[0]+KELVIN) 
	* (energy->T[0]+KELVIN) * (energy->T[0]+KELVIN);
    }
    else {
      rad[0]       = net_short[0] + longwave 
	- STEFAN_B * (air_temp+KELVIN) * (air_temp+KELVIN) 
	* (air_temp+KELVIN) * (air_temp+KELVIN);
    }
    snow->MELTING        = FALSE;
  }

  energy->shortwave = (*net_short);

  return(melt);

}

