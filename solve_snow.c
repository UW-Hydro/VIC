#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double solve_snow(snow_data_struct    *snow,
		  layer_data_struct   *layer_wet,
		  layer_data_struct   *layer_dry,
		  veg_var_struct      *veg_var_wet,
		  veg_var_struct      *veg_var_dry,
		  dmy_struct          *dmy,
		  energy_bal_struct   *energy,
		  soil_con_struct      soil_con,
		  char                 overstory,
		  char                *SNOW,
		  char                *BARE,
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
		  double              *Le,
		  double              *Ls,
		  double              *aero_resist,
		  double              *tmp_wind,
		  double              *net_short,
		  double              *out_short,
		  double              *rad,
		  double              *Evap,
		  double              *new_T1,
		  double              *Tsurf, /** half time step (e bal) **/
		  double              *Tgrnd, /** half time step (e bal) **/
		  double              *Tend_surf, /** end of time step **/
		  double              *Tend_grnd, /** end of time step **/
		  double              *tmp_snow_energy,
		  double              *snow_inflow,
		  double              *ppt,
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

*********************************************************************/

  extern option_struct   options;
  extern debug_struct    debug;
  extern veg_lib_struct *veg_lib;

  char                ErrStr[MAXSTRING];
  char                FIRST_SOLN[1];
  double              rainonly;
  double              canopy_temp;
  double              tmp_energy_val;
  double              tmp_swq;
  double              tmp_Wdew[2];
  double              tmp_air_temp;
  double              melt;
  double              Melt;
  double              rainfall[2];
  double              snowfall[2];
  double              snow_coverage;
  double              grnd_temp;
  double              tmp_rain;
  double              surf_long;
  energy_bal_struct   tmp_energy;

  melt           = 0.;
  Melt           = 0.;
  ppt[WET]       = 0.; 
  ppt[DRY]       = 0.; 
  tmp_energy_val = 0.;

  /** Check if Snow Model is Needed, if so Then Solve **/
  SNOW[0]   = FALSE;

  /** Reset Temporary Storage Structures **/
  tmp_energy   = energy[0];
  tmp_air_temp = air_temp + soil_con.Tfactor[band];

  /** Calculate Fraction of Precipitation that falls as Rain **/
  prec          *= soil_con.Pfactor[band];
  rainonly       = calc_rainonly(tmp_air_temp,prec,
				 MAX_SNOW_TEMP,MIN_RAIN_TEMP,mu);
  snowfall[WET]  = prec - rainonly;
  rainfall[WET]  = rainonly;
  snowfall[DRY]  = 0.;
  rainfall[DRY]  = 0.;
  if(snowfall[WET] < 1e-5) snowfall[WET] = 0.;
 
  Le[0] = (2.501 - 0.002361 * tmp_air_temp) * 1.0e6;
  Ls[0] = 0.;
  
  if(options.SNOW_MODEL && (snow->snow || snowfall[0] > 0.
			    || (snow->snow_canopy>0. && overstory))) {
    if(mu!=1 && options.FULL_ENERGY) {
      sprintf(ErrStr,"Snow model cannot be used if mu (%lf) is not equal to 1.\n\tsolve_snow.c: record = %i,\t vegetation type = %i",
	      mu,rec,iveg);
      vicerror(ErrStr);
    }
    else if(mu!=1) {
      fprintf(stderr,"WARNING: Snow is falling, but mu not equal to 1 (%lf)\n",
	      mu);
      fprintf(stderr,"\trec = %i, veg = %i, hour = %i\n",rec,iveg,hour);
    }
  }

  if(options.SNOW_MODEL && (snow->snow || snowfall[0] > 0.
			    || (snow->snow_canopy>0. && overstory)) && mu==1) {
    
    /************************************************
      Snow is Present or Falling 
    ************************************************/

    SNOW[0]            = TRUE;
    tmp_snow_energy[0] = 0.;
    if(!overstory) surf_atten = 1.;
    snow_coverage = snow->coverage;
      
    /** Compute Snow Pack Albedo **/
    if(snow->snow || snowfall[0] > 0.) {
      snow->albedo   = snow_albedo(snowfall[0],snow->swq,
				   snow->coldcontent,dt,
				   snow->last_snow);
      energy->albedo = snow->albedo /* * snow->coverage  */
/* 	+ bare_albedo * (1. - snow->coverage) */;
    }
    else {
      snow->albedo   = NEW_SNOW_ALB;
      energy->albedo = bare_albedo;
    }
    
    /** Age Snow Pack **/
    if(snowfall[0]>0)
      snow->last_snow=1;
    else snow->last_snow++;
    
    /** Compute Radiation Balance over Snow **/ 
    out_short[0] = energy->albedo * shortwave;
    net_short[0] = (1.0 - energy->albedo) * shortwave;
    if(snow->snow) {
      rad[0]       = net_short[0] + longwave 
	- STEFAN_B * (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN) 
	* (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN);
    }
    else {
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
    }

    
    if(iveg!=Nveg) {
      
      /****************************************
	Check Vegetation for Intercepted Snow
      ****************************************/
      
      if(overstory) {
	if(snowfall[0]>0. || snow->snow_canopy>0.) {
	  
	  /** Compute Canopy Interception, if Snow and Canopy **/
	  if(snow->snow) 
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
	      /* surf_long = longwave; */
	  }
	  snow_intercept((double)dt,1.,
			 veg_lib[veg_class].LAI[dmy[rec].month-1],
			 veg_lib[veg_class].Wdmax[dmy[rec].month-1],
			 aero_resist[1],density,vp,Le[0],shortwave,
			 longwave+surf_long,pressure,tmp_air_temp,vpd,
			 tmp_wind[1],rainfall,snowfall,
			 &veg_var_wet->Wdew,
			 &snow->snow_canopy,
			 &snow->tmp_int_storage,
			 &snow->canopy_vapor_flux,&canopy_temp,
			 &tmp_energy_val,dmy[rec].month,rec,hour);
	  energy->longwave = STEFAN_B * (canopy_temp+KELVIN) 
	    * (canopy_temp+KELVIN) * (canopy_temp+KELVIN) 
	    * (canopy_temp+KELVIN);
	  veg_var_wet->throughfall = rainfall[0] + snowfall[0];
	}
	else {
	  /** Compute Canopy Evaporation, if Canopy and No Snow **/
	  tmp_Wdew[WET] = veg_var_wet->Wdew;
	  if(options.DIST_PRCP) tmp_Wdew[DRY] = 0.;
	  Evap[0] = canopy_evap(layer_wet,layer_dry,veg_var_wet,
				veg_var_dry,
				FALSE,veg_class,dmy[rec].month,mu,tmp_Wdew,
				tmp_air_temp,(double)dt,rad[0],vpd,
				net_short[0],
				tmp_air_temp,aero_resist[0],displacement,
				roughness,ref_height,
				(double)soil_con.elevation,
				rainfall,soil_con.depth,soil_con.Wcr,
				soil_con.Wpwp,root);
	  rainfall[WET] = veg_var_wet->throughfall;
	  energy->longwave = STEFAN_B * (tmp_air_temp+KELVIN) 
	    * (tmp_air_temp+KELVIN) * (tmp_air_temp+KELVIN) 
	    * (tmp_air_temp+KELVIN);
	  if(options.DIST_PRCP) 
	    rainfall[DRY] = veg_var_dry->throughfall;
	}
      }
      else if(snowfall[0] > 0. && veg_var_wet->Wdew > 0.) {
	/** If No Overstory, Empty Vegetation of Stored Water **/
	tmp_rain                       
	  = calc_rainonly(tmp_air_temp, veg_var_wet->Wdew,
			  MAX_SNOW_TEMP,MIN_RAIN_TEMP,mu);
	rainfall[WET]                 += tmp_rain;
	snowfall[WET]                 += veg_var_wet->Wdew - tmp_rain;
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
    
    tmp_swq = snow->swq;
    
    if(snow->swq>0.0 || snowfall[0] > 0) {
      
      /******************************
	Snow Pack Present on Ground
      ******************************/

      if(overstory && snowfall[0] > 0.) {
	/** recompute surface properties if overstory drops snow **/
	snow->albedo    = NEW_SNOW_ALB;
	energy->albedo  = snow->albedo;
	snow->last_snow = 1;
	out_short[0]    = energy->albedo * shortwave;
	net_short[0]    = (1.0 - energy->albedo) * shortwave;
	rad[0]          = net_short[0] + longwave 
	  - STEFAN_B * (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN) 
	  * (snow->surf_temp+KELVIN) * (snow->surf_temp+KELVIN);
      }
    
      snow_inflow[0] += rainfall[WET] + snowfall[WET];
      if(options.FULL_ENERGY || options.FROZEN_SOIL) grnd_temp = energy->T[0];
      else grnd_temp = air_temp;
      snow_melt(soil_con,rec,iveg,wind_h+soil_con.snow_rough,aero_resist[2],
		Le[0],snow,(double)dt,0.00,soil_con.snow_rough,
		surf_atten,rainfall[WET],snowfall[WET],tmp_wind[2],
		grnd_temp,tmp_air_temp,net_short[0],
		energy->longwave,density,pressure,vpd,vp,&melt,
		&energy->advection,&energy->deltaCC,
		&energy->grnd_flux,&energy->latent,
		&energy->sensible,&energy->error,
		&energy->refreeze_energy);
      ppt[WET]             += melt;
      energy->albedo   = snow->albedo;
      tmp_snow_energy[0] = energy->advection 
	- energy->deltaCC + energy->refreeze_energy;
      Ls[0] = (677. - 0.07 * snow->surf_temp) * 4.1868 * 1000;
      Tend_surf[0] = snow->surf_temp;
      
      /** Compute Snow Parameters **/
      if(snow->swq > 0.) {
	snow->snow = TRUE;
	
	/** Calculate Snow Density **/
	snow->density = snow_density(dmy[rec].day_in_year,
				     snowfall[0],tmp_air_temp,
				     tmp_swq,snow->depth,
				     snow->coldcontent,(double)dt,
				     snow->surf_temp);
	
	/** Calculate Snow Depth (H.B.H. 7.2.1) **/
	snow->depth = 1000. * snow->swq 
	  / snow->density; 
	
	/** Check for Thin Snowpack which only Partially Covers Grid Cell **/
	if(snow->swq < MAX_FULL_COVERAGE_SWQ) {
	  BARE[0]          = FALSE; 
	  snow->coverage = 1. / MAX_FULL_COVERAGE_SWQ  
	    * snow->swq; 
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
	  BARE[0]          = FALSE;
	  snow->coverage = 1.;
	}
	
	if(options.FULL_ENERGY  || options.FROZEN_SOIL) {
	  /** Compute Surface Energy Balance with Snow Present **/
	  if(options.CALC_SNOW_FLUX) {
	    /** Compute Thermal Flux Through the Snow Pack **/
	    FIRST_SOLN[0] = TRUE;
	    Tend_grnd[0] 
	      = calc_snow_ground_flux(dt,Nnodes,rec,iveg,mu,dp,moist,ice0,
				      snow->surf_temp,new_T1,&energy->grnd_flux,
				      &energy->deltaH,&energy->snow_flux,
				      energy,snow,layer_wet,layer_dry,
				      soil_con,FIRST_SOLN);
	  }
	  else if(options.FULL_ENERGY || options.FROZEN_SOIL) {
	    /** Ignore Thermal Flux Through the Snow Pack **/
	    Tend_grnd[0] = snow->surf_temp;
	    new_T1[0] = energy->T[1];
	  }
	}
      }
      else {
	/** Reset Snow Pack Variables after Complete Melt **/
	snow->snow       = FALSE;
	snow->density    = 0.;
	snow->depth      = 0.;
	snow->surf_water = 0;
	snow->pack_water = 0;
	snow->surf_temp  = 0;
	snow->pack_temp  = 0;
	snow->coverage   = 0;
	BARE[0]          = TRUE;

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
      BARE[0] = TRUE;
    }
    
    /** Store Energy and Water Components for the Current Elevation Band **/
    Melt                  += melt;
    
  }
  else {
    
    /*****************************
      No Snow Present or Falling
    *****************************/

    /** ppt will be set outside this routine **/
    
    BARE[0]            = TRUE;
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

  }

  return(Melt);

}




