#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id: surface_fluxes.c,v 4.2.2.6 2006/09/27 16:45:49 vicadmin Exp $";

void surface_fluxes(char                 overstory,
		    int                  rec,
		    int                  band,
		    int                  veg_class,
		    int                  iveg,
		    int                  Nveg,
		    int                  Ndist,
		    int                  Nbands,
		    int                  Nlayers,
		    int                  dp,
		    double               mu,
		    double               ice0,
		    double               moist,
		    double               surf_atten,
		    double               height,
		    double               displacement,
		    double               roughness,
		    double               ref_height,
		    double               bare_albedo,
		    double              *aero_resist,
		    double              *aero_resist_used,
		    double              *baseflow_wet,
		    double              *baseflow_dry,
		    double              *runoff_wet,
		    double              *runoff_dry,
		    double              *out_prec,
		    double              *out_rain,
		    double              *out_snow,
		    double              *wind,
		    double              *Le,
		    double              *Ls,
		    double              *Melt,
		    double              *inflow_wet,
		    double              *inflow_dry,
		    double              *snow_inflow,
		    double              *gauge_correction,
		    double              ra_lake,
		    float               *root,
		    double              local_water,
		    atmos_data_struct   *atmos,
		    soil_con_struct     *soil_con,
		    dmy_struct          *dmy,
		    global_param_struct *gp,
		    energy_bal_struct   *energy,
		    snow_data_struct    *snow,
		    layer_data_struct   *layer_wet,
		    layer_data_struct   *layer_dry,
		    veg_var_struct      *veg_var_wet,
		    veg_var_struct      *veg_var_dry,
		    int flag_irr)
/**********************************************************************
	surface_fluxes	Keith Cherkauer		February 29, 2000

  Formerally a part of full_energy.c this routine computes all surface
  fluxes, and solves the snow accumulation and ablation algorithm.
  Solutions are for the current snow band and vegetation type (these
  are defined in full_energy before the routine is called).

  modifications:
  02-07-03 fixed indexing problem for sub-daily snow model within
           daily water balance VIC: hour (now hidx) is incremented
           by 1 rather than the sub-daily time step, so the atmospheric
           forcing data is now properly indexed.                KAC
  04-23-03 Indexing fix sent SNOW_STEP to calc_surf_energy_bal rather
           than the model time step, meaning that without snow the 
           evaporation was computed for SNOW_STEP hours rather than a
           full day.  This was fixed by introducing step_inc to 
           index the arrays, while step_dt keeps track of the correct
           time step.                                            KAC
  07-May-04 Fixed initialization of canopyevap to initialize for every
	    value of dist, rather than just dist 0.		TJB
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    actually used in flux calculations.			TJB
  2006-Sep-11 Implemented flexible output configuration; moved tracking
              of rain and snow for output to this function.  TJB
  2006-Sep-26 Moved tracking of out_rain and out_snow to solve_snow.c.  TJB

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  double                 total_store_moist[3];
  double                 step_store_moist[3];

  int                    dist;
  int                    lidx;
  int                    endhidx;
  int                    hidx;
  int                    step_inc;
  int                    step_dt;
  int                    N_steps;
  double                 Tsurf;
  double                 air_temp;
  double                 Evap;
  double                 T0;
  double                 step_rad;
  double                 step_net_short;
  double                 ppt[2];
  double                 rainfall[2];
  double                 step_ppt[2];
  double                 step_snow_energy;
  double                 step_out_prec;
  double                 step_out_rain;
  double                 step_out_snow;
  double                 step_out_short;
  double                 step_Evap;
  double                 step_melt;
  double                 step_prec[2];
  double                 store_throughfall[2];
  double                 store_melt;
  double                 store_vapor_flux;
  double                 store_canopy_vapor_flux;
  double                 store_canopyevap[2];
  double                 store_layerevap[2][MAX_LAYERS];
  double                 store_potevap_lake[2]; //ingjerd dec 2008
  double                 store_potevap[2]; //ingjerd oct 2009
  double                 store_ppt[2];
  double                 store_shortwave;
  double                 store_longwave;
  double                 store_sensible;
  double                 store_latent;
  double                 store_grnd_flux;
  double                 store_deltaH;
  double                 store_advection; 
  double                 store_deltaCC; 
  double                 store_snow_flux; 
  double                 store_refreeze_energy; 
  double                 store_albedo;
  layer_data_struct      step_layer[2][MAX_LAYERS];
  energy_bal_struct      snow_energy;
  energy_bal_struct      bare_energy;
  snow_data_struct       step_snow;
  veg_var_struct         snow_veg_var[2];
  veg_var_struct         bare_veg_var[2];

  /***********************************************************************
    Set temporary variables - preserves original values until iterations
    are completed
  ***********************************************************************/

  energy->advection       = 0; 
  energy->deltaCC         = 0; 
  energy->snow_flux       = 0; 
  energy->refreeze_energy = 0; 
  snow_energy             = (*energy);
  bare_energy             = (*energy);
  snow_veg_var[WET]       = (*veg_var_wet);
  snow_veg_var[DRY]       = (*veg_var_dry);
  bare_veg_var[WET]       = (*veg_var_wet);
  bare_veg_var[DRY]       = (*veg_var_dry);
  step_snow               = (*snow);
  for(lidx=0;lidx<Nlayers;lidx++) {
    step_layer[WET][lidx] = layer_wet[lidx];
    step_layer[DRY][lidx] = layer_dry[lidx];
  }

  /********************************
    Set-up sub-time step controls
    (May eventually want to set this up so that it is also true 
    if frozen soils are present)
  ********************************/

  if(snow->swq > 0 || snow->snow_canopy > 0 || atmos->snowflag[NR]) {
    hidx      = 0;
    endhidx   = hidx + NF;
    step_inc  = 1;
    step_dt   = options.SNOW_STEP;
  }
  else {
    hidx      = NR;
    endhidx   = hidx + gp->dt;
    step_inc  = step_dt = gp->dt;
  }

  /*******************************************
    Initialize sub-model time step variables
  *******************************************/

  for(dist = 0; dist < Ndist; dist++) {
    store_throughfall[dist] = 0.;
    store_canopyevap[dist]  = 0.;
    store_potevap_lake[dist]  = 0.; /* ingjerd dec 2008 */
    store_potevap[dist]  = 0.; /* ingjerd dec 2008 */
    for(lidx = 0; lidx < options.Nlayer; lidx++) {
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
  store_shortwave         = 0;
  store_longwave          = 0;
  store_sensible          = 0;
  store_latent            = 0;
  store_grnd_flux         = 0;
  store_deltaH            = 0;
  store_advection         = 0; 
  store_deltaCC           = 0; 
  store_snow_flux         = 0; 
  store_refreeze_energy   = 0; 
  store_albedo            = 0;
  N_steps                 = 0;
      
  /*************************
    Compute surface fluxes 
  *************************/

  do {

    /*********************************************************
      Solve for snow interception, accumulation and ablation
    *********************************************************/

    air_temp = atmos->air_temp[hidx] + soil_con->Tfactor[band];
    step_prec[WET] = atmos->prec[hidx] / mu * soil_con->Pfactor[band];

    for(dist = 0; dist < Ndist; dist ++) {
      snow_veg_var[dist].canopyevap = 0;
      bare_veg_var[dist].canopyevap = 0;
      bare_veg_var[dist].potevap_lake = 0.; //ingjerd dec 2008
      bare_veg_var[dist].potevap = 0.; //ingjerd oct 2009
      for(lidx = 0; lidx < Nlayers; lidx ++) 
	step_layer[dist][lidx].evap = 0;
    }
    step_snow.canopy_vapor_flux = 0;
    step_snow.vapor_flux = 0;

    if ( options.GRND_FLUX ) T0 = energy->T[0];
    else T0 = air_temp;
    //if(rec<10) printf("surface_fluxes rec%d vp %f vpd %f flag_irr %d\n",rec,atmos->vp[hidx],atmos->vpd[hidx],flag_irr);    
    step_melt = solve_snow( &(step_snow), step_layer[WET], step_layer[DRY],
			    &(snow_veg_var[WET]), &(snow_veg_var[DRY]), 
			    dmy[rec].month, dmy[rec].day_in_year,
			    &(snow_energy), soil_con, overstory, 
			    step_dt, 
			    rec, veg_class, iveg, Nveg, band, dmy[rec].hour, 
			    options.Nnode, atmos->shortwave[hidx], 
			    atmos->longwave[hidx], air_temp, 
			    step_prec[WET], atmos->density[hidx], 
			    atmos->vp[hidx], atmos->vpd[hidx], 
			    atmos->pressure[hidx], mu, roughness, 
			    displacement, ref_height, surf_atten, 
			    gp->MAX_SNOW_TEMP, gp->MIN_RAIN_TEMP, 
			    gp->wind_h, moist, ice0, dp, bare_albedo, 
			    rainfall, &step_out_prec, &step_out_rain, &step_out_snow,
			    Le, Ls, aero_resist, aero_resist_used,
			    wind, ra_lake, &step_net_short, &step_out_short, &step_rad, 
			    &step_Evap, &step_snow_energy, snow_inflow, 
			    step_ppt, gauge_correction, root, flag_irr);
      
    /**********************************************************
      Solve Energy Balance Components for Ground Free of Snow
    **********************************************************/
    //printf("surface_fluxes hidx %d longwave %f\n",hidx,atmos->longwave[hidx]);

    Tsurf = calc_surf_energy_bal(rec, iveg, options.Nlayer, Nveg, 
				 step_dt, 
				 options.Nnode, veg_class, band, hidx, 
				 ice0, moist, dp, surf_atten, T0, 
				 atmos->shortwave[hidx], 
				 snow_energy.longwave, atmos->longwave[hidx],
				 air_temp, (*Le), 
				 (*Ls), mu, displacement, roughness, 
				 ref_height, step_snow_energy, 
				 step_snow.vapor_flux,  bare_albedo, 
				 aero_resist, aero_resist_used, wind, 
				 rainfall, step_ppt, root, 
				 atmos, &bare_veg_var[WET], 
				 &bare_veg_var[DRY], &bare_energy, 
				 &(step_snow), step_layer[WET], 
				 step_layer[DRY], soil_con, &(dmy[rec]),flag_irr,
				 ra_lake);  
	      	      
    //printf("surface_fluxes hidx %d longwave %f bare_energylongwave %f\n",hidx,atmos->longwave[hidx],bare_energy.longwave);
    /**************************************
      Store sub-model time step variables 
    **************************************/

    for(dist = 0; dist < Ndist; dist++) {

      if(iveg < Nveg) {
	if(step_snow.snow) {
	  store_throughfall[dist] += snow_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += snow_veg_var[dist].canopyevap;
	  bare_veg_var[dist].Wdew  = snow_veg_var[dist].Wdew;
	  store_potevap_lake[dist]  += 0; //ingjerd jan 2010. assume lake evap = 0 if snow
	  store_potevap[dist]  += 0; //ingjerd jan 2010. assume potevap=0 if snow
	}
	else {
	  store_throughfall[dist] += bare_veg_var[dist].throughfall;
	  store_canopyevap[dist]  += bare_veg_var[dist].canopyevap;
	  snow_veg_var[dist].Wdew  = bare_veg_var[dist].Wdew;
	  store_potevap_lake[dist] += bare_veg_var[dist].potevap_lake; //ingjerd dec 2008
	  store_potevap[dist]     += bare_veg_var[dist].potevap; //ingjerd dec 2008
	}
      }

      //     printf("surface_fluxes veg=%d %f %f %f %f\n",iveg,store_potevap[dist],bare_veg_var[dist].potevap,store_canopyevap[dist],bare_veg_var[dist].canopyevap);

      for(lidx = 0; lidx < options.Nlayer; lidx++)
	store_layerevap[dist][lidx] += step_layer[dist][lidx].evap;

      store_ppt[dist] += step_ppt[dist];
    }


    if(iveg < Nveg) 
      store_canopy_vapor_flux += step_snow.canopy_vapor_flux;
    store_vapor_flux += step_snow.vapor_flux;
      
    store_melt  += step_melt;
    out_prec[0] += step_out_prec * mu;
    out_rain[0] += step_out_rain * mu;
    out_snow[0] += step_out_snow * mu;
   
    if(overstory)
      store_shortwave += step_snow.coverage * snow_energy.shortwave 
	* surf_atten + (1. - step_snow.coverage) * bare_energy.shortwave; 
    else
      store_shortwave += step_snow.coverage * snow_energy.shortwave 
	+ (1. - step_snow.coverage) * bare_energy.shortwave; 
    store_longwave  += step_snow.coverage * snow_energy.longwave 
      + (1. - step_snow.coverage) * bare_energy.longwave; 
    store_latent    += step_snow.coverage * snow_energy.latent 
      + (1. - step_snow.coverage) * bare_energy.latent; 
    store_sensible  += step_snow.coverage * snow_energy.sensible 
      + (1. - step_snow.coverage) * bare_energy.sensible; 
    store_grnd_flux += bare_energy.grnd_flux; 
    store_deltaH    += bare_energy.deltaH; 
    store_albedo    += step_snow.coverage * snow_energy.albedo
      + (1. - step_snow.coverage) * bare_energy.albedo;
    if (step_snow.snow) {
      store_advection       += snow_energy.advection; 
      store_deltaCC         += snow_energy.deltaCC; 
      store_snow_flux       += bare_energy.snow_flux; 
      store_refreeze_energy += snow_energy.refreeze_energy; 
    }

    /* increment time step */
    N_steps ++;
    hidx += step_inc;

  } while (hidx < endhidx);

  /************************************************
    Store snow variables for sub-model time steps 
  ************************************************/

  (*snow) = step_snow;
  snow->vapor_flux = store_vapor_flux;
  snow->canopy_vapor_flux = store_canopy_vapor_flux;
  (*Melt) = store_melt;
  for(dist = 0; dist < 2; dist++) ppt[dist] = store_ppt[dist];

  /******************************************************
    Store energy flux averages for sub-model time steps 
  ******************************************************/

  (*energy) = bare_energy;
  if(overstory)
    energy->shortwave = store_shortwave / (double)N_steps; 
  else
    energy->shortwave = store_shortwave / (double)N_steps; 
  energy->longwave    = store_longwave / (double)N_steps; 
  energy->latent      = store_latent / (double)N_steps; 
  energy->sensible    = store_sensible / (double)N_steps; 
  energy->grnd_flux   = store_grnd_flux / (double)N_steps; 
  energy->deltaH      = store_deltaH / (double)N_steps; 
  energy->albedo      = store_albedo / (double)N_steps;
  if (snow->snow) {
    energy->advection       = store_advection / (double)N_steps; 
    energy->deltaCC         = store_deltaCC / (double)N_steps; 
    energy->snow_flux       = store_snow_flux / (double)N_steps; 
    energy->refreeze_energy = store_refreeze_energy / (double)N_steps; 
  }
	  
  /**********************************************************
    Store vegetation variable sums for sub-model time steps 
  **********************************************************/

  if(iveg < Nveg) {
    veg_var_wet->throughfall = store_throughfall[WET];
    veg_var_dry->throughfall = store_throughfall[DRY];
    veg_var_wet->canopyevap  = store_canopyevap[WET];
    veg_var_dry->canopyevap  = store_canopyevap[DRY];
    veg_var_wet->potevap_lake   = store_potevap_lake[WET]; //ingjerd dec 2008
    veg_var_dry->potevap_lake   = store_potevap_lake[DRY]; //ingjerd dec 2008
    veg_var_wet->potevap   = store_potevap[WET]; //ingjerd dec 2008
    veg_var_dry->potevap   = store_potevap[DRY]; //ingjerd dec 2008

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
	 baseflow_wet, baseflow_dry, ppt, mu, gp->dt, options.Nnode, 
	 band, rec, iveg);
	    
  //printf("surface_fluxes rec%d iveg%d %lf %lf %lf %lf %lf\n",
//	 rec,iveg,*runoff_wet,*runoff_dry,*baseflow_wet,*baseflow_dry,local_water);
}

