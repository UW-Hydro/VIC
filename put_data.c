#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
	      veg_con_struct    *veg_con,
              outfiles_struct    outfiles,
              double            *depth,
	      double            *AreaFract,
	      dmy_struct        *dmy,
              int                rec)
/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed presipitation data structures KAC

**********************************************************************/
{
  extern option_struct    options;
  extern debug_struct     debug;

  out_data_struct        *out_data;

  int                     veg, vegnum;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
  double                  tmp_ice;
  double                  rad_temp = 0.0;
  double                  inflow, outflow, storage;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = options.SNOW_BAND;

  out_data = (out_data_struct *) calloc(1,sizeof(out_data_struct));
  out_data->swq             = (double *)calloc(Nbands+1,sizeof(double));
  out_data->snow_canopy     = (double *)calloc(Nbands+1,sizeof(double));
  out_data->snow_depth      = (double *)calloc(Nbands+1,sizeof(double));
  out_data->advection       = (double *)calloc(Nbands+1,sizeof(double));
  out_data->deltaCC         = (double *)calloc(Nbands+1,sizeof(double));
  out_data->snow_flux       = (double *)calloc(Nbands+1,sizeof(double));
  out_data->refreeze_energy = (double *)calloc(Nbands+1,sizeof(double));
  out_data->coverage        = (double *)calloc(Nbands+1,sizeof(double));

  out_data->prec = atmos->prec;
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp[0].cell;
  veg_var = prcp[0].veg_var;
  snow    = prcp[0].snow;
  energy  = prcp[0].energy;
  
  /**************************************
    Store Output for all Vegetation Types
    **************************************/
  for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++) {
    
    /** record total evaporation **/
    for ( dist = 0; dist < Ndist; dist++ ) {
      if(dist==0) 
	mu = prcp[0].mu[veg];
      else 
	mu = 1. - prcp[0].mu[veg];

      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0.) {
	  tmp_evap = 0.0;
	  for(index=0;index<options.Nlayer;index++)
	    tmp_evap += cell[dist][veg][band].layer[index].evap;
	  if(options.FULL_ENERGY || options.SNOW_MODEL) {
	    tmp_evap += snow[veg][band].vapor_flux * 1000.;
	    tmp_evap += snow[veg][band].canopy_vapor_flux * 1000.;
	  }
	  tmp_evap += veg_var[dist][veg][band].canopyevap;
	  out_data->evap += tmp_evap * veg_con[veg].Cv * mu * AreaFract[band]; 
	  
	  /** record runoff **/
	  out_data->runoff   += cell[dist][veg][band].runoff 
	    * veg_con[veg].Cv * mu * AreaFract[band];
	  
	  /** record baseflow **/
	  out_data->baseflow += cell[dist][veg][band].baseflow 
	    * veg_con[veg].Cv * mu * AreaFract[band]; 
	  
	  /** record inflow **/
	  out_data->inflow += (cell[dist][veg][band].inflow 
			       + veg_var[dist][veg][band].canopyevap) 
	    * veg_con[veg].Cv * mu * AreaFract[band];
	  
	  /** record canopy interception **/
	  out_data->Wdew += veg_var[dist][veg][band].Wdew * veg_con[veg].Cv 
	    * mu * AreaFract[band];
	  
	  /** record aerodynamic resistance **/
	  out_data->aero_resist += cell[WET][veg][0].aero_resist[0] 
	    * veg_con[veg].Cv * mu * AreaFract[band];
	  
	  /** recored layer moistures **/
	  for(index=0;index<options.Nlayer;index++) {
	    tmp_moist 
	      = find_total_layer_moisture(cell[dist][veg][band].layer[index], 
					  depth[index]);
	    tmp_ice 
	      = find_total_layer_ice(cell[dist][veg][band].layer[index], 
					  depth[index]);
	    tmp_moist -= tmp_ice;
	    if(options.MOISTFRACT) {
	      tmp_moist /= depth[index] * 1000.;
	      tmp_ice /= depth[index] * 1000.;
	    }
	    tmp_moist *= veg_con[veg].Cv * mu * AreaFract[band];
	    tmp_ice *= veg_con[veg].Cv * mu * AreaFract[band];
	    out_data->moist[index] += tmp_moist;
	    out_data->ice[index]   += tmp_ice;
	  }
	}
      }
    }
    
    for(band=0;band<Nbands;band++) {
      if(AreaFract[band]>0.) {

	/** record freezing and thawing front depths **/
	if(options.FROZEN_SOIL) {
	  for(index=0;index<2;index++)
	    out_data->fdepth[index] += energy[veg][band].fdepth[index] 
	      * veg_con[veg].Cv * 100. * AreaFract[band];
	}
	
	/** record surface fluxes **/
	if(options.FULL_ENERGY) {
	  out_data->net_short += energy[veg][band].shortwave
	    * (1.0 - energy[veg][band].albedo)
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->net_long  += energy[veg][band].longwave
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->albedo    += energy[veg][band].albedo
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->latent    -= energy[veg][band].latent
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->sensible  -= energy[veg][band].sensible
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->grnd_flux -= (energy[veg][band].grnd_flux
				  + energy[veg][band].deltaH)
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->latent_pet -= (energy[veg][band].grnd_flux)
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->latent_pet_mm -= (energy[veg][band].deltaH)
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->deltaH    -= energy[veg][band].deltaH
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->energy_error += energy[veg][band].error
	    * veg_con[veg].Cv * AreaFract[band];
	
	  /** record radiative effective temperature [K], emissivities set = 1.0  **/
	  rad_temp            += pow((273.15+energy[veg][band].T[0]),4.0) 
	    * veg_con[veg].Cv * AreaFract[band];
	  
	  /** record mean surface temperature [C]  **/
	  out_data->surf_temp += energy[veg][band].T[0] * veg_con[veg].Cv 
	    * AreaFract[band];
	  
	}
	if(options.FULL_ENERGY || options.SNOW_MODEL) {
	  out_data->swq[0]         
	    += snow[veg][band].swq * veg_con[veg].Cv 
	    * 1000. * AreaFract[band];
	  out_data->snow_depth[0]  
	    += snow[veg][band].depth * veg_con[veg].Cv
	    * 100. * AreaFract[band];
	  out_data->snow_canopy[0] 
	    += (snow[veg][band].snow_canopy) 
	    * veg_con[veg].Cv * 1000. * AreaFract[band];
	  out_data->coverage[0]    
	    += snow[veg][band].coverage 
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->deltaCC[0]           
	    += energy[veg][band].deltaCC 
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->advection[0]         
	    += energy[veg][band].advection 
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->snow_flux[0]         
	    += energy[veg][band].snow_flux 
	    * veg_con[veg].Cv * AreaFract[band];
	  out_data->refreeze_energy[0]   
	    += energy[veg][band].refreeze_energy 
	    * veg_con[veg].Cv * AreaFract[band];
	  if(options.PRT_SNOW_BAND) {
	    out_data->swq[band+1]         
	      += snow[veg][band].swq * veg_con[veg].Cv  * 1000.;
	    out_data->snow_depth[band+1]  
	      += snow[veg][band].depth * veg_con[veg].Cv * 100.;
	    out_data->snow_canopy[band+1] 
	      += (snow[veg][band].snow_canopy) 
	      * veg_con[veg].Cv * 1000.;
	    out_data->coverage[band+1]    
	      += snow[veg][band].coverage 
	      * veg_con[veg].Cv;
	    out_data->deltaCC[band+1]           
	      += energy[veg][band].deltaCC 
	      * veg_con[veg].Cv;
	    out_data->advection[band+1]         
	      += energy[veg][band].advection 
	      * veg_con[veg].Cv;
	    out_data->snow_flux[band+1]         
	      += energy[veg][band].snow_flux 
	      * veg_con[veg].Cv;
	    out_data->refreeze_energy[band+1]   
	      += energy[veg][band].refreeze_energy 
	      * veg_con[veg].Cv;
	  }
	}
      }
    }
  }
  
  /***********************
    Store Bare Soil Output
    ***********************/

  /** record evaporation for bare soil **/
  vegnum = veg_con[0].vegetat_type_num;
  for ( dist = 0; dist < Ndist; dist++ ) {
    if(dist==0) 
      mu = prcp[0].mu[veg];
    else 
      mu = 1. - prcp[0].mu[veg];

    for(band=0;band<Nbands;band++) {
      if(AreaFract[band]>0.) {
	tmp_evap=0.;
	for(index=0;index<options.Nlayer;index++)
	  tmp_evap += cell[dist][vegnum][band].layer[index].evap;
	if(options.FULL_ENERGY || options.SNOW_MODEL)
	  tmp_evap += snow[vegnum][band].vapor_flux * 1000.;
	out_data->evap += tmp_evap * (1.0 - veg_con[0].Cv_sum) * mu 
	  * AreaFract[band];
	
	/** record runoff for bare soil **/
	out_data->runoff   += cell[dist][vegnum][band].runoff 
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
	/** record baseflow for bare soil **/
	out_data->baseflow += cell[dist][vegnum][band].baseflow 
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
	/** record inflow for bare soil **/
	out_data->inflow += cell[dist][vegnum][band].inflow 
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
	/** record aerodynamic resistance **/
	out_data->aero_resist += cell[dist][vegnum][band].aero_resist[0]
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
	/** record layer moisture for bare soil **/
	for(index=0;index<options.Nlayer;index++) {
	  tmp_moist 
	    = find_total_layer_moisture(cell[dist][vegnum][band].layer[index], 
					depth[index]);
	  tmp_ice 
	    = find_total_layer_ice(cell[dist][vegnum][band].layer[index], 
				   depth[index]);
	  tmp_moist -= tmp_ice;
	  if(options.MOISTFRACT) {
	    tmp_moist /= depth[index] * 1000.;
	    tmp_ice /= depth[index] * 1000.;
	  }
	  tmp_moist              *= (1.0 - veg_con[0].Cv_sum) * mu 
	    * AreaFract[band];
	  tmp_ice                *= (1.0 - veg_con[0].Cv_sum) * mu 
	    * AreaFract[band];
	  out_data->moist[index] += tmp_moist;
	  out_data->ice[index]   += tmp_ice;
	}
      }
    }
  }

  for(band=0;band<Nbands;band++) {
    if(AreaFract[band]>0.) {

      /** record freezing and thawing depths for bare soil **/
      if(options.FROZEN_SOIL) {
	for(index=0;index<2;index++) 
	  out_data->fdepth[index] += energy[vegnum][band].fdepth[index] 
	    * (1.0 - veg_con[0].Cv_sum) * 100. * AreaFract[band];
      }
      
      /** record surface fluxes for bare soil **/
      if(options.FULL_ENERGY) {
	out_data->net_short += energy[vegnum][band].shortwave 
	  * (1.0 - energy[vegnum][band].albedo)
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->net_long  += energy[vegnum][band].longwave 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->albedo    += energy[vegnum][band].albedo
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->latent    -= energy[vegnum][band].latent 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->sensible  -= energy[vegnum][band].sensible 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->grnd_flux -= (energy[vegnum][band].grnd_flux 
				+ energy[vegnum][band].deltaH)
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->latent_pet -= (energy[vegnum][band].grnd_flux)
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->latent_pet_mm -= (energy[vegnum][band].deltaH)
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->deltaH    -= energy[vegnum][band].deltaH 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->energy_error += energy[vegnum][band].error 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->surf_temp += energy[vegnum][band].T[0] 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	rad_temp            += pow((273.15+energy[vegnum][band].T[0]),4.0) 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	
      }
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
	out_data->swq[0]             += snow[vegnum][band].swq 
	  * (1.0 - veg_con[0].Cv_sum) * 1000. * AreaFract[band];
	out_data->snow_depth[0]      += snow[vegnum][band].depth 
	  * (1.0 - veg_con[0].Cv_sum) * 100. * AreaFract[band];
	out_data->coverage[0]        += snow[vegnum][band].coverage 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->deltaCC[0] += energy[vegnum][band].deltaCC 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->advection[0]       += energy[vegnum][band].advection   
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->snow_flux[0]       += energy[vegnum][band].snow_flux   
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	out_data->refreeze_energy[0] += energy[vegnum][band].refreeze_energy 
	  * (1.0 - veg_con[0].Cv_sum) * AreaFract[band];
	if(options.PRT_SNOW_BAND) {
	  out_data->swq[band+1]             += snow[vegnum][band].swq 
	    * (1.0 - veg_con[0].Cv_sum) * 1000.;
	  out_data->snow_depth[band+1]      += snow[vegnum][band].depth 
	    * (1.0 - veg_con[0].Cv_sum) * 100.;
	  out_data->coverage[band+1]        += snow[vegnum][band].coverage 
	    * (1.0 - veg_con[0].Cv_sum);
	  out_data->deltaCC[band+1] += energy[vegnum][band].deltaCC 
	    * (1.0 - veg_con[0].Cv_sum);
	  out_data->advection[band+1]       
	    += energy[vegnum][band].advection   
	    * (1.0 - veg_con[0].Cv_sum);
	  out_data->snow_flux[band+1]       
	    += energy[vegnum][band].snow_flux   
	    * (1.0 - veg_con[0].Cv_sum);
	  out_data->refreeze_energy[band+1] 
	    += energy[vegnum][band].refreeze_energy 
	    * (1.0 - veg_con[0].Cv_sum);
	}
      }
    }
  }


  out_data->rad_temp = pow(rad_temp,0.25);
  out_data->r_net    = out_data->net_short + out_data->net_long;

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data->prec;
  outflow = out_data->evap+out_data->runoff+out_data->baseflow;
  storage = 0.;
  for(index=0;index<options.Nlayer;index++)
    storage += out_data->moist[index] + out_data->ice[index];
  storage += out_data->swq[0] + out_data->snow_canopy[0] + out_data->Wdew;
  calc_water_balance_error(rec,inflow,outflow,storage);
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec,out_data->net_short+out_data->net_long,
			      out_data->latent,out_data->sensible,
			      out_data->grnd_flux,
			      out_data->advection[0]-out_data->deltaCC[0]
			      -out_data->snow_flux[0]
			      +out_data->refreeze_energy[0]);

  write_data(out_data, outfiles, dmy);

  free((char *)out_data->swq);
  free((char *)out_data->snow_canopy);
  free((char *)out_data->snow_depth);
  free((char *)out_data->advection);
  free((char *)out_data->deltaCC);
  free((char *)out_data->snow_flux);
  free((char *)out_data->refreeze_energy);
  free((char *)out_data->coverage);
  free((char *)out_data); 
}
