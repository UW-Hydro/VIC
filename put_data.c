#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
	      veg_con_struct    *veg_con,
              outfiles_struct   *outfiles,
              double            *depth,
	      double            *dz,
	      double             dp,
	      double            *AreaFract,
	      char              *AboveTreeLine,
	      dmy_struct        *dmy,
              int                rec,
	      int                dt,
	      int                Nnodes,
	      int                skipyear)
/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed presipitation data structures KAC
  01-20-00 modified to deal with simplified frozen soil moisture layers
           and frost depth / thaw depth accounting                 KAC
  03-08-00 modified to eliminate extra lines for storing bare
           soil variables.                                         KAC
  03-12-03 modified to add additional energy balance variable storage 
           when output of snow bands is selected.                  KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  28-Sep-04 Replaced aero_resist[0] with aero_resist_used, the aerodynamic
	    resistance that was actually used in flux calculations.	TJB

**********************************************************************/
{
  extern veg_lib_struct  *veg_lib;
  extern option_struct    options;
#if LINK_DEBUG
  extern debug_struct     debug;
#endif

#if OPTIMIZE
  static int              prtdt;
  static double           runoff;
  static double           baseflow;
#endif

  out_data_struct        *out_data;

  int                     veg;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  double                  Cv;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
  double                  tmp_ice;
  double                  rad_temp;
  double                  surf_temp;
  double                  inflow;
  double                  outflow;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  // Compute treeline adjustment factors
  for ( band = 0; band < Nbands; band++ ) {
    if ( AboveTreeLine[band] ) {
      Cv = 0;
      for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++ ) {
	if ( veg_lib[veg_con[veg].veg_class].overstory )
	  Cv += veg_con[veg].Cv;
      }
      TreeAdjustFactor[band] = 1. / ( 1. - Cv );
    }
    else TreeAdjustFactor[band] = 1.;
    if ( TreeAdjustFactor[band] != 1 && rec == 0 )
      fprintf( stderr, "WARNING: Tree adjust factor for band %i is equal to %f.\n", band, TreeAdjustFactor[band] );
  }

  /** Initialize Output Data Array **/
  out_data = (out_data_struct *) calloc(1,sizeof(out_data_struct)); 
  
  out_data->snow_canopy[0] = 0.0;
  out_data->prec           = atmos->out_prec;
  out_data->wind           = atmos->wind[NR];
  out_data->air_temp       = atmos->air_temp[NR];
  out_data->rel_humid      = 100.*atmos->vp[NR]/(atmos->vp[NR]+atmos->vpd[NR]);
  out_data->surf_temp      = 0.;
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
  /****************************************
    Store Output for all Vegetation Types
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {
    
    if ( veg < veg_con[0].vegetat_type_num ) 
      Cv = veg_con[veg].Cv;
    else
      Cv = (1.0 - veg_con[0].Cv_sum);

    if ( Cv > 0 ) {

      /*******************************************************
        Compute Average Variables from Wet and Dry Fractions
      *******************************************************/
      for ( dist = 0; dist < Ndist; dist++ ) {
	if(dist==0) 
	  mu = prcp[0].mu[veg];
	else 
	  mu = 1. - prcp[0].mu[veg];

	/*********************************
          Record Water Balance Variables 
	*********************************/

	/** record total evaporation **/
	for(band=0;band<Nbands;band++) {
	  if(AreaFract[band] > 0. && ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !veg_lib[veg_con[veg].veg_class].overstory)))) {
	    
#if !OPTIMIZE
	    
	    tmp_evap = 0.0;
	    for(index=0;index<options.Nlayer;index++)
	      tmp_evap += cell[dist][veg][band].layer[index].evap;
	    if ( veg < veg_con[0].vegetat_type_num )
	      out_data->evap_veg += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    else 
	      out_data->evap_bare += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];

	    tmp_evap += snow[veg][band].vapor_flux * 1000.;
	    out_data->sub_snow += snow[veg][band].vapor_flux * 1000. 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    if ( veg <= veg_con[0].vegetat_type_num ) {
	      tmp_evap += snow[veg][band].canopy_vapor_flux * 1000.;
	      out_data->sub_canop += snow[veg][band].canopy_vapor_flux 
		* 1000. * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    }
	    if ( veg < veg_con[0].vegetat_type_num ) {
	      tmp_evap += veg_var[dist][veg][band].canopyevap;
	      out_data->evap_canop += veg_var[dist][veg][band].canopyevap 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    }
	    out_data->evap += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 

#endif
	  
	    /** record runoff **/
	    out_data->runoff   += cell[dist][veg][band].runoff 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** record baseflow **/
	    out_data->baseflow += cell[dist][veg][band].baseflow 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    
#if ! OPTIMIZE
	    
	    /** record inflow **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data->inflow += (cell[dist][veg][band].inflow 
				   + veg_var[dist][veg][band].canopyevap) 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    else 
	      out_data->inflow += (cell[dist][veg][band].inflow) 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** record canopy interception **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data->Wdew += veg_var[dist][veg][band].Wdew 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	  
	    /** record aerodynamic resistance **/
	    out_data->aero_resist += cell[WET][veg][0].aero_resist_used
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** recored layer moistures **/
	    for(index=0;index<options.Nlayer;index++) {
	      tmp_moist = cell[dist][veg][band].layer[index].moist;
	      tmp_ice   = cell[dist][veg][band].layer[index].ice;
	      tmp_moist -= tmp_ice;
	      if(options.MOISTFRACT) {
		tmp_moist /= depth[index] * 1000.;
		tmp_ice /= depth[index] * 1000.;
	      }
	      tmp_moist *= Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      tmp_ice *= Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      out_data->moist[index] += tmp_moist;
	      out_data->ice[index]   += tmp_ice;
	    }
#endif
	  }
	}
      }

#if !OPTIMIZE
      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0. && ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !veg_lib[veg_con[veg].veg_class].overstory)))) {

	  /**********************************
	    Record Frozen Soil Variables
	  **********************************/

	  /** record freezing and thawing front depths **/
	  if(options.FROZEN_SOIL) {
	    for(index = 0; index < MAX_FRONTS; index++) {
	      if(energy[veg][band].fdepth[index] != MISSING)
		out_data->fdepth[index] += energy[veg][band].fdepth[index] 
		  * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	      if(energy[veg][band].tdepth[index] != MISSING)
		out_data->tdepth[index] += energy[veg][band].tdepth[index] 
		  * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	    }
	  }
	  
	  /**********************************
            Record Energy Balance Variables
	  **********************************/

	  /** record surface radiative temperature **/
	  if(snow[veg][band].swq>0) 
	    rad_temp = snow[veg][band].surf_temp + KELVIN;
	  else
	    rad_temp = energy[veg][band].T[0] + KELVIN;
	  
	  /** record soil surface temperature **/
	  surf_temp = energy[veg][band].T[0];
	  // MODIFIED FOR USE WITH ROSEMOUNT SIMULATIONS
	  //surf_temp = (energy[veg][band].T[0] + energy[veg][band].T[1])/2.;
	  
	  /** record net shortwave radiation **/
	  out_data->net_short += energy[veg][band].shortwave
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->swband[band] = energy[veg][band].shortwave
	      * Cv;
	  
	  /** record net longwave radiation **/
	  out_data->net_long  += energy[veg][band].longwave
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->lwband[band] = energy[veg][band].longwave
	      * Cv;
	  
	  /** record incoming longwave radiation **/
	  out_data->in_long  += ((energy[veg][band].longwave + STEFAN_B 
				  * (rad_temp) * (rad_temp)
				  * (rad_temp) * (rad_temp))
				 * Cv * AreaFract[band] * TreeAdjustFactor[band]);
	  
	  /** record albedo **/
	  out_data->albedo    += energy[veg][band].albedo
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->albedoband[band] = energy[veg][band].albedo
	      * Cv;
	  
	  /** record latent heat flux **/
	  out_data->latent    -= energy[veg][band].latent
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->latentband[band] = energy[veg][band].latent
	      * Cv;

	  /** record sensible heat flux **/
	  out_data->sensible  -= energy[veg][band].sensible
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->sensibleband[band] = energy[veg][band].sensible
	      * Cv;

	  /** record ground heat flux (+ heat storage) **/
	  out_data->grnd_flux -= (energy[veg][band].grnd_flux
				  + energy[veg][band].deltaH)
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  if(options.PRT_SNOW_BAND) 
	    out_data->grndband[band] = -1.*(energy[veg][band].grnd_flux
					    + energy[veg][band].deltaH)
	      * Cv;

	  /** record heat storage **/
	  out_data->deltaH    -= energy[veg][band].deltaH
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record energy balance error **/
	  out_data->energy_error += energy[veg][band].error
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record radiative effective temperature [K], 
	      emissivities set = 1.0  **/
	  out_data->rad_temp += ((rad_temp) * (rad_temp) 
				 * (rad_temp) * (rad_temp)) 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record mean surface temperature [C]  **/
	  out_data->surf_temp += surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /*****************************
	    Record Snow Pack Variables 
	  *****************************/
	  
	  /** record snow water equivalence **/
	  out_data->swq[0]         
	    += snow[veg][band].swq * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack depth **/
	  out_data->snow_depth[0]  
	    += snow[veg][band].depth * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record canopy intercepted snow **/
	  if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	    out_data->swq[0] 
#else
	      out_data->snow_canopy[0] 
#endif
	      += (snow[veg][band].snow_canopy) 
	      * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow cover fraction **/
	  out_data->coverage[0]    
	    += snow[veg][band].coverage * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack cold content **/
	  out_data->deltaCC[0]           
	    += energy[veg][band].deltaCC * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack advection **/
	  out_data->advection[0]         
	    += energy[veg][band].advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow energy flux **/
	  out_data->snow_flux[0]         
	    += energy[veg][band].snow_flux * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record refreeze energy **/
	  out_data->refreeze_energy[0]   
	    += energy[veg][band].refreeze_energy 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** if snow elevation bands are to be printed separately **/
	  if(options.PRT_SNOW_BAND) {
	    
	    /** record band snow water equivalent **/
	    out_data->swq[band+1]         
	      += snow[veg][band].swq * Cv  * 1000.;
	    
	    /** record band snowpack depth **/
	    out_data->snow_depth[band+1]  
	      += snow[veg][band].depth * Cv * 100.;
	    
	    /** record band canopy intercepted snow **/
	    if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	      out_data->swq[band+1]
#else
		out_data->snow_canopy[band+1]  
#endif
		+= (snow[veg][band].snow_canopy) * Cv * 1000.;
	    
	    /** record band snow coverage **/
	    out_data->coverage[band+1]    
	      += snow[veg][band].coverage * Cv;
	    
	    /** record band cold content **/
	    out_data->deltaCC[band+1]           
	      += energy[veg][band].deltaCC * Cv;
	    
	    /** record band advection **/
	    out_data->advection[band+1]         
	      += energy[veg][band].advection * Cv;
	    
	    /** record band snow flux **/
	    out_data->snow_flux[band+1]         
	      += energy[veg][band].snow_flux * Cv;
	    
	    /** record band refreeze energy **/
	    out_data->refreeze_energy[band+1]   
	      += energy[veg][band].refreeze_energy * Cv;
	    
	  }
	}
      }

#endif /* not OPTIMIZE */

    }
  }
  
#if !OPTIMIZE

  /** record radiative temperature **/
  out_data->rad_temp = pow(out_data->rad_temp,0.25);

  /** record net radiation **/
  out_data->r_net    = out_data->net_short + out_data->net_long;

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data->prec;
  outflow = out_data->evap+out_data->runoff+out_data->baseflow;
  storage = 0.;
  for(index=0;index<options.Nlayer;index++)
    if(options.MOISTFRACT)
      storage += (out_data->moist[index] + out_data->ice[index]) 
	* depth[index] * 1000;
    else
      storage += out_data->moist[index] + out_data->ice[index];
  storage += out_data->swq[0] + out_data->snow_canopy[0] + out_data->Wdew;
  calc_water_balance_error(rec,inflow,outflow,storage);
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data->net_short + out_data->net_long,
			      out_data->latent, out_data->sensible,
			      out_data->grnd_flux, out_data->advection[0] 
			      - out_data->deltaCC[0] - out_data->snow_flux[0]
			      + out_data->refreeze_energy[0]);

  /*************
    Write Data
  *************/

  if(rec >= skipyear)
    write_data(out_data, outfiles, dmy, dt);

#else

  if ( rec == 0 ) prtdt = 0;
  if ( prtdt == 0 ) {
    runoff = out_data->runoff;
    baseflow = out_data->baseflow;
    prtdt ++;
  }
  else {
    runoff += out_data->runoff;
    baseflow += out_data->baseflow;
    prtdt ++;
  }
  if ( prtdt == 24 / dt ) {
    out_data->runoff = runoff;
    out_data->baseflow = baseflow;
    write_data(out_data, outfiles, dmy, dt);
    prtdt = 0;
  } 

#endif

  free((char *)out_data); 

}
