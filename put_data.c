#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
	      veg_con_struct    *veg_con,
              outfiles_struct   outfiles,
              double            *depth,
	      dmy_struct        *dmy,
              int               rec)

/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct  options;
  extern debug_struct   debug;

  static out_data_struct *out_data;

  int veg, vegnum;
  int index;
  int Ndist;
  int dist;
  int dt;
  double mu;
  double tmp_evap;
  double tmp_moist;
  double tmp_ice;
  double rad_temp = 0.0;
  double inflow, outflow, storage;

  cell_data_struct *cell;
  snow_data_struct *snow;
  energy_bal_struct *energy;
  veg_var_struct *veg_var;

  out_data = (out_data_struct *) calloc(1,sizeof(out_data_struct));

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;

  out_data->prec = atmos->prec;
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  for ( dist = 0; dist < Ndist; dist++ ) {
    if(dist==0) 
      mu = prcp[0].mu;
    else 
      mu = 1. - prcp[0].mu;

    cell = prcp[0].dist[dist].cell;
    snow = prcp[0].dist[dist].snow;
    energy = prcp[0].dist[dist].energy;
    veg_var = prcp[0].dist[dist].veg_var;

    /**************************************
      Store Output for all Vegetation Types
      **************************************/
    for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++) {

      /** record total evaporation **/
      tmp_evap = 0.0;
      for(index=0;index<options.Nlayer;index++)
        tmp_evap += cell[veg].layer[index].evap;
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        tmp_evap += snow[veg].vapor_flux * 1000.;
        tmp_evap += snow[veg].canopy_vapor_flux * 1000.;
      }
      out_data->evap += (tmp_evap + veg_var[veg].canopyevap)
                      * veg_con[veg].Cv * mu; 

      /** record runoff **/
      out_data->runoff   += cell[veg].runoff * veg_con[veg].Cv * mu;

      /** record baseflow **/
      out_data->baseflow += cell[veg].baseflow * veg_con[veg].Cv * mu;
  
      /** record inflow **/
      out_data->inflow += (cell[veg].inflow + veg_var[veg].canopyevap) 
                         * veg_con[veg].Cv * mu;

      /** record canopy interception **/
      out_data->Wdew   += veg_var[veg].Wdew * veg_con[veg].Cv * mu;

      /** record aerodynamic resistance **/
      out_data->aero_resist += cell[veg].aero_resist[0] * veg_con[veg].Cv * mu;

      /** recored layer moistures **/
      for(index=0;index<options.Nlayer;index++) {
        tmp_moist = (cell[veg].layer[index].moist_thaw 
                    * cell[veg].layer[index].tdepth
                    / depth[index]);
        tmp_moist += (cell[veg].layer[index].moist_froz 
                     * (cell[veg].layer[index].fdepth
                     - cell[veg].layer[index].tdepth)
                     / depth[index]);
        tmp_ice = (cell[veg].layer[index].ice 
                  * (cell[veg].layer[index].fdepth
                  - cell[veg].layer[index].tdepth)
                  / depth[index]);
        tmp_moist += (cell[veg].layer[index].moist 
                     * (depth[index]
                     - cell[veg].layer[index].fdepth)
                     / depth[index]);
        if(options.MOISTFRACT) {
          tmp_moist /= depth[index] * 1000.;
          tmp_ice /= depth[index] * 1000.;
        }
        tmp_moist *= veg_con[veg].Cv * mu;
        tmp_ice *= veg_con[veg].Cv * mu;
        out_data->moist[index] += tmp_moist;
        out_data->ice[index]   += tmp_ice;
      }
  
      /** record freezing and thawing front depths **/
      if(options.FROZEN_SOIL) {
        for(index=0;index<2;index++)
          out_data->fdepth[index] += energy[veg].fdepth[index] 
                                    * veg_con[veg].Cv * mu * 100.;
      }

      /** record surface fluxes **/
      if(options.FULL_ENERGY) {
        out_data->net_short += energy[veg].shortwave
                          * (1.0 - energy[veg].albedo)
                          * veg_con[veg].Cv * mu;
        out_data->net_long  += energy[veg].longwave
                          * veg_con[veg].Cv * mu;
        out_data->albedo    += energy[veg].albedo
                          * veg_con[veg].Cv * mu;
        out_data->latent    += energy[veg].latent
                          * veg_con[veg].Cv * mu;
        out_data->sensible  += energy[veg].sensible
                          * veg_con[veg].Cv * mu;
        out_data->grnd_flux += (energy[veg].grnd_flux
                          + energy[veg].deltaH)
                          * veg_con[veg].Cv * mu;
        out_data->deltaH    += energy[veg].deltaH
                          * veg_con[veg].Cv * mu;
        out_data->energy_error += energy[veg].error
                          * veg_con[veg].Cv * mu;
        out_data->surf_temp += energy[veg].T[0]
                          * veg_con[veg].Cv * mu;

        /** record radiative effective temperature [K], emissivities set = 1.0  **/
        rad_temp            += pow((273.15+energy->T[0]),4.0) * veg_con[veg].Cv * mu;

        /** record mean surface temperature [C]  **/
        out_data->surf_temp += energy->T[0] * veg_con[veg].Cv * mu;

      }
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        out_data->swq         += snow[veg].swq * veg_con[veg].Cv * mu * 1000.;
        out_data->snow_depth  += snow[veg].depth * veg_con[veg].Cv
	                       * mu * 1000.;
        out_data->snow_canopy += (snow[veg].snow_canopy
                               + snow[veg].tmp_int_storage) * veg_con[veg].Cv
	                       * mu * 1000.;
        out_data->deltaCC     += energy[veg].deltaCC * veg_con[veg].Cv * mu;
        out_data->advection   += energy[veg].advection * veg_con[veg].Cv * mu;
        out_data->snow_flux   += energy[veg].snow_flux * veg_con[veg].Cv * mu;
        out_data->refreeze_energy   += energy[veg].refreeze_energy 
	                             * veg_con[veg].Cv * mu;
      }
    }
  
    /***********************
      Store Bare Soil Output
      ***********************/

    /** record evaporation for bare soil **/
    vegnum = veg_con[0].vegetat_type_num;
    tmp_evap=0.;
    for(index=0;index<options.Nlayer;index++)
      tmp_evap += cell[vegnum].layer[index].evap;
    if(options.FULL_ENERGY || options.SNOW_MODEL)
      tmp_evap += snow[veg].vapor_flux * 1000.;
    out_data->evap += tmp_evap * (1.0 - veg_con[0].Cv_sum) * mu;
  

    /** record runoff for bare soil **/
    out_data->runoff   += cell[vegnum].runoff * (1.0 - veg_con[0].Cv_sum) * mu;

    /** record baseflow for bare soil **/
    out_data->baseflow += cell[vegnum].baseflow * (1.0 - veg_con[0].Cv_sum) * mu;

    /** record inflow for bare soil **/
    out_data->inflow += cell[vegnum].inflow * (1.0 - veg_con[0].Cv_sum) * mu;

    /** record aerodynamic resistance **/
    out_data->aero_resist += cell[vegnum].aero_resist[0]
                           * (1.0 - veg_con[veg].Cv_sum) * mu;

    /** record layer moisture for bare soil **/
    for(index=0;index<options.Nlayer;index++) {
      tmp_moist  = (cell[vegnum].layer[index].moist_thaw 
                    * cell[vegnum].layer[index].tdepth / depth[index]);
      tmp_moist += (cell[vegnum].layer[index].moist_froz 
                    * (cell[vegnum].layer[index].fdepth 
                    - cell[vegnum].layer[index].tdepth) / depth[index]);
      tmp_moist += (cell[vegnum].layer[index].moist * (depth[index]
                    - cell[vegnum].layer[index].fdepth) / depth[index]);
      tmp_ice    = (cell[vegnum].layer[index].ice 
                    * (cell[vegnum].layer[index].fdepth
                    - cell[vegnum].layer[index].tdepth)
                    / depth[index]);
      if(options.MOISTFRACT) {
        tmp_moist /= depth[index] * 1000.;
        tmp_ice /= depth[index] * 1000.;
      }
      tmp_moist              *= (1.0 - veg_con[0].Cv_sum) * mu;
      tmp_ice                *= (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->moist[index] += tmp_moist;
      out_data->ice[index]   += tmp_ice;
    }
  
    /** record freezing and thawing depths for bare soil **/
    if(options.FROZEN_SOIL) {
      for(index=0;index<2;index++) 
	out_data->fdepth[index] += energy[vegnum].fdepth[index] 
                                   * (1.0 - veg_con[0].Cv_sum) * mu * 100.;
    }
  
    /** record surface fluxes for bare soil **/
    if(options.FULL_ENERGY) {
      out_data->net_short += energy[vegnum].shortwave 
	                    * (1.0 - energy[vegnum].albedo)
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->net_long  += energy[vegnum].longwave 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->albedo    += energy[vegnum].albedo
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->latent    += energy[vegnum].latent 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->sensible  += energy[vegnum].sensible 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->grnd_flux += (energy[vegnum].grnd_flux 
                            + energy[vegnum].deltaH)
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->deltaH    += energy[vegnum].deltaH 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->energy_error += energy[vegnum].error 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->surf_temp += energy[vegnum].T[0] 
                            * (1.0 - veg_con[0].Cv_sum) * mu;
    }
    if(options.FULL_ENERGY || options.SNOW_MODEL) {
      out_data->swq         += snow[vegnum].swq * (1.0 - veg_con[0].Cv_sum) 
                              * mu * 1000.;
      out_data->snow_depth  += snow[vegnum].depth 
                              * (1.0 - veg_con[0].Cv_sum) * mu * 1000.;
      out_data->deltaCC += energy[vegnum].deltaCC 
                              * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->advection   += energy[vegnum].advection   
                              * (1.0 - veg_con[0].Cv_sum) * mu;
      out_data->snow_flux   += energy[vegnum].snow_flux   
                              * (1.0 - veg_con[0].Cv_sum) * mu;
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
  storage += out_data->swq + out_data->snow_canopy + out_data->Wdew;
  calc_water_balance_error(rec,inflow,outflow,storage);
  calc_energy_balance_error(rec,out_data->net_short+out_data->net_long,
			    out_data->latent,out_data->sensible,
			    out_data->grnd_flux,
			    out_data->advection-out_data->deltaCC
			    -out_data->snow_flux+out_data->refreeze_energy);

  write_data(out_data, outfiles, dmy);

  free((char *) out_data); 
}
