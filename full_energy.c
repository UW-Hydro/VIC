#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#define OVERSTORY_ATTENUATION 3.5
#define TRUNK_RATIO 0.5

void full_energy(int rec,
                 atmos_data_struct *atmos,
                 soil_con_struct soil_con,
                 veg_con_struct *veg_con,
                 prcp_var_struct *dist,
                 dmy_struct *dmy,
                 global_param_struct gp,
                 double mu,
                 int prcpdist,
                 int gridcell,
                 char NEWCELL)
/**********************************************************************
	full_energy	Keith Cherkauer		January 8, 1997

  This subroutine controls the model core, it solves both the energy
  and water balance models, as well as frozen soils.  

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern debug_struct debug;

  static double tmax[3];
  static double tmin[3];

  char   RUN_TYPE;
  char   DONE;
  char   snow_present = FALSE;
  char   overstory;
  int    SNOW;
  int    iter;
  int    i;
  int    iveg;
  int    Nveg;
  int    veg_class;
  int    hour;
  double ppt;		/** Precipitation **/
  double surf_temp;
  double tmp_surf_temp;
  double last_T1;
  double T0, new_T0;
  double out_short;	/* outgoing shortwave radiation (W/m^2) */
  double inshort;
  double inlong;
  double dp;
  double last_moist;
  double grnd_flux;	/* temporary storage for ground flux under snow pack */
  double deltaH;
  double ice0, ice, moist;
  double surf_atten;
  double threshold = 0.001;
  double tmp_swq;
  double x;
  double rainfall;
  double snowfall;
  double height;
  double displacement;
  double roughness;
  double ref_height;
  double dT;
  double *null_ptr;
  double dummy;
  double snow_depth = 0.;
  double Ls;
  double canopy_temp;
  double tmp_rain;
  double C1, C2, C3;
  double new_T1;
  double kappa_snow;
  double snow_flux;
  double throughfall;
  double Evap;
  double tmp_throughfall;
  double tmp_shortwave;
  double tmp_longwave;
  double tmp_rad;
  double tmp_rainfall;
  double tmp_snowfall;
  double tmp_wind[3];
  double tmp_aero_resist;
  double tmp_melt;
  double tmp_vapor_flux;
  double tmp_canopy_vapor_flux;
  double tmp_canopyevap;
  double tmp_snowinflow;
  double tmp_snow_surf_temp;
  double tmp_snow_energy;
  double *tmp_layerevap;
  atmos_data_struct tmp_atmos;
  layer_data_struct *tmp_layer;
  veg_var_struct    tmp_veg_var;
  cell_data_struct  *cell;
  veg_var_struct    *veg_var;
  energy_bal_struct *energy;
  energy_bal_struct *tmp_energy;
  snow_data_struct  *snow;
  snow_data_struct  *tmp_snow;
  veg_var_struct    *tmp_veg;

  cell = dist[0].cell;
  veg_var = dist[0].veg_var;
  snow = dist[0].snow;
  energy = dist[0].energy;

  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    tmp_layer = (layer_data_struct *)calloc((options.Nlayer),
        sizeof(layer_data_struct));
  }

  if(gp.dt == 24 && options.SNOW_MODEL) {
    store_max_min_temp(atmos,tmax,tmin,rec,gp.nrecs,prcpdist);
  }
      
  Nveg = veg_con[0].vegetat_type_num;
  tmp_longwave = atmos->longwave;
  tmp_atmos = atmos[0];

  /**************************************************
    Solve Energy and/or Water Balance for Each
    Vegetation Type
  **************************************************/
  for(iveg=0;iveg<=Nveg;iveg++){

    /** Solve Veg Type only if Coverage Greater than 0% **/
    if(iveg<Nveg) {
      if(veg_con[iveg].Cv>0.0) RUN_TYPE = 1;
      else RUN_TYPE = 0;
    }
    else if(iveg==Nveg && veg_con[0].Cv_sum<1.) RUN_TYPE = 1;
    else RUN_TYPE = 0;
    
    if(RUN_TYPE) {

#ifdef ROSEMOUNT
      if(dmy[rec].day == 26 && dmy[rec].month == 10 && dmy[rec].year == 1995
         && dmy[rec].hour == 0) {
        cell[iveg].layer[0].moist = 0.3071 * soil_con.depth[0] * 1000.;
        cell[iveg].layer[1].moist = 0.4249 * soil_con.depth[1] * 1000.;
      }
#endif

      /**************************************************
        Initialize Model Parameters
      **************************************************/
      for(i=0;i<options.Nlayer;i++) cell[iveg].layer[i].evap = 0.;
      if(iveg<Nveg) veg_var[iveg].canopyevap = 0.;

      tmp_rain = 0.;
  
      x = 1.;

      if(iveg<Nveg) veg_class = veg_con[iveg].veg_class;
      else veg_class = 0;

      if(iveg!=Nveg && (!options.SNOW_MODEL
			|| (!snow[iveg].snow || veg_lib[veg_class].overstory)))
	  surf_atten=exp(-0.50*veg_lib[veg_class].LAI[dmy[rec].month-1]);
      else surf_atten=1;
        
      if(options.SNOW_MODEL) {
        if(snow[iveg].snow) Ls = (677. - 0.07 * snow[iveg].surf_temp)
                               * 4.1868 * 1000;
        else Ls = 0.;
      }
      else Ls = 0.;

      tmp_snowinflow = 0.;

      /**************************************************
        Store Water Balance Terms for Debugging
      **************************************************/
      if(debug.DEBUG || debug.PRT_MOIST || debug.PRT_BALANCE) {
        /** Compute current total moisture for water balance check **/
        debug.store_moist[options.Nlayer+2] = 0.;
        if(iveg<Nveg) {
          debug.store_moist[0] = veg_var[iveg].Wdew;
          if(options.FULL_ENERGY || options.SNOW_MODEL)
            debug.store_moist[0] += snow[iveg].snow_canopy*1000.;
          debug.store_moist[options.Nlayer+2] += debug.store_moist[0];
        }
        if(options.FULL_ENERGY || options.SNOW_MODEL) {
          debug.store_moist[1] = snow[iveg].swq*1000.;
          debug.store_moist[options.Nlayer+2] += debug.store_moist[1];
        }
        for(i=0;i<options.Nlayer;i++) {
          debug.store_moist[i+2] = cell[iveg].layer[i].moist_thaw
                      * cell[iveg].layer[i].tdepth 
               / soil_con.depth[i];
          debug.store_moist[i+2] += cell[iveg].layer[i].moist_froz
                      * (cell[iveg].layer[i].fdepth
                      - cell[iveg].layer[i].tdepth)
                      / soil_con.depth[i];
          debug.store_moist[i+2] += cell[iveg].layer[i].ice
                      * (cell[iveg].layer[i].fdepth 
               - cell[iveg].layer[i].tdepth) / soil_con.depth[i];
          debug.store_moist[i+2] += cell[iveg].layer[i].moist * (soil_con.depth[i] 
               - cell[iveg].layer[i].fdepth) / soil_con.depth[i];
          debug.store_moist[options.Nlayer+2] += debug.store_moist[i+2];
        }
      }

      /**************************************************
        Initialize Energy Balance Terms
      **************************************************/
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
	/** Clear Energy Balance Terms for Snow Pack **/
	energy[iveg].deltaCC         = 0.;
	energy[iveg].snow_flux       = 0.;
	energy[iveg].refreeze_energy = 0.;
	energy[iveg].advection       = 0.;

        /** Set Damping Depth **/
        dp=soil_con.dp;

        /** Compute Total Moisture of Surface Layer **/
        moist = cell[iveg].layer[0].moist_thaw * cell[iveg].layer[0].tdepth 
              / soil_con.depth[0];
        moist += cell[iveg].layer[0].moist_froz * (cell[iveg].layer[0].fdepth 
               - cell[iveg].layer[0].tdepth) / soil_con.depth[0];
        moist += cell[iveg].layer[0].ice * (cell[iveg].layer[0].fdepth 
               - cell[iveg].layer[0].tdepth) / soil_con.depth[0];
        moist += cell[iveg].layer[0].moist * (soil_con.depth[0] 
               - cell[iveg].layer[0].fdepth) / soil_con.depth[0];
        moist /= soil_con.depth[0]*1000.;
        if(options.FROZEN_SOIL){
          if((energy[iveg].T[0]+energy[iveg].T[1])/2.<0.) {
            ice0 = moist - maximum_unfrozen_water((energy[iveg].T[0]
	        + energy[iveg].T[1])/2.,
	          soil_con.max_moist[0]/(soil_con.depth[0]*1000.),
	          soil_con.bubble,soil_con.expt[0]);
            if(ice0<0.) ice0=0.;
            if(ice0>soil_con.max_moist[0]/(soil_con.depth[0]*1000.)) 
	      ice0=soil_con.max_moist[0]/(soil_con.depth[0]*1000.);
          }
          else ice0=0.;
        }
        else ice0=0.;

        if(options.FULL_ENERGY && !options.FROZEN_SOIL) {
          /** Compute Soil Thermal Properties if not using Frozen Soils **/
          null_ptr = NULL;
          soil_thermal_calc(soil_con, cell[iveg].layer, energy[iveg], null_ptr,
              null_ptr, null_ptr, null_ptr, null_ptr, options.Nlayer,
              gp.Ulayer+gp.Llayer+2);
 
          /** Save Thermal Conductivities for Energy Balance **/
          energy[iveg].kappa[0] = cell[iveg].layer[0].kappa;
          energy[iveg].Cs[0] = cell[iveg].layer[0].Cs;
          energy[iveg].kappa[1] = cell[iveg].layer[1].kappa;
          energy[iveg].Cs[1] = cell[iveg].layer[1].Cs;
        }

      }
  
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        snow[iveg].vapor_flux = 0.;
        snow[iveg].canopy_vapor_flux = 0.;
      }
      veg_var[iveg].throughfall = 0.;

      /** Compute latent heat of vaporization (J/kg) **/
      cell[iveg].Le = (2.501 - 0.002361*atmos->air_temp) * 1.0e6;
   
      /** Compute the aerodynamic resistance **/
      if(options.FULL_ENERGY && options.SNOW_MODEL) {
        snow_depth = snow[iveg].depth;
        snow_present = snow[iveg].snow;
      }
      if(iveg<Nveg) {
        displacement = veg_lib[veg_class].displacement[dmy[rec].month-1];
        roughness = veg_lib[veg_class].roughness[dmy[rec].month-1];
        overstory = veg_lib[veg_class].overstory;
      }
      if(iveg==Nveg || roughness==0) {
        displacement = 0.;
        roughness = soil_con.rough;
        overstory = FALSE;
      }
      tmp_wind[0] = atmos->wind;
 
      /** for PILPS comparison only **/
      height = calc_veg_height(displacement);
      if(!overstory) {
	if(displacement<gp.wind_h) ref_height = gp.wind_h;
	else ref_height = displacement + gp.wind_h + roughness;
      }
      else {
	ref_height = veg_lib[veg_class].wind_h;
	if(ref_height<=0.) ref_height = calc_veg_height(displacement)
			             + gp.wind_h;
      }
      CalcAerodynamic(overstory,iveg,Nveg,OVERSTORY_ATTENUATION,
                      height,soil_con.rough,soil_con.snow_rough,
                      &displacement,&roughness,&ref_height,TRUNK_RATIO,
                      tmp_wind,cell[iveg].aero_resist);
  
      /** Check for Presence of Snow Pack, or New Snow **/
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
  
        /**************************************************
           Full Energy Balance Model
        **************************************************/
  
        rainfall = atmos->rainonly;
        snowfall = atmos->prec - atmos->rainonly;

        if(options.SNOW_MODEL && (snow[iveg].snow || snowfall > 0.
            || (snow[iveg].snow_canopy>0. && overstory))) {

          /************************************************
             New Snow, or Old Snow on Ground or in Canopy 
          ************************************************/

          SNOW=TRUE;

          /***** Compute Radiation Over Snow Pack *****/
          if(snow[iveg].snow || snowfall > 0.)
            snow[iveg].albedo = snow_albedo(atmos->prec-atmos->rainonly,
                                snow[iveg].swq,snow[iveg].coldcontent,gp.dt,
                                snow[iveg].last_snow);
          else
            snow[iveg].albedo = snow_albedo(atmos->prec-atmos->rainonly,
                                snow[iveg].snow_canopy,atmos->air_temp,gp.dt,
                                snow[iveg].last_snow);
          if(atmos->prec-atmos->rainonly==0)
            snow[iveg].last_snow++;
          else snow[iveg].last_snow=1;
          out_short = snow[iveg].albedo * atmos->shortwave;
          atmos->net_short = (1.0 - snow[iveg].albedo) * atmos->shortwave;
          atmos->rad = atmos->net_short + atmos->longwave 
                     - STEFAN_B * pow(snow[iveg].surf_temp+KELVIN,4.0);

          if(iveg!=Nveg) {
            /***** Check if Vegetation Covered by Snow *****/
            if(snow[iveg].depth>height) surf_atten=1.;

            /***** Compute Canopy Interception, if Overstory Present *****/
            if(overstory && (snowfall>0. || snow[iveg].snow_canopy>0.)) {
              if(snow[iveg].snow) surf_temp = snow[iveg].surf_temp;
              else surf_temp = energy[iveg].T[0];
              snow_intercept((double)gp.dt,1.,
			     veg_lib[veg_class].LAI[dmy[rec].month-1],
			     veg_lib[veg_class].Wdmax[dmy[rec].month-1],
			     cell[iveg].aero_resist[1],atmos->density,
			     atmos->vp,cell[iveg].Le,atmos->shortwave,
			     atmos->longwave,atmos->pressure,atmos->air_temp,
			     atmos->vpd,tmp_wind[1],&rainfall,&snowfall,
			     &veg_var[iveg].Wdew,&snow[iveg].snow_canopy,
			     &snow[iveg].tmp_int_storage,
                             &snow[iveg].canopy_vapor_flux,&canopy_temp,
                             &energy[iveg].refreeze_energy);
              atmos->longwave = STEFAN_B * pow(canopy_temp+KELVIN,4.0);
              veg_var[iveg].throughfall = rainfall + snowfall;
            }
            else if(overstory) {
              Evap = canopy_evap(cell[iveg].layer,&veg_var[iveg],
				 FALSE,veg_class,
				 dmy[rec].month,veg_var[iveg].Wdew,
				 atmos->air_temp,(double)gp.dt,
				 atmos->rad,atmos->vpd,atmos->net_short,
				 atmos->air_temp,
				 cell[iveg].aero_resist[0],rainfall,
				 displacement,roughness,ref_height,
				 soil_con.elevation,soil_con.depth,
				 soil_con.Wcr,soil_con.Wpwp);
              rainfall = veg_var[iveg].throughfall;
              atmos->longwave = STEFAN_B * pow(atmos->air_temp+KELVIN,4.0);
            }
            else if(snowfall > 0.) {
              tmp_rain = calc_rainonly(atmos->air_temp,veg_var[iveg].Wdew,
				       gp.MAX_SNOW_TEMP,gp.MIN_RAIN_TEMP);
              rainfall += tmp_rain;
              snowfall += veg_var[iveg].Wdew - tmp_rain;
              tmp_rain = veg_var[iveg].Wdew;
              veg_var[iveg].Wdew = 0.;
              if(debug.PRT_BALANCE) debug.store_moist[0] = 0.;
            }
          }

          tmp_swq = snow[iveg].swq;
          if(tmp_swq>0.) tmp_snow_surf_temp = snow[iveg].surf_temp;
          else tmp_snow_surf_temp = 0.;
    
          /** Snowpack on Ground **/
          if(snow[iveg].swq>0.0 || snowfall > 0) {
            snow_melt(atmos,soil_con,gp.wind_h+soil_con.snow_rough,
                      cell[iveg].aero_resist[2],
                      cell[iveg].Le,&snow[iveg],(double)gp.dt,0.00,
                      soil_con.snow_rough,surf_atten,rainfall,snowfall,
		      tmp_wind[2],energy[iveg].T[0],
                      &energy[iveg].advection,&energy[iveg].deltaCC,
                      &energy[iveg].grnd_flux,&energy[iveg].latent,
                      &energy[iveg].sensible,&snow[iveg].Qnet,
		      &energy[iveg].Trad[0],&energy[iveg].refreeze_energy);
            ppt=atmos->melt;
            energy[iveg].albedo = snow[iveg].albedo;
	    tmp_snow_energy = energy[iveg].advection - energy[iveg].deltaCC
	                    + energy[iveg].refreeze_energy;

            /** Compute Snow Parameters **/
            if(snow[iveg].swq > 0.) {
              snow[iveg].snow = TRUE;
	      snow[iveg].density = snow_density(dmy[rec].day_in_year,
                                   snowfall,atmos->air_temp,
                                   tmp_swq,snow[iveg].depth,
                                   snow[iveg].coldcontent,(double)gp.dt);
              snow[iveg].depth = 1000. * snow[iveg].swq 
                               / snow[iveg].density;  /* HBH 7.2.1 */
            }
            else {
              snow[iveg].snow = FALSE;
              snow[iveg].density = 0.;
              snow[iveg].depth = 0.;
	      snow[iveg].surf_water=0;
	      snow[iveg].pack_water=0;
	      snow[iveg].surf_temp=0;
	      snow[iveg].pack_temp=0;
            }
  
          }
          else {
            ppt = rainfall;
            tmp_snow_energy = 0.;
          }
        }
        else {

          /******************
            No Snow Present
          ******************/

          SNOW=FALSE;
          ppt = rainfall;
          if(!options.SNOW_MODEL) {
            ppt += snowfall;
            rainfall += snowfall;
            snowfall=0.;
          }
	  tmp_snow_energy = 0.;
        }
  
        /** Compute Surface Temperature, and Update Snowpack Variables **/
        if(options.SNOW_MODEL && snow[iveg].swq>0.0) {
    
          /************************************************
            Energy Balance Case 1: Ground Snow is Present
          ************************************************/
         
	  for(i=0;i<options.Nlayer;i++) tmp_layer[i] = cell[iveg].layer[i];
    
          if(snow[iveg].depth >= SNOWTHRES) {
  
            /***********************************************
              Ground Snow Present Case 1a: Thick Snow Pack
            ***********************************************/

            if(debug.DEBUG) fprintf(stdout,"\nFIRST SNOW ITER- %i/%i/%i at %i[%i]: air temp = %lf\n",
                dmy[rec].month,dmy[rec].day,dmy[rec].year,dmy[rec].hour,
                rec,atmos->air_temp); 
  
            /** Solve Surface Temperature Through Snow Pack **/
            /** replace displacement with veg height later **/
  
            if(options.FROZEN_SOIL || options.CALC_SNOW_FLUX) {
              surf_temp = calc_snow_ground_flux(dp,moist,ice0,
						energy[iveg].Trad[0],&last_T1,
						&energy[iveg],&snow[iveg],
						soil_con,gp);

              if(options.FROZEN_SOIL) {
	        tmp_surf_temp=surf_temp;
	        frozen_soil(soil_con,cell[iveg].layer,&energy[iveg],rec,
	                    iveg,gp.dt,tmp_surf_temp,gp.Ulayer,gp.Llayer);
              }
	      else {
		energy[iveg].T[0] = surf_temp;
		energy[iveg].T[1] = last_T1;
	      }
            }
            else {
              surf_temp=energy[iveg].T[0];
              last_T1=energy[iveg].T[1];

              energy[iveg].T[0] = surf_temp;
              energy[iveg].T[1] = last_T1;
/*****
              surf_temp = snow[iveg].surf_temp;
              energy[iveg].T[0] = surf_temp;
	      C1 = energy[iveg].Cs[1] * dp / soil_con.depth[1] 
	          * ( 1. - exp(-soil_con.depth[1]/dp));
	      C2 = - ( 1. - exp(soil_con.depth[0]/dp) ) 
	          * exp(-soil_con.depth[1]/dp);
	      C3 = energy[iveg].kappa[0]/soil_con.depth[0] 
	          - energy[iveg].kappa[1]/soil_con.depth[0] 
	          + energy[iveg].kappa[1]/soil_con.depth[0]
	          * exp(-soil_con.depth[0]/dp);
	      energy[iveg].T[1] = (energy[iveg].kappa[0]/2./soil_con.depth[0]
                   / soil_con.depth[1] * surf_temp 
                  + C1/((double)(gp.dt)*3600.)*energy[iveg].T[1] 
                  + (2.*C2-1.+exp(-soil_con.depth[0]/dp))
                  * energy[iveg].kappa[1]/2./soil_con.depth[0]/soil_con.depth[1]
                  * energy[iveg].T[gp.Ulayer+gp.Llayer+1])
                  / (C1/((double)(gp.dt)*3600.) 
                  + energy[iveg].kappa[1]/soil_con.depth[0]
                  / soil_con.depth[1]*C2 + C3/2./soil_con.depth[1]);
*****/
            }
  
          } /** End Thick Snowpack **/

          else { /** Snow not deeper than threshold, use bare soil balance **/
    
            /***********************************************
              Ground Snow Present Case 1b: Thin Snow Pack
            ***********************************************/

            /** assume veg height exceeds thin soil pack depth **/
            surf_temp = snow[iveg].surf_temp;
    
            /**fprintf(stdout,"\nTHIN SNOW- %i/%i/%i at %i[%i]: air temp = %lf\n",
               dmy[rec].month,dmy[rec].day,dmy[rec].year,dmy[rec].hour,
               rec,atmos->air_temp);**/
    
	    if(options.FROZEN_SOIL) {
	      tmp_surf_temp=surf_temp;
	      frozen_soil(soil_con,cell[iveg].layer,&energy[iveg],rec,
	                  iveg,gp.dt,tmp_surf_temp, gp.Ulayer, gp.Llayer);
	    }
	    else {
              energy[iveg].T[0] = surf_temp;
	      C1 = energy[iveg].Cs[1] * dp / soil_con.depth[1] 
	          * ( 1. - exp(-soil_con.depth[1]/dp));
	      C2 = - ( 1. - exp(soil_con.depth[0]/dp) ) 
	          * exp(-soil_con.depth[1]/dp);
	      C3 = energy[iveg].kappa[0]/soil_con.depth[0] 
	          - energy[iveg].kappa[1]/soil_con.depth[0] 
	          + energy[iveg].kappa[1]/soil_con.depth[0]
	          * exp(-soil_con.depth[0]/dp);
	      energy[iveg].T[1] = (energy[iveg].kappa[0]/2./soil_con.depth[0]
                   / soil_con.depth[1] * surf_temp 
                  + C1/((double)(gp.dt)*3600.)*energy[iveg].T[1] 
                  + (2.*C2-1.+exp(-soil_con.depth[0]/dp))
                  * energy[iveg].kappa[1]/2./soil_con.depth[0]/soil_con.depth[1]
                  * energy[iveg].T[gp.Ulayer+gp.Llayer+1])
                  / (C1/((double)(gp.dt)*3600.) 
                  + energy[iveg].kappa[1]/soil_con.depth[0]
                  / soil_con.depth[1]*C2 + C3/2./soil_con.depth[1]);
	    }
  
	  }  /** End thin snowpack **/
  
        }  /** End Ground Snow Present **/
      
        else {
    
          /************************************************
            Energy Balance Case 2: No Ground Snow Present
          ************************************************/
         
	  if(debug.DEBUG) fprintf(stdout,"\nSNOW MELTED- %i/%i/%i at %i[%i]: air temp = %lf\n",
              dmy[rec].month,dmy[rec].day,dmy[rec].year,dmy[rec].hour,
              rec,atmos->air_temp);
  
          if(snow[iveg].snow_canopy > 0.) {

            /*************************************
              Canopy Has Snow, but No Snow on Ground
            *************************************/

            if(atmos->prec-atmos->rainonly==0)
              snow[iveg].last_snow++;
            else snow[iveg].last_snow=1;

          }
          else if(snow[iveg].snow) snow[iveg].snow = FALSE;
  
          /***** Compute Surface Energy Balance *****/

          for(i=0;i<options.Nlayer;i++) tmp_layer[i] = cell[iveg].layer[i];
          if(iveg<Nveg) {
            tmp_veg_var = veg_var[iveg];
          }
          else {
            tmp_veg_var = veg_var[0];
          }
          roughness = soil_con.rough;
          T0 = calc_surf_energy_bal(rec,iveg,options.Nlayer,Nveg,SNOW,
				    ice0,moist,dp,
				    surf_atten,energy[iveg].T[0],&last_T1,
				    cell[iveg].aero_resist,
				    cell[iveg].Le,Ls,mu,&x,0.,roughness,
				    ref_height+roughness,
				    tmp_snow_energy,
				    tmp_wind,snow[iveg].canopy_vapor_flux
				    +snow[iveg].vapor_flux,
				    rainfall,atmos,&tmp_veg_var,
				    &energy[iveg],gp,tmp_layer,soil_con,
				    veg_class,dmy[rec]);  

          energy[iveg].T[0] = T0;
          for(i=0;i<options.Nlayer;i++) cell[iveg].layer[i]=tmp_layer[i];
          if(iveg<Nveg) veg_var[iveg] = tmp_veg_var;

          if(!SNOW && iveg!=Nveg && veg_lib[veg_class].LAI[dmy[rec].month-1])
            ppt=veg_var[iveg].throughfall;
          else if(!SNOW) ppt=rainfall;
	  
          if(options.FROZEN_SOIL) {
            tmp_surf_temp=energy[iveg].T[0];
            frozen_soil(soil_con,cell[iveg].layer,&energy[iveg],rec ,
		        iveg,gp.dt,tmp_surf_temp,gp.Ulayer,gp.Llayer);
          }
          else if(options.FULL_ENERGY) {
            energy[iveg].T[1] = last_T1;
          }
	  /*****
	    if(debug.PRT_BALANCE && debug.DEBUG) {
	    printf("\nAfter Frozen Soil\n");
	    write_layer(cell[iveg].layer,iveg,options.Nlayer,soil_con.depth);
	    }
	    *****/
        } /** End No Ground Snow **/
 
        /** Compute Arno Runoff *****/

        cell[iveg].inflow = ppt;
        runoff(cell[iveg].layer,&energy[iveg],soil_con,&cell[iveg].runoff,
               &cell[iveg].baseflow,mu,x,ppt,gp.dt,gp.Ulayer+gp.Llayer+2);
  
        out_short = atmos->albedo * atmos->shortwave;
        atmos->net_short = (1.0 - atmos->albedo) * atmos->shortwave;
  
      } /** End Full Energy Balance Model **/

      else {

        /******************************************************
	  Water Balance Model
        ******************************************************/

        if(options.SNOW_MODEL) {

          ppt = 0.;
          tmp_melt = 0.;
          throughfall = 0.;
          tmp_throughfall = 0.;
          tmp_rad = atmos->rad;
          tmp_canopy_vapor_flux = 0.;
          tmp_vapor_flux = 0.;
          tmp_canopyevap = 0.;
          energy[iveg].refreeze_energy=0.;
          tmp_layerevap = (double*)calloc(options.Nlayer,sizeof(double));

          if(snow[iveg].snow || (atmos->prec-atmos->rainonly) > 0.
            || (snow[iveg].snow_canopy>0. && overstory)) {

	    /*********************************************
	      New Snow, or Old Snow on Ground or in Canopy
              *********************************************/

            SNOW = TRUE;

            if(iveg<Nveg && !overstory && veg_var[iveg].Wdew>0.) {
	      /** If canopy not present, empty vegetation storage since
		vegetation is assumed to be covered by ground snow **/
              tmp_rain = calc_rainonly(atmos->air_temp,veg_var[iveg].Wdew,
				       gp.MAX_SNOW_TEMP,gp.MIN_RAIN_TEMP);
              rainfall += tmp_rain;
              snowfall += veg_var[iveg].Wdew - tmp_rain;
              tmp_rain = veg_var[iveg].Wdew;
              veg_var[iveg].Wdew = 0.;
              if(debug.PRT_BALANCE) debug.store_moist[0] = 0.;
            }
            else tmp_rain = 0.;

            for(hour=0;hour<4;hour++) {

	      /** Snow melt model is solved using a 6 hour time step,
	        Initialize parameters for the shorter time step **/

              atmos->air_temp = calc_air_temperature(tmax,tmin,hour*6);
              atmos->rainonly = calc_rainonly(atmos->air_temp,atmos->prec,
					      gp.MAX_SNOW_TEMP,
					      gp.MIN_RAIN_TEMP);
              atmos->density = 3.486*atmos->pressure/(275.0
                             + atmos->air_temp);
              rainfall = (atmos->rainonly + tmp_rain)/4.;
              snowfall = (atmos->prec - atmos->rainonly)/4.;
              snow[iveg].vapor_flux = 0.;
              snow[iveg].canopy_vapor_flux = 0.;
	      for(i=0;i<options.Nlayer;i++)
		cell[iveg].layer[i].evap += 0.;
              if(snow[iveg].snow || snowfall > 0.)
                snow[iveg].albedo = snow_albedo(snowfall,
                                    snow[iveg].swq,snow[iveg].coldcontent,
                                    gp.dt/4,snow[iveg].last_snow);
              else
                snow[iveg].albedo = snow_albedo(snowfall,
                                    snow[iveg].snow_canopy,atmos->air_temp,
                                    gp.dt/4,snow[iveg].last_snow);
              calc_long_shortwave(&inshort,
                   &inlong,&atmos->tskc,atmos->air_temp,
                   atmos->vp,(double)soil_con.time_zone_lng,
                   (double)soil_con.lng,(double)soil_con.lat,
                   (double)dmy[rec].day_in_year,(double)hour*6.,
                   FALSE,FALSE,TRUE);
              if((overstory && snow[iveg].snow_canopy>0.)
                  || !overstory || snowfall>0.) {
                energy[iveg].albedo = snow[iveg].albedo;
              }
              else {
                energy[iveg].albedo = veg_lib[veg_class].albedo[dmy[rec].month-1];
              }
              out_short = energy[iveg].albedo * inshort;
              atmos->net_short = inshort * (1.0 - energy[iveg].albedo);
              atmos->longwave = inlong;
	      atmos->rad = atmos->net_short + inlong
                         - STEFAN_B * pow(energy[iveg].Trad[0]+KELVIN,4.0);

              if(iveg!=Nveg) {
                /***** Check if Vegetation Covered by Snow *****/
                if(snow[iveg].depth>height) surf_atten=1.;

                /***** Compute Canopy Interception, if Overstory Present *****/
                if(overstory
                    && (snowfall>0. || snow[iveg].snow_canopy>0.)) {
		  /** Canopy intercepts or stores snow **/
                  if(snow[iveg].snow) surf_temp = snow[iveg].surf_temp;
                  else surf_temp = atmos->air_temp;
                  snow_intercept((double)gp.dt,1.,
				 veg_lib[veg_class].LAI[dmy[rec].month-1],
				 veg_lib[veg_class].Wdmax[dmy[rec].month-1],
				 cell[iveg].aero_resist[1],atmos->density,
				 atmos->vp,cell[iveg].Le,inshort,
				 atmos->longwave,atmos->pressure,
				 atmos->air_temp,
				 atmos->vpd,tmp_wind[0],&rainfall,&snowfall,
				 &veg_var[iveg].Wdew,&snow[iveg].snow_canopy,
				 &snow[iveg].tmp_int_storage,
				 &snow[iveg].canopy_vapor_flux,&canopy_temp,
				 &energy[iveg].refreeze_energy);
                  atmos->longwave = STEFAN_B * pow(canopy_temp+KELVIN,4.0);
                  tmp_throughfall = rainfall + snowfall;
                }
                else if(overstory) {
		  /** Canopy intercept rain only, and has no snow storage **/
                  Evap = canopy_evap(cell[iveg].layer,&veg_var[iveg],
				     TRUE,veg_class,dmy[rec].month,
				     veg_var[iveg].Wdew,atmos->air_temp,
				     (double)gp.dt,
				     atmos->rad,atmos->vpd,atmos->net_short,
				     atmos->air_temp,
				     cell[iveg].aero_resist[0],rainfall,
				     displacement,roughness,ref_height,
				     soil_con.elevation,soil_con.depth,
				     soil_con.Wcr,soil_con.Wpwp);
                  atmos->longwave = (1. - surf_atten) * STEFAN_B
                                  * pow(atmos->air_temp+KELVIN,4.0)
                                  + surf_atten * atmos->longwave;
                  out_short = snow[iveg].albedo * inshort;
                  atmos->net_short = inshort * (1.0 - snow[iveg].albedo);
                  rainfall = veg_var[iveg].throughfall;
                  tmp_throughfall = veg_var[iveg].throughfall;
                  tmp_canopyevap += veg_var[iveg].canopyevap;
		  for(i=0;i<options.Nlayer;i++)
		    tmp_layerevap[i] += cell[iveg].layer[i].evap;
                }
              }

              tmp_swq = snow[iveg].swq;
              throughfall += tmp_throughfall;

              if(snow[iveg].swq>0.0 || snowfall > 0) {
		tmp_snowinflow += rainfall + snowfall;
                snow_melt(atmos,soil_con,gp.wind_h+soil_con.snow_rough,
                          cell[iveg].aero_resist[2],
                          cell[iveg].Le,&snow[iveg],(double)gp.dt/4,0.00,
                          soil_con.snow_rough,surf_atten,rainfall,
			  snowfall,tmp_wind[2],0.,&dummy,&dummy,
                          &dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
                tmp_melt+=atmos->melt;

                if(snow[iveg].swq>0.) {
		  /** Recompute snowpack parameters, if snowpack still
		    present **/
                  snow[iveg].snow = TRUE;
                  snow[iveg].density = snow_density(dmy[rec].day_in_year,
						    snowfall,atmos->air_temp,
						    tmp_swq,snow[iveg].depth,
						    snow[iveg].coldcontent,
						    (double)gp.dt/4.);
  
                  snow[iveg].depth = 1000. * snow[iveg].swq
                                   / snow[iveg].density;  /* HBH 7.2.1 */
                }
                else {
		  /** Reset snowpack parameters, when snowpack has melted **/
                  snow[iveg].snow = FALSE;
                  snow[iveg].density = 0.;
                  snow[iveg].depth = 0.;
                }
                if(atmos->prec-atmos->rainonly==0)
                  snow[iveg].last_snow++;
                else snow[iveg].last_snow=1;

              }
              else {

                /****************
                  No Snow Present
		  ****************/

                if(overstory) ppt += tmp_throughfall;
                else ppt += rainfall+snowfall;
                if(overstory && snow[iveg].snow_canopy>0.) {
                  snow[iveg].snow = FALSE;
                  snow[iveg].depth=0.;
                  snow[iveg].density=0.;
                  snow[iveg].surf_water=0.;
                  snow[iveg].pack_water=0.;
                  snow[iveg].surf_temp=0.;
                  snow[iveg].pack_temp=0.;
                }
                else {
                  snow[iveg].snow = FALSE;
                  snow[iveg].last_snow=0;
                  snow[iveg].depth=0.;
                  snow[iveg].density=0.;
                  snow[iveg].surf_water=0.;
                  snow[iveg].pack_water=0.;
                  snow[iveg].surf_temp=0.;
                  snow[iveg].pack_temp=0.;
                  snow[iveg].snow_canopy=0.;
                  snow[iveg].tmp_int_storage=0.;
                }
              }
              tmp_canopy_vapor_flux += snow[iveg].canopy_vapor_flux;
              tmp_vapor_flux += snow[iveg].vapor_flux;
            }
            atmos->longwave = tmp_longwave;
            atmos->rad = tmp_rad;
            atmos->melt = tmp_melt;
            snow[iveg].canopy_vapor_flux = tmp_canopy_vapor_flux;
            snow[iveg].vapor_flux = tmp_vapor_flux;
            veg_var[iveg].canopyevap = tmp_canopyevap;
            veg_var[iveg].throughfall = throughfall;
            ppt += tmp_melt;
	    for(i=0;i<options.Nlayer;i++)
	      cell[iveg].layer[i].evap += tmp_layerevap[i];

          }
          else {
            SNOW=FALSE;
            rainfall=atmos->prec;
          }

        }
        else {
          SNOW = FALSE;
          rainfall = atmos->prec - atmos->melt;
          if(rainfall<0.) rainfall=0.;
        }
 
        if(!SNOW) {
          if(iveg!=Nveg) {
            /** Compute Evaporation from Vegetation **/
            if(veg_lib[veg_class].LAI[dmy[rec].month-1] > 0.0) {
              Evap = canopy_evap(cell[iveg].layer,&veg_var[iveg],
				 TRUE,veg_class,dmy[rec].month,
				 veg_var[iveg].Wdew,atmos->air_temp,
				 (double)gp.dt,
				 atmos->rad,atmos->vpd,atmos->net_short,
				 atmos->air_temp,
				 cell[iveg].aero_resist[0],rainfall,
				 displacement,roughness,ref_height,
				 soil_con.elevation,soil_con.depth,
				 soil_con.Wcr,soil_con.Wpwp);
              ppt = veg_var[iveg].throughfall;
            }
            else {
              Evap = arno_evap(cell[iveg].layer, atmos->rad,
			       atmos->air_temp, atmos->vpd, atmos->net_short,
			       soil_con.depth[0], soil_con.max_moist[0],
			       soil_con.elevation, soil_con.b_infilt, 
			       soil_con.max_infil, atmos->air_temp,
			       displacement, roughness, ref_height, 
			       cell[iveg].aero_resist[0], (double)gp.dt);
              ppt = rainfall;
            }
          }

          /** Compute Evaporation from Bare Soil **/
          else Evap = arno_evap(cell[iveg].layer, atmos->rad,
				atmos->air_temp, atmos->vpd, 
				atmos->net_short, soil_con.depth[0],
				soil_con.max_moist[0], soil_con.elevation,
				soil_con.b_infilt, soil_con.max_infil, 
				atmos->air_temp, displacement, roughness, 
				ref_height, cell[iveg].aero_resist[0], 
				(double)gp.dt);

        }
        if(!options.SNOW_MODEL) ppt += atmos->melt;

        /** Compute Arno Runoff **/
        cell[iveg].inflow = ppt;

        runoff(cell[iveg].layer, &energy[iveg],soil_con, &cell[iveg].runoff,
               &cell[iveg].baseflow,mu,x,ppt,gp.dt,gp.Ulayer+gp.Llayer+2);
 

      } /** End Water Balance Model **/
  
      /** Store Energy Balance Data **/
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
        energy[iveg].shortwave = atmos->shortwave;
        energy[iveg].longwave = atmos->longwave;
	energy[iveg].longwave -= STEFAN_B 
	                       * pow(energy[iveg].Trad[0]+KELVIN,4.0);
      }

      /***** Reset modified Paramters *****/
      atmos->longwave = tmp_longwave;

      /***** Debug Option Output *****/

      if(iveg<Nveg) tmp_veg = &veg_var[iveg];
      else tmp_veg = NULL;
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        tmp_energy = &energy[iveg];
        tmp_snow = &snow[iveg];
      }
      if(debug.PRT_BALANCE) {
        debug.inflow[options.Nlayer+2] = atmos->prec;
        if(iveg<Nveg && (!SNOW || overstory)) {
          debug.inflow[0] = atmos->prec + tmp_rain;
          debug.inflow[1] = veg_var[iveg].throughfall;
          debug.outflow[0] = veg_var[iveg].throughfall;
        }
        else {
          debug.inflow[0] = 0.;
          debug.inflow[1] = atmos->prec + tmp_rain;
          debug.outflow[0] = 0.;
        }
        if((options.FULL_ENERGY || options.SNOW_MODEL) 
            && (snow[iveg].snow || debug.store_moist[1]>0 || atmos->melt>0.)) {
	  if(!options.FULL_ENERGY) debug.inflow[1] = tmp_snowinflow;
          debug.outflow[1] = atmos->melt;
        }
        else if(iveg<Nveg) {
          debug.inflow[1] = 0.;
          debug.inflow[2] = veg_var[iveg].throughfall - cell[iveg].runoff;
          debug.outflow[1] = 0.;
        }
        else {
          debug.inflow[1] = 0.;
          debug.inflow[2] = atmos->prec + tmp_rain - cell[iveg].runoff;
          debug.outflow[1] = 0.;
        }
	if(!options.FULL_ENERGY && !options.SNOW_MODEL) {
	  debug.inflow[2] += atmos->melt;
        }
        if(veg_var[iveg].throughfall>0. && atmos->prec==0.
            && debug.store_moist[0]==0.) {
          debug.inflow[options.Nlayer+2] = veg_var[iveg].throughfall;
          debug.inflow[0] = veg_var[iveg].throughfall;
        }
          
      }
      if(iveg!=Nveg) {
        write_debug(atmos[0],soil_con,cell[iveg],
            tmp_energy,tmp_snow,tmp_veg,dmy[rec],gp,
            out_short,mu,Nveg,
            iveg,rec,gridcell,prcpdist,NEWCELL); 
      }
      else {
        write_debug(atmos[0],soil_con,cell[iveg],
            tmp_energy,tmp_snow,tmp_veg,dmy[rec],gp,
            out_short,mu,Nveg,iveg,rec,gridcell,
            prcpdist,NEWCELL); 
      }
  
    }

    atmos[0] = tmp_atmos;

  } /** end of vegetation loop **/

  /** END Temperature and Moisture Profile Debugging Output **/
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    free ((char *)tmp_layer);
  }

}

void store_max_min_temp(atmos_data_struct *atmos,
                        double *tmax,
                        double *tmin,
                        int rec,
                        int Nrecs,
                        char prcpdist) {
/**********************************************************************
  This subroutine sets prepares arrays of maximum and minimum 
  daily air temperature.
**********************************************************************/

  static int last_rec;

  if(rec==0 && prcpdist==0) {
    tmax[0] = tmax[1] = atmos[0].tmax;
    tmin[0] = tmin[1] = atmos[0].tmin;
    tmax[2] = atmos[1].tmax;
    tmin[2] = atmos[1].tmin;
    last_rec = rec;
  }
  else if(rec!=last_rec) {
    last_rec = rec;
    if(rec==Nrecs-1) {
      tmax[0] = tmax[1];
      tmax[1] = tmax[2];
      tmin[0] = tmin[1];
      tmin[1] = tmin[2];
    }
    else {
      tmax[0] = tmax[1];
      tmax[1] = tmax[2];
      tmax[2] = atmos[1].tmax;
      tmin[0] = tmin[1];
      tmin[1] = tmin[2];
      tmin[2] = atmos[1].tmin;
    }
  }
}
