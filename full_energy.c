#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

#define OVERSTORY_ATTENUATION 0.5  /* attenutation coefficient for overstory */
#define TRUNK_RATIO           0.2  /* Fraction of Tree Height that is trunk */
#define SNOW_STEP             1    /* Time step in hours to use when solving
				      snow model in water balance mode */

void full_energy(int rec,
                 atmos_data_struct   *atmos,
                 soil_con_struct      soil_con,
                 veg_con_struct      *veg_con,
                 dist_prcp_struct    *prcp,
                 dmy_struct          *dmy,
                 global_param_struct  gp,
                 int                  gridcell,
                 char                 NEWCELL,
		 double              *prec,
		 double              *rainonly,
		 double              *melt)
/**********************************************************************
	full_energy	Keith Cherkauer		January 8, 1997

  This subroutine controls the model core, it solves both the energy
  and water balance models, as well as frozen soils.  

  modifications:
  07-98 resturctured to fix problems with distributed precipitation, 
        and to add the ability to solve the snow model at different 
	elevation bands within a single grid cell.                 KAC

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
  extern debug_struct    debug;

  static double          tmax[3];
  static double          tmin[3];
  static int             tmax_hour[3];
  static int             tmin_hour[3];

  char                   RUN_TYPE;
  char                   overstory;
  char                   ANY_SNOW;
  char                   SNOW;
  char                   SOLVE_SURF_ENERGY;
  int                    i, j;
  int                    Ndist;
  int                    dist;
  int                    iveg;
  int                    sveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands;
  int                    hour;
  double                *ppt;
  double                 tmp_surf_temp;
  double                 last_T1;
  double                 out_short;
  double                 inshort;
  double                 inlong;
  double                 dp;
  double                 ice0, moist;
  double                 surf_atten;
  double                 Tsurf;
  double                 Tgrnd;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 Ttemp;
  double                 rainfall[2]; 
/*   double                 snowfall[2]; */
  double                 height;
  double                 displacement;
  double                 roughness;
  double                 ref_height;
  double                 Le;
  double                 Ls;
  double                 Evap;
  double                *Melt;
  double                 air_temp[24];
  double                 bare_albedo;
  double                *snow_inflow;
  double                 step_air_temp;
  double                 step_vp;
  double                 step_vpd;
  double                 step_rad;
  double                 step_density;
  double                 step_net_short;
  double                 tmp_aero_resist;
  double                *tmp_throughfall[2];
  double                 tmp_rad;
  double                 tmp_wind[3];
  double                *tmp_melt;
  double                *tmp_ppt;
  double                *tmp_vapor_flux;
  double                *tmp_canopy_vapor_flux;
  double                *tmp_canopyevap[2];
  double                 tmp_snow_energy;
  double                 tmp_vp;
  double                 tmp_vpd;
  double                 tmp_Wdew[2];
  double                 tmp_mu;
  double                *tmp_layerevap[2];
  atmos_data_struct      tmp_atmos;
  layer_data_struct     *tmp_layer[2];
  veg_var_struct         tmp_veg_var[2];
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;
  energy_bal_struct    **energy;
  energy_bal_struct     *tmp_energy;
  energy_bal_struct     *ptr_energy;
  energy_bal_struct     *band_energy;
  snow_data_struct     **snow;
  snow_data_struct      *tmp_snow;
  veg_var_struct        *tmp_veg[2];

  cell    = prcp[0].cell;
  veg_var = prcp[0].veg_var;
  snow    = prcp[0].snow;
  energy  = prcp[0].energy;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = options.SNOW_BAND;

  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    for(i=0;i<2;i++) {
      tmp_layer[i] = (layer_data_struct *)calloc((options.Nlayer),
						 sizeof(layer_data_struct));
    }
    tmp_energy = (energy_bal_struct *)calloc(1,sizeof(energy_bal_struct));
  }
  if(options.SNOW_MODEL || options.FULL_ENERGY) {
    snow_inflow = (double *)calloc(Nbands,sizeof(double));
  }
  band_energy       = (energy_bal_struct *)calloc(1,
						  sizeof(energy_bal_struct));
  ppt               = (double *)calloc(2*Nbands,sizeof(double));
  Melt              = (double *)calloc(2*Nbands,sizeof(double));
  if(gp.dt == 24 && options.SNOW_MODEL) {
    store_max_min_temp(atmos,tmax,tmax_hour,tmin,tmin_hour,rec,
		       gp.nrecs);
  }
      
  Nveg         = veg_con[0].vegetat_type_num;
/*****  tmp_longwave = atmos->longwave;*****/
  tmp_atmos    = atmos[0];

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

      /**************************************************
        Initialize Model Parameters
      **************************************************/
      for(dist=0;dist<Ndist;dist++) {
	for(band=0;band<Nbands;band++) {
	  if(soil_con.AreaFract[band]>0) {
	    for(i=0;i<options.Nlayer;i++) 
	      cell[dist][iveg][band].layer[i].evap = 0.;
	    if(iveg<Nveg) {
	      veg_var[dist][iveg][band].canopyevap  = 0.;
	      veg_var[dist][iveg][band].throughfall = 0.;
	    }
	  }
	}
      }

      if(iveg<Nveg) veg_class = veg_con[iveg].veg_class;
      else veg_class = 0;

      /** Compute Surface Attenuation due to Vegetation Coverage **/
      if(iveg<Nveg)
	surf_atten=exp(-0.50*veg_lib[veg_class].LAI[dmy[rec].month-1]);
      else surf_atten = 1.;
        
      /**************************************************
        Store Water Balance Terms for Debugging
      **************************************************/
      if(debug.DEBUG || debug.PRT_MOIST || debug.PRT_BALANCE) {
        /** Compute current total moisture for water balance check **/
	store_moisture_for_debug(iveg,Nveg,prcp[0].mu,cell,
				 veg_var,snow,soil_con);
	if(debug.PRT_BALANCE) {
	  for(j=0;j<Ndist;j++) {
	    for(band=0;band<Nbands;band++) {
	      if(soil_con.AreaFract[band]>0) {
		for(i=0;i<options.Nlayer+3;i++) {
		  debug.inflow[j][band][i]  = 0;
		  debug.outflow[j][band][i] = 0;
		}
	      }
	    }
	  }
	}
      }

      /**************************************************
        Initialize Energy Balance Terms
      **************************************************/

      if(options.FULL_ENERGY || options.SNOW_MODEL) {
	/** Clear Energy Balance Terms **/
	for(band=0;band<Nbands;band++) {
	  if(soil_con.AreaFract[band]>0) {
	    energy[iveg][band].shortwave       = atmos->shortwave;
	    energy[iveg][band].longwave        = 0.;
	    energy[iveg][band].latent          = 0.;
	    energy[iveg][band].sensible        = 0.;
	    energy[iveg][band].grnd_flux       = 0.;
	    energy[iveg][band].deltaH          = 0.;
	    energy[iveg][band].albedo          = 0.;
	    energy[iveg][band].error           = 0.;
	    energy[iveg][band].deltaCC         = 0.;
	    energy[iveg][band].snow_flux       = 0.;
	    energy[iveg][band].refreeze_energy = 0.;
	    energy[iveg][band].advection       = 0.;
	    energy[iveg][band].Trad[0]         = 0.;
	    energy[iveg][band].Trad[1]         = 0.;
	  }
	}
      }
      if(options.SNOW_MODEL) {
	/** Initialize Snow Vapor Flux Terms **/
	for(band=0;band<Nbands;band++) {
	  if(soil_con.AreaFract[band]>0) {
	    snow[iveg][band].vapor_flux        = 0.;
	    snow[iveg][band].canopy_vapor_flux = 0.;
	    snow_inflow[band]                  = 0.;
	    Melt[band*2]                       = 0.;
	  }
	}
      }

      if(options.FULL_ENERGY) {
        /** Set Damping Depth **/
        dp=soil_con.dp;

        /** Compute Total Moisture of Surface Layer **/
	prepare_full_energy(iveg, Nveg, gp.Nnodes, prcp, 
			    soil_con, &moist, &ice0);
      }

      /** Compute Bare Soil (free of snow) Albedo **/
      if(iveg!=Nveg) bare_albedo 
		       = veg_lib[veg_class].albedo[dmy[rec].month-1];
      else bare_albedo = atmos->albedo;
      
      /** Compute the aerodynamic resistance **/
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
                      tmp_wind,cell[WET][iveg][0].aero_resist);
  
      /** Check for Presence of Snow Pack, or New Snow **/
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
  
        /**************************************************
           Full Energy Balance Model
         **************************************************/
  
/* 	for(j=0;j<Ndist;j++) { */
/* 	  rainfall[j] = rainonly[j+iveg*2]; */
/* 	  snowfall[j] = prec[j+iveg*2] - rainonly[j+iveg*2]; */
/* 	} */

	/*******************************************************************
          Solve Snow Accumulation and Melt on the Ground and in the Canopy 
        *******************************************************************/

	for(band=0;band<Nbands;band++) {
	  if(soil_con.AreaFract[band]>0) {
	    band_energy[0] = energy[iveg][band];
	    if(iveg==Nveg) sveg=0;
	    else sveg=iveg;
	    Melt[band*2] 
	      = solve_snow(&snow[iveg][band],cell[WET][iveg][band].layer,
			   cell[DRY][iveg][band].layer,
			   &veg_var[WET][sveg][band],
			   &veg_var[DRY][sveg][band],dmy,band_energy,soil_con,
			   overstory,&SNOW,&SOLVE_SURF_ENERGY,gp.dt,rec,
			   veg_class,iveg,Nveg,band,dmy[rec].hour,
			   gp.Nnodes,
			   atmos->shortwave,atmos->longwave,atmos->air_temp,
			   prec[iveg*2],atmos->density,atmos->vp,atmos->vpd,
			   atmos->pressure,prcp[0].mu[iveg],roughness,
			   displacement,ref_height,surf_atten,gp.MAX_SNOW_TEMP,
			   gp.MIN_RAIN_TEMP,gp.wind_h,moist,ice0,dp,
			   bare_albedo,&Le,&Ls,cell[WET][iveg][0].aero_resist,
			   tmp_wind,&step_net_short,&out_short,&step_rad,
			   &Evap,&last_T1,&Tsurf,&Tgrnd,&Tend_surf,&Tend_grnd,
			   &tmp_snow_energy,&snow_inflow[band],&ppt[band*2]);
	    
	    tmp_surf_temp = 0.;
 
	    if(!options.SNOW_MODEL || 
	       (options.SNOW_MODEL && SOLVE_SURF_ENERGY)) {
    
	      /**********************************************************
	        Solve Energy Balance Components for Ground Free of Snow
	      **********************************************************/
	      
	      /** Make Copies of Vegetation, Soil Layer, and Energy 
		  Balance Variables before Calculations **/
	      for(j=0;j<Ndist;j++) {
		if(iveg<Nveg) {
		  tmp_veg_var[j] = veg_var[j][iveg][band];
		}
		else {
		  tmp_veg_var[j] = veg_var[j][0][band];
		}
		for(i=0;i<options.Nlayer;i++) 
		  tmp_layer[j][i] = cell[j][iveg][band].layer[i];
	      }
	      tmp_energy[0]  = band_energy[0];
	      roughness      = soil_con.rough;
	      rainfall[WET]  = prec[iveg*2] * soil_con.Pfactor[band];
	      rainfall[DRY]  = prec[iveg*2+1] * soil_con.Pfactor[band];
	      Tend_grnd 
		= calc_surf_energy_bal(!SNOW,rec,iveg,options.Nlayer,
				       Nveg,
				       gp.dt,gp.Nnodes,veg_class,
				       band,ice0,moist,dp,surf_atten,
				       energy[iveg][band].T[0],
				       atmos->shortwave,tmp_energy->longwave,
				       Le,Ls,prcp[0].mu[iveg],
				       0.,roughness,ref_height+roughness,
				       tmp_snow_energy,
				       snow[iveg][band].vapor_flux,
				       &last_T1,cell[WET][iveg][0].aero_resist,
				       tmp_wind,rainfall,atmos,
				       &tmp_veg_var[WET],&tmp_veg_var[DRY],
				       tmp_energy,tmp_layer[WET],
				       tmp_layer[DRY],soil_con,dmy[rec]);  
	      
	      tmp_energy[0].longwave  
		-= STEFAN_B * pow(tmp_energy[0].Trad[0]+KELVIN,4.0);  
	      tmp_energy[0].Trad[0] = Tsurf;
	      tmp_energy[0].Trad[1] = Tgrnd;
	      
	      for(j=0;j<Ndist;j++) {
		for(i=0;i<options.Nlayer;i++) 
		  cell[j][iveg][band].layer[i]=tmp_layer[j][i];
		if(!SNOW) {
		  if(iveg<Nveg) {
		    veg_var[j][iveg][band] = tmp_veg_var[j];
		    ppt[band*2+j]      = veg_var[j][iveg][band].throughfall; 
		  }
		  else ppt[band*2+j] = prec[iveg*2+j] 
			 * soil_con.Pfactor[band]; 
		}
	      }
	    
	      band_energy[0] = tmp_energy[0];
	    } /* finished computing energy balance for no snow cover */	
	    else { 
	      band_energy[0].longwave  
		-= STEFAN_B * pow(Tsurf+KELVIN,4.0); 
	      band_energy[0].albedo = snow[iveg][band].albedo;
	    } 

	    energy[iveg][band] = band_energy[0];
	    if(overstory && snow[iveg][band].snow)
	      energy[iveg][band].shortwave     *= surf_atten; 
	    
	    /*****************************************
	      Compute Ground Thermal Fluxes Based on 
	      Average Soil Surface Temperature
	    *****************************************/

	    if(options.FULL_ENERGY && !options.FROZEN_SOIL) { 
	      /** Use Xu Liang's Equations **/
	      energy[iveg][band].T[1] = last_T1; 
	    } 
	    energy[iveg][band].T[0] = Tend_grnd;	  
	  
	    /********************************************************
	      Compute Runoff, Baseflow, and Soil Moisture Transport
	    ********************************************************/

	    for(j=0;j<Ndist;j++) 
	      cell[j][iveg][band].inflow = ppt[band*2+j];
	    runoff(cell[WET][iveg][band].layer,cell[DRY][iveg][band].layer,
		   &energy[iveg][band],soil_con,&cell[WET][iveg][band].runoff,
		   &cell[DRY][iveg][band].runoff,
		   &cell[WET][iveg][band].baseflow,
		   &cell[DRY][iveg][band].baseflow,&ppt[band*2],
		   prcp[0].mu[iveg],
		   gp.dt,gp.Nnodes,band,rec,iveg);
	    
	    out_short = atmos->albedo * atmos->shortwave;
	    step_net_short = (1.0 - atmos->albedo) * atmos->shortwave;
	  }
	} /** End Loop Through Elevation Bands **/
      } /** End Full Energy Balance Model **/

      else {

        /******************************************************
	  Water Balance Model
         ******************************************************/

	ANY_SNOW = FALSE;
	for(band=0;band<options.SNOW_BAND;band++)
	  if(snow[iveg][band].snow 
	     || snow[iveg][band].snow_canopy>0.) ANY_SNOW = TRUE;
        if(options.SNOW_MODEL && (ANY_SNOW || atmos->prec-atmos->rainonly>0)) {

          for(j=0;j<Ndist;j++) {
	    tmp_layerevap[j]           = (double*)calloc(options.Nlayer,
							 sizeof(double));
	    tmp_throughfall[j]         = (double*)calloc(Nbands,
							 sizeof(double));
	    tmp_canopyevap[j]          = (double*)calloc(Nbands,
							 sizeof(double));
	  }
	  tmp_canopy_vapor_flux        = (double*)calloc(Nbands,
							 sizeof(double));
	  tmp_vapor_flux               = (double*)calloc(Nbands,
							 sizeof(double));
          tmp_melt                     = (double*)calloc(Nbands*2,
							 sizeof(double));
          tmp_ppt                      = (double*)calloc(Nbands*2,
							 sizeof(double));

	  if (tmp_wind[0] > 0.0)
	    tmp_aero_resist 
	      = cell[WET][iveg][0].aero_resist[0] 
	      / StabilityCorrection(ref_height, displacement, 
				    atmos->air_temp, atmos->air_temp, 
				    tmp_wind[0], roughness);
	  else
	    tmp_aero_resist = HUGE_RESIST;
	  
	  if(debug.PRT_FLUX && rec==0) {
	    fprintf(debug.fg_energy,"DATE\tVEG TYPE\tPRCPDIST\tNET SHT\t");
	    fprintf(debug.fg_energy,"NET LNG\tLATENT\tSENSBL\tADVEC\t");
	    fprintf(debug.fg_energy,"del CC\tMELT\t");
	    fprintf(debug.fg_energy,"T_AIR\tT_SNOW\tWIND\n");
	  }
	  
	  /***********************************************
	    New Snow, or Old Snow on Ground or in Canopy
	  ***********************************************/

	  HourlyT(1,tmax_hour,tmax,tmin_hour,tmin,air_temp);
	  
	  for(hour=0;hour<gp.dt/SNOW_STEP;hour++) {
	    
	    /** Snow melt model is solved using a time step of SNOW_STEP **/
	    
	    /** Initialize parameters for the shorter time step **/
	    step_air_temp = air_temp[hour];
/* 	    if(step_air_temp<tmin[1] && hour>tmax_hour[1]) { */
/* 	      tmp_vp     = atmos->vp; */
/* 	      tmp_vpd    = atmos->vpd; */
/* 	      step_vp  = atmos[1].vp; */
/* 	      step_vpd = atmos[1].vpd; */
/* 	    } */
/* 	    else if(step_air_temp<tmin[1] && hour<tmin_hour[1]) { */
/* 	      tmp_vp     = atmos->vp; */
/* 	      tmp_vpd    = atmos->vpd; */
/* 	      step_vp  = atmos[-1].vp; */
/* 	      step_vpd = atmos[-1].vpd; */
/* 	    } */
/* 	    else if(step_air_temp<tmin[1]) { */
/* 	      vicerror("Air temperature falls below daily minimum temperature.  Please report error."); */
/* 	    } */
/* 	    else { */
/* 	      step_vp  = atmos->vp; */
/* 	      step_vpd = atmos->vpd; */
/* 	    } */
	    step_vp  = svp(step_air_temp) - atmos->vpd; 
	    step_vpd = atmos->vpd; 
	    
	    step_density = 3.486*atmos->pressure/(275.0
						    + step_air_temp);
	    for(band=0;band<Nbands;band++) {
	      if(soil_con.AreaFract[band]>0) {
		snow[iveg][band].vapor_flux         = 0.;
		snow[iveg][band].canopy_vapor_flux  = 0.;
		if(iveg<Nveg) {
		  veg_var[WET][iveg][band].canopyevap = 0.;
		  veg_var[DRY][iveg][band].canopyevap = 0.;
		}
		for(i=0;i<Ndist;i++) 
		  for(j=0;j<options.Nlayer;j++) 
		    cell[i][iveg][band].layer[j].evap = 0.;
	      }
	    }
	    calc_long_shortwave(&inshort,&inlong,&atmos->tskc,
				step_air_temp,step_vp,
				(double)soil_con.time_zone_lng,
				(double)soil_con.lng,(double)soil_con.lat,
				(double)dmy[rec].day_in_year,
				(double)hour*(double)SNOW_STEP,
				FALSE,FALSE,TRUE);

	    for(band=0;band<Nbands;band++) {	      
	      if(soil_con.AreaFract[band]>0) {
		energy[iveg][band].refreeze_energy = 0.;
		if(iveg==Nveg) sveg=0;
		else sveg=iveg;
		tmp_melt[band*2] 
		  = solve_snow(&snow[iveg][band],cell[WET][iveg][band].layer,
			       cell[DRY][iveg][band].layer,
			       &veg_var[WET][sveg][band],
			       &veg_var[DRY][sveg][band],dmy,
			       &energy[iveg][band],
			       soil_con,overstory,&SNOW,&SOLVE_SURF_ENERGY,
			       SNOW_STEP,rec,veg_class,iveg,Nveg,band,hour,
			       gp.Nnodes,inshort,inlong,
			       step_air_temp,prec[iveg*2]*(double)SNOW_STEP
			       /(double)gp.dt,step_density,step_vp,
			       step_vpd,atmos->pressure,prcp[0].mu[iveg],
			       roughness,displacement,ref_height,surf_atten,
			       gp.MAX_SNOW_TEMP,gp.MIN_RAIN_TEMP,gp.wind_h,
			       moist,ice0,dp,bare_albedo,&Le,&Ls,
			       cell[WET][iveg][0].aero_resist,tmp_wind,
			       &step_net_short,&out_short,&step_rad,
			       &Evap,&last_T1,&Tsurf,&Tgrnd,&Tend_surf,
			       &Tend_grnd,&tmp_snow_energy,
			       &snow_inflow[band],&ppt[band*2]);
		
		if(!SNOW) {
		  
		  /**********************************************
	            Compute Evaporation for Ground Free of Snow
		  **********************************************/

		  if(iveg<Nveg) {
		    if(veg_lib[veg_class].LAI[dmy[rec].month-1] > 0.0) {
		      for(j=0;j<Ndist;j++) {
			tmp_Wdew[j] = veg_var[j][iveg][band].Wdew;
			rainfall[j] 
			  = prec[j+iveg*2] * soil_con.Pfactor[band] 
			  * (double)SNOW_STEP / (double)gp.dt; 
		      }
		      Evap = canopy_evap(cell[WET][iveg][band].layer,
					 cell[DRY][iveg][band].layer,
					 &veg_var[WET][iveg][band],
					 &veg_var[DRY][iveg][band],
					 TRUE,veg_class,dmy[rec].month,
					 prcp[0].mu[iveg],tmp_Wdew,
					 step_air_temp
					 + soil_con.Tfactor[band],
					 (double)SNOW_STEP,step_rad,
					 step_vpd,step_net_short,step_air_temp
					 + soil_con.Tfactor[band],
					 tmp_aero_resist,
					 displacement,roughness,ref_height,
					 soil_con.elevation,rainfall,
					 soil_con.depth,
					 soil_con.Wcr,soil_con.Wpwp);
		      for(j=0;j<Ndist;j++) 
			ppt[band*2+j] += veg_var[j][iveg][band].throughfall; 
		    }
		    else {
		      Evap = arno_evap(cell[WET][iveg][band].layer, 
				       cell[DRY][iveg][band].layer, 
				       step_rad, step_air_temp, 
				       step_vpd, 
				       step_net_short, soil_con.depth[0], 
				       soil_con.max_moist[0], 
				       soil_con.elevation, 
				       soil_con.b_infilt, step_air_temp 
				       + soil_con.Tfactor[band],
				       displacement, roughness, ref_height, 
				       tmp_aero_resist, (double)SNOW_STEP,
				       prcp[0].mu[iveg]);
		      for(j=0;j<Ndist;j++) {
			ppt[band*2+j] = prec[iveg*2+j] 
			  * soil_con.Pfactor[band] 
			  * (double)SNOW_STEP / (double)gp.dt;
			veg_var[j][iveg][band].throughfall = ppt[band*2+j]; 
		      }
		    }  /** end if vegetation not growing **/
		  }  /** end if vegetation **/
	    
		  /** Compute Evaporation from Bare Soil **/
		  else {
		    Evap = arno_evap(cell[WET][iveg][band].layer, 
				     cell[DRY][iveg][band].layer, 
				     step_rad, step_air_temp, step_vpd, 
				     step_net_short, soil_con.depth[0],
				     soil_con.max_moist[0], soil_con.elevation,
				     soil_con.b_infilt, 
				     step_air_temp, displacement, roughness, 
				     ref_height, tmp_aero_resist, 
				     (double)SNOW_STEP,prcp[0].mu[iveg]);
		    for(j=0;j<Ndist;j++) {
		      ppt[band*2+j] = prec[iveg*2+j] * soil_con.Pfactor[band]
			* (double)SNOW_STEP / (double)gp.dt;
		    }
		  }
		} /* End Computations for Bare Ground */
		
		/** Store Huorly Moisture Flux Terms **/
		if(iveg<Nveg) {
		  tmp_throughfall[WET][band]  
		    += veg_var[WET][iveg][band].throughfall;
		  tmp_canopyevap[WET][band]   
		    += veg_var[WET][iveg][band].canopyevap;
		  tmp_canopy_vapor_flux[band] 
		    += snow[iveg][band].canopy_vapor_flux;
		}
		tmp_vapor_flux[band]        
		  += snow[iveg][band].vapor_flux;
		if(options.DIST_PRCP && hour==0) {
		  tmp_throughfall[DRY][band] 
		    += veg_var[DRY][iveg][band].throughfall;
		  tmp_canopyevap[DRY][band] 
		    += veg_var[DRY][iveg][band].canopyevap;
		}
		
		for(j=0;j<Ndist;j++) {
		  Melt[band*2+j] += tmp_melt[band*2+j];
		  tmp_ppt[band*2+j]  += ppt[band*2+j];
		}
		
	      }
	    } /* End Loop Through Elevation Bands */

/* 	    if(atmos->air_temp<tmin[1] && hour>tmax_hour[1]) { */
/* 	      atmos->vp  = tmp_vp; */
/* 	      atmos->vpd = tmp_vpd; */
/* 	    } */
	    
	  } /* End Hourly Snow Solution */
	  cell[WET][iveg][0].aero_resist[0] = tmp_aero_resist;
	  for(band=0;band<Nbands;band++) {
	    if(soil_con.AreaFract[band]>0) {
	      snow[iveg][band].canopy_vapor_flux = tmp_canopy_vapor_flux[band];
	      snow[iveg][band].vapor_flux = tmp_vapor_flux[band];
	    }
	  }
	  for(j=0;j<Ndist;j++) {
	    for(band=0;band<Nbands;band++) {
	      if(soil_con.AreaFract[band]>0) {
		if(iveg < Nveg) {
		  veg_var[j][iveg][band].canopyevap  
		    = tmp_canopyevap[j][band];
		  veg_var[j][iveg][band].throughfall 
		    = tmp_throughfall[j][band];
		}
		ppt[band*2+j] = tmp_ppt[band*2+j];
		for(i=0;i<options.Nlayer;i++)
		  cell[j][iveg][band].layer[i].evap += tmp_layerevap[j][i];
	      }
	    }  /** loop through elevation bands **/
	  }  /** Loop through wet and dry fractions **/

	  for(j=0;j<Ndist;j++) {
	    free((char*)tmp_layerevap[j]);
	    free((char*)tmp_throughfall[j]);
	    free((char*)tmp_canopyevap[j]);
	  }
	  free((char*)tmp_canopy_vapor_flux);
	  free((char*)tmp_vapor_flux);
	  free((char*)tmp_melt);
	  free((char*)tmp_ppt);
	  
	  ANY_SNOW = TRUE;
	  
        }
        else {
          SNOW = FALSE;
        }
	
	if(!ANY_SNOW) {
	  
	  /******************************************
	    No Snow on Ground Or Snow Model not Run
	  ******************************************/
	  
	  for(j=0;j<Ndist;j++) 
	    for(band=0;band<Nbands;band++) ppt[band*2+j] = 0.;
	  if (tmp_wind[0] > 0.0)
	    cell[j][iveg][0].aero_resist[0] 
	      /= StabilityCorrection(ref_height, displacement, 
				     atmos->air_temp, atmos->air_temp, 
				     tmp_wind[0], roughness);
	  else
	    cell[j][iveg][0].aero_resist[0] = HUGE_RESIST;

	  for(band=0;band<Nbands;band++) {
	    if(soil_con.AreaFract[band]>0) {

	      /** Compute Net Radiation at Surface **/
	      step_net_short = (1. - bare_albedo) * atmos->shortwave; 
	      step_rad = step_net_short + atmos->longwave;

	      if(iveg!=Nveg) {
		/** Compute Evaporation from Vegetation **/
		if(veg_lib[veg_class].LAI[dmy[rec].month-1] > 0.0) {
		  for(j=0;j<Ndist;j++) {
		    tmp_Wdew[j] = veg_var[j][iveg][band].Wdew;
		    rainfall[j] = prec[iveg*2+j] * soil_con.Pfactor[band];
		  }
		  Evap = canopy_evap(cell[WET][iveg][band].layer,
				     cell[DRY][iveg][band].layer,
				     &veg_var[WET][iveg][band],
				     &veg_var[DRY][iveg][band],
				     TRUE,veg_class,dmy[rec].month,
				     prcp[0].mu[iveg],tmp_Wdew,atmos->air_temp
				     + soil_con.Tfactor[band],
				     (double)gp.dt,step_rad,atmos->vpd,
				     step_net_short,atmos->air_temp
				     + soil_con.Tfactor[band],
				     cell[WET][iveg][0].aero_resist[0],
				     displacement,roughness,ref_height,
				     soil_con.elevation,rainfall,
				     soil_con.depth,
				     soil_con.Wcr,soil_con.Wpwp);
		  for(j=0;j<Ndist;j++) 
		    ppt[band*2+j] += veg_var[j][iveg][band].throughfall; 
		}
		else {
		  Evap = arno_evap(cell[WET][iveg][band].layer, 
				   cell[DRY][iveg][band].layer, 
				   step_rad, atmos->air_temp, atmos->vpd, 
				   step_net_short, soil_con.depth[0], 
				   soil_con.max_moist[0], soil_con.elevation, 
				   soil_con.b_infilt, atmos->air_temp 
				   + soil_con.Tfactor[band],
				   displacement, roughness, ref_height, 
				   cell[WET][iveg][0].aero_resist[0], 
				   (double)gp.dt,
				   prcp[0].mu[iveg]);
		  for(j=0;j<Ndist;j++) {
		    ppt[band*2+j] = prec[iveg*2+j] * soil_con.Pfactor[band];
		    veg_var[j][iveg][band].throughfall = ppt[band*2+j]; 
		  }
		}
	      }
	    
	      /** Compute Evaporation from Bare Soil **/
	      else {
		Evap = arno_evap(cell[WET][iveg][band].layer, 
				 cell[DRY][iveg][band].layer, 
				 step_rad, atmos->air_temp, atmos->vpd, 
				 step_net_short, soil_con.depth[0],
				 soil_con.max_moist[0], soil_con.elevation,
				 soil_con.b_infilt, 
				 atmos->air_temp, displacement, roughness, 
				 ref_height, 
				 cell[WET][iveg][0].aero_resist[0], 
				 (double)gp.dt,prcp[0].mu[iveg]);
		for(j=0;j<Ndist;j++) {
		  ppt[band*2+j] = prec[iveg*2+j] * soil_con.Pfactor[band];
		}
	      }
	    }
	  } /** End loop through elevation bands **/
        }
	
	for(band=0;band<Nbands;band++) {
	  if(soil_con.AreaFract[band]>0) {
	    if(!options.SNOW_MODEL) 
	      for(j=0;j<Ndist;j++) ppt[band*2+j] += melt[j];
	    
	    /********************************************************
	      Compute Runoff, Baseflow, and Soil Moisture Transport
	    ********************************************************/
	    
	    for(j=0;j<Ndist;j++) cell[j][iveg][band].inflow = ppt[band*2+j];
	    runoff(cell[WET][iveg][band].layer, cell[DRY][iveg][band].layer, 
		   &energy[iveg][band],soil_con, 
		   &cell[WET][iveg][band].runoff, 
		   &cell[DRY][iveg][band].runoff, 
		   &cell[WET][iveg][band].baseflow,
		   &cell[DRY][iveg][band].baseflow,&ppt[band*2],
		   prcp[0].mu[iveg],gp.dt,
		   gp.Nnodes,band,rec,iveg);
	  }
	}  /** End loop through elevation bands **/
      } /** End Water Balance Model **/
  
      /****************************
	Controls Debugging Output
      ****************************/

      for(j=0;j<Ndist;j++) {
	if(iveg<Nveg) tmp_veg[j] = veg_var[j][iveg];
	else tmp_veg[j] = NULL;
      }
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
        ptr_energy = energy[iveg];
        tmp_snow = snow[iveg];
      }
      for(j=0;j<Ndist;j++) {
	if(j==0) tmp_mu = prcp[0].mu[iveg]; 
	else tmp_mu = 1. - prcp[0].mu[iveg]; 
	/** for debugging water balance: [0] = vegetation, 
	    [1] = ground snow, [2..Nlayer+1] = soil layers **/
	if(debug.PRT_BALANCE) {
	  for(band=0;band<Nbands;band++) {
	    if(soil_con.AreaFract[band]>0) {
	      debug.inflow[j][band][options.Nlayer+2] 
		+= prec[j+iveg*2] * soil_con.Pfactor[band];
	      debug.inflow[j][band][0]  = 0.;
	      debug.inflow[j][band][1]  = 0.;
	      debug.outflow[j][band][0] = 0.;
	      debug.outflow[j][band][1] = 0.;
	      if(iveg<Nveg) {
		/** Vegetation Present **/
		debug.inflow[j][band][0] += prec[j+iveg*2] 
		  * soil_con.Pfactor[band]; 
		debug.outflow[j][band][0] 
		  += veg_var[j][iveg][band].throughfall;
	      }
	      if(j==0)
		debug.inflow[j][band][1] += snow_inflow[band];
	      if(options.FULL_ENERGY || options.SNOW_MODEL) {
		debug.outflow[j][band][1] += Melt[band*2+j];
	      }
	    }
	  }  /** End loop through elevation bands **/
	}
         
	if(iveg!=Nveg) {
	  write_debug(atmos[0],soil_con,cell[j][iveg],
		      ptr_energy,tmp_snow,tmp_veg[j],dmy[rec],gp,
		      out_short,tmp_mu,Nveg,
		      iveg,rec,gridcell,j,NEWCELL); 
	}
	else {
	  write_debug(atmos[0],soil_con,cell[j][iveg],
		      ptr_energy,tmp_snow,tmp_veg[j],dmy[rec],gp,
		      out_short,tmp_mu,Nveg,iveg,rec,gridcell,
		      j,NEWCELL); 
	}
      }
    } /** end if RUN_TYPE **/

    atmos[0] = tmp_atmos;

  } /** end of vegetation loop **/

  /** END Temperature and Moisture Profile Debugging Output **/
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    for(j=0;j<2;j++) {
      free ((char *)tmp_layer[j]);
    }
    free ((char *)tmp_energy);
  }
  if(options.SNOW_MODEL || options.FULL_ENERGY) {
    free((char*)snow_inflow);
  }
  free((char*)band_energy);
  free((char*)ppt);
  free((char*)Melt);

}

void store_max_min_temp(atmos_data_struct *atmos,
                        double *tmax,
			int    *tmax_hour,
                        double *tmin,
			int    *tmin_hour,
                        int rec,
                        int Nrecs) {
/**********************************************************************
  This subroutine sets prepares arrays of maximum and minimum 
  daily air temperature.
**********************************************************************/

  static int last_rec;

  if(rec==0) {
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
  tmax_hour[0] = tmax_hour[1] = tmax_hour[2] = 15;
  tmin_hour[0] = tmin_hour[1] = tmin_hour[2] = 5;
}

void store_moisture_for_debug(int                 iveg,
		              int                 Nveg,
			      double             *mu,
			      cell_data_struct ***cell,
			      veg_var_struct   ***veg_var,
			      snow_data_struct  **snow,
			      soil_con_struct     soil_con) {
/****************************************************************
  This subroutine was written to save the current water storage
  terms for use in calculating the model water balance error
****************************************************************/

  extern option_struct options;
  extern debug_struct  debug;

  int               Ndist;
  int               i;
  int               band;
  int               dist;
  int               Nbands;
  layer_data_struct tmp_layer;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = options.SNOW_BAND;

  for(band=0;band<Nbands;band++) {
    if(soil_con.AreaFract[band]>0) {
      for(dist=0;dist<Ndist;dist++) 
	for(i=0;i<options.Nlayer+3;i++)
	  debug.store_moist[dist][band][i] = 0.;
      if(iveg<Nveg) {
	for(dist=0;dist<Ndist;dist++)
	  debug.store_moist[dist][band][0] 
	    += veg_var[dist][iveg][band].Wdew;
	if(options.FULL_ENERGY || options.SNOW_MODEL)
	  debug.store_moist[WET][band][0] += (snow[iveg][band].snow_canopy) 
	    * 1000.;
      }
      for(dist=0;dist<Ndist;dist++)
	debug.store_moist[dist][band][options.Nlayer+2] 
	  += debug.store_moist[dist][band][0];
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
	debug.store_moist[WET][band][1] += (snow[iveg][band].swq*1000.);
	for(dist=0;dist<Ndist;dist++)
	  debug.store_moist[dist][band][options.Nlayer+2] 
	    += debug.store_moist[dist][band][1];
      }
      for(i=0;i<options.Nlayer;i++) {
	for(dist=0;dist<Ndist;dist++) {
	  debug.store_moist[dist][band][i+2] 
	    = find_total_layer_moisture(cell[dist][iveg][band].layer[i], 
					soil_con.depth[i]);
	  debug.store_moist[dist][band][options.Nlayer+2] 
	    += debug.store_moist[dist][band][i+2];
	}
      }
    }
  }
}

void prepare_full_energy(int               iveg,
			 int               Nveg,
			 int               Nnodes,
			 dist_prcp_struct *prcp,
			 soil_con_struct   soil_con,
			 double           *moist,
			 double           *ice0) {
/*******************************************************************
  This subroutine prepares variables for use with the full energy
  balance model.
*******************************************************************/

  extern option_struct options;

  int                i, band;
  double            *null_ptr;
  layer_data_struct  layer[options.Nlayer];

  for(band=0;band<options.SNOW_BAND;band++) {
    for(i=0;i<options.Nlayer;i++) 
      layer[i] = find_average_layer(prcp[0].cell[WET][iveg][band].layer[i],
				    prcp[0].cell[DRY][iveg][band].layer[i],
				    soil_con.depth[i], prcp[0].mu[iveg]);
    
    moist[0] = find_total_layer_moisture(layer[0],soil_con.depth[0]);
    moist[0] /= soil_con.depth[0]*1000.;
    if(options.FROZEN_SOIL){
      if((prcp[0].energy[iveg][band].T[0]+prcp[0].energy[iveg][band].T[1])/2.<0.) {
	ice0[0] = moist[0] 
	  - maximum_unfrozen_water((prcp[0].energy[iveg][band].T[0]
				    + prcp[0].energy[iveg][band].T[1])/2.,
				   soil_con.max_moist[0]/(soil_con.depth[0]*1000.),
				   soil_con.bubble,soil_con.expt[0]);
	if(ice0[0]<0.) ice0[0]=0.;
	if(ice0[0]>soil_con.max_moist[0]/(soil_con.depth[0]*1000.)) 
	  ice0[0]=soil_con.max_moist[0]/(soil_con.depth[0]*1000.);
      }
      else ice0[0]=0.;
    }
    else ice0[0]=0.;

    if(options.FULL_ENERGY && !options.FROZEN_SOIL) {
      /** Compute Soil Thermal Properties if not using Frozen Soils **/
      null_ptr = NULL;
      soil_thermal_calc(soil_con, layer, prcp[0].energy[iveg][band], null_ptr,
			null_ptr, null_ptr, options.Nlayer,
			Nnodes);
      
      for(i=0;i<options.Nlayer;i++) {
	prcp[0].cell[WET][0][band].layer[i].kappa = layer[i].kappa;
	prcp[0].cell[DRY][0][band].layer[i].kappa = layer[i].kappa;
	prcp[0].cell[WET][0][band].layer[i].Cs = layer[i].Cs;
	prcp[0].cell[DRY][0][band].layer[i].Cs = layer[i].Cs;
      }
      
      /** Save Thermal Conductivities for Energy Balance **/
      prcp[0].energy[iveg][band].kappa[0] = layer[0].kappa;
      prcp[0].energy[iveg][band].Cs[0] = layer[0].Cs;
      prcp[0].energy[iveg][band].kappa[1] = layer[1].kappa;
      prcp[0].energy[iveg][band].Cs[1] = layer[1].Cs;
    }
  }
}
