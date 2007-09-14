#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  full_energy(char                 NEWCELL,
		 int                  gridcell,
                 int                  rec,
                 atmos_data_struct   *atmos,
                 dist_prcp_struct    *prcp,
                 dmy_struct          *dmy,
                 global_param_struct *gp,
		 lake_con_struct     *lake_con,
                 soil_con_struct     *soil_con,
                 veg_con_struct      *veg_con)
/**********************************************************************
	full_energy	Keith Cherkauer		January 8, 1997

  This subroutine controls the model core, it solves both the energy
  and water balance models, as well as frozen soils.  

  modifications:
  07-98 restructured to fix problems with distributed precipitation, 
        and to add the ability to solve the snow model at different 
	elevation bands within a single grid cell.                 KAC
  01-19-00 modified to work with the new atmosphere data structure 
           implemented when the radiation forcing routines were 
	   updated.  Also modified to use the new simplified
	   soil moisture storage for the frozen soil algorithm.    KAC
  12-01-00 modified to include the lakes and wetlands algorithm.   KAC
  11-18-02 Modified to handle blowing snow.  Also added debugging
           output for lake model.                                  LCB
  05-27-03 Updated passing of veg_con parameters for blowing snow
           to surface_fluxes.  Original did not account for the fact
           that veg_con is not allocated for veg = Nveg (bare soil)
           case.  This eliminates a memory error.                  KAC
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    used in flux calculations.				TJB
  2006-Sep-23 Implemented flexible output configuration; now computation
	      of soil wetness and root zone soil moisture happens here. TJB
  2006-Nov-07 Removed LAKE_MODEL option. TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
           main subroutine, rather than ending the simulation.   GCT/KAC
  2007-May-01 Added case of SPATIAL_FROST = TRUE in modifications
              from 2006-Sep-23. GCT
  2007-Aug-10 Added features for EXCESS_ICE option.  JCA
               Including calculating subsidence for each layer and
               updating soil depth, effective porosity,
               bulk density, and soil moisture and fluxes by calling
               runoff function if subsidence occurs.
  2007-Sep-7 No longer resets ice content to previous time-step ice content if
               subsidence has occurred.  JCA
**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  char                   overstory;
  int                    i, j;
  int                    lidx;
  int                    Ndist;
  int                    dist;
  int                    iveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands;
  int                    ErrorFlag;
#if SPATIAL_FROST
  int                    frost_area;
#endif // SPATIAL_FROST
  double                 out_prec[2*MAX_BANDS];
  double                 out_rain[2*MAX_BANDS];
  double                 out_snow[2*MAX_BANDS];
  double                 out_short=0;
  double                 dp;
  double                 ice0;
  double                 moist;
  double                 surf_atten;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 wind_h;
  double                 height;
  double                 displacement[3];
  double                 roughness[3];
  double                 ref_height[3];
  double                 Cv;
  double                 Le;
  double                 Melt[2*MAX_BANDS];
  double                 bare_albedo;
  double                 snow_inflow[MAX_BANDS];
  double                 lake_prec;
  double                 rainonly;
  double                 sum_runoff;
  double                 sum_baseflow;
  double                 tmp_wind[3];
  double                 tmp_mu;
  double                 tmp_total_moist;
  double                 gauge_correction[2];
  float 	         lag_one;
  float 	         sigma_slope;
  float  	         fetch;
  lake_var_struct       *lake_var;
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;
  energy_bal_struct    **energy;
  energy_bal_struct     *ptr_energy;
  snow_data_struct     **snow;
  snow_data_struct      *tmp_snow;
  veg_var_struct        *tmp_veg[2];
  veg_var_struct        *wet_veg_var;
  veg_var_struct        *dry_veg_var;
  veg_var_struct         empty_veg_var;
#if EXCESS_ICE
  int                    SubsidenceUpdate = 0;
  int                    index;
  int                    MaxVeg;
  char                   ErrStr[MAXSTRING];
  double                 max_ice_layer; //mm/mm
  double                 ave_ice_fract; //mm/mm
  double                 ave_ice, tmp_ice; //mm
  double                 ice_layer; //mm
  double                 subsidence[MAX_LAYERS]; //mm
  double                 total_subsidence; //m
  double                 total_meltwater; //mm
  double                 tmp_depth, tmp_depth_prior; //m
  double                 ppt[2]; 
  double                 moist_prior[2][MAX_VEG][MAX_BANDS][MAX_LAYERS]; //mm
  double                 evap_prior[2][MAX_VEG][MAX_BANDS][MAX_LAYERS]; //mm
#endif //EXCESS_ICE

  /* set local pointers */
  cell    = prcp->cell;
  energy  = prcp->energy;
  lake_var = &prcp->lake_var;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;

  /* set variables for distributed precipitation */
  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  /* Set number of vegetation types */
  Nveg      = veg_con[0].vegetat_type_num;

  /** Set Damping Depth **/
  dp        = soil_con->dp;

  /* Compute gauge undercatch correction factors 
     - this assumes that the gauge is free of vegetation effects, so gauge
     correction is constant for the entire grid cell */
  if( options.CORRPREC && atmos->prec[NR] > 0 ) 
    correct_precip(gauge_correction, atmos->wind[NR], gp->wind_h, 
		   soil_con->rough, soil_con->snow_rough);
  else {
    gauge_correction[0] = 1;
    gauge_correction[1] = 1;
  }
  atmos->out_prec = 0;
  atmos->out_rain = 0;
  atmos->out_snow = 0;

  /**************************************************
    Solve Energy and/or Water Balance for Each
    Vegetation Type
  **************************************************/
  for(iveg = 0; iveg <= Nveg; iveg++){

    /** Solve Veg Type only if Coverage Greater than 0% **/
    if ((iveg <  Nveg && veg_con[iveg].Cv  > 0.) || 
	(iveg == Nveg && veg_con[0].Cv_sum < 1.)) {

      if ( iveg < Nveg ) Cv = veg_con[iveg].Cv;
      else Cv = 1. - veg_con[0].Cv_sum;

      /**************************************************
        Initialize Model Parameters
      **************************************************/

      for(band = 0; band < Nbands; band++) {
	if(soil_con->AreaFract[band] > 0) {

	  /* Initialize energy balance variables */
	  energy[iveg][band].shortwave = 0;
	  energy[iveg][band].longwave  = 0.;

	  /* Initialize snow variables */
	  snow[iveg][band].vapor_flux        = 0.;
	  snow[iveg][band].canopy_vapor_flux = 0.;
	  snow_inflow[band]                  = 0.;
	  Melt[band*2]                       = 0.;

	}
      }

      /* Initialize precipitation storage */
      for ( j = 0; j < 2*MAX_BANDS; j++ ) {
        out_prec[j] = 0;
        out_rain[j] = 0;
        out_snow[j] = 0;
      }
    
      /** Define vegetation class number **/
      if (iveg < Nveg) 
	veg_class = veg_con[iveg].veg_class;
      else 
	veg_class = 0;
      if ( iveg < Nveg ) wind_h = veg_lib[veg_class].wind_h;
      else wind_h = gp->wind_h;

      /** Compute Surface Attenuation due to Vegetation Coverage **/
      if ( iveg < Nveg )
	surf_atten = exp(-veg_lib[veg_class].rad_atten 
			 * veg_lib[veg_class].LAI[dmy[rec].month-1]);
      else 
	surf_atten = 1.;
        
      /* Initialize soil thermal properties for the top two layers */
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
	prepare_full_energy(iveg, Nveg, options.Nnode, prcp, 
			    soil_con, &moist, &ice0);
      }

      /** Compute Bare Soil (free of snow) Albedo **/
      if(iveg != Nveg) 
	bare_albedo = veg_lib[veg_class].albedo[dmy[rec].month-1];
      else 
	bare_albedo = BARE_SOIL_ALBEDO;
      
      /*************************************
	Compute the aerodynamic resistance 
      *************************************/

      /* Set surface descriptive variables */
      overstory = FALSE;
      if ( iveg < Nveg ) {
        displacement[0] = veg_lib[veg_class].displacement[dmy[rec].month-1];
        roughness[0]    = veg_lib[veg_class].roughness[dmy[rec].month-1];
        overstory       = veg_lib[veg_class].overstory;
	if ( roughness[0] == 0 ) roughness[0] = soil_con->rough;
      }
      if(iveg == Nveg || roughness == 0) {
        displacement[0] = 0.;
        roughness[0]    = soil_con->rough;
        overstory       = FALSE;
      }

      /* Initialize wind speeds */
      tmp_wind[0] = atmos->wind[NR];
      tmp_wind[1] = -999.;
      tmp_wind[2] = -999.;
 
      /* Estimate vegetation height */
      height = calc_veg_height(displacement[0]);

      /* Estimate reference height */
      if(displacement[0] < wind_h) 
	ref_height[0] = wind_h;
      else 
	ref_height[0] = displacement[0] + wind_h + roughness[0];

      /* Compute aerodynamic resistance over various surface types */
      ErrorFlag = CalcAerodynamic(overstory, height, veg_lib[veg_class].trunk_ratio, 
                      soil_con->snow_rough, soil_con->rough, 
		      veg_lib[veg_class].wind_atten,
		      gp->wind_h, cell[WET][iveg][0].aero_resist, tmp_wind, 
                      displacement, ref_height, roughness, 
                      Nveg, iveg);
      if ( ErrorFlag == ERROR ) return ( ERROR );  

      /**************************************************
        Store Water Balance Terms for Debugging
      **************************************************/

#if LINK_DEBUG
      if(debug.DEBUG || debug.PRT_MOIST || debug.PRT_BALANCE) {
        /** Compute current total moisture for water balance check **/
	store_moisture_for_debug(iveg, Nveg, prcp->mu, cell,
				 veg_var, snow, soil_con);
	if(debug.PRT_BALANCE) {
	  for(j=0; j<Ndist; j++) {
	    for(band=0; band<Nbands; band++) {
	      if(soil_con->AreaFract[band] > 0) {
		for(i=0; i<options.Nlayer+3; i++) {
		  debug.inflow[j][band][i]  = 0;
		  debug.outflow[j][band][i] = 0;
		}
	      }
	    }
	  }
	}
      }
#endif // LINK_DEBUG

      /******************************
        Solve ground surface fluxes 
      ******************************/
  
      for ( band = 0; band < Nbands; band++ ) {
	if( soil_con->AreaFract[band] > 0 ) {

	  if ( iveg < Nveg ) {
	    wet_veg_var = &(veg_var[WET][iveg][band]);
	    dry_veg_var = &(veg_var[DRY][iveg][band]);
	    lag_one     = veg_con[iveg].lag_one;
	    sigma_slope = veg_con[iveg].sigma_slope;
	    fetch       = veg_con[iveg].fetch;
	  }
	  else {
	    wet_veg_var = &(empty_veg_var);
	    dry_veg_var = &(empty_veg_var);
	    lag_one     = veg_con[0].lag_one;
	    sigma_slope = veg_con[0].sigma_slope;
	    fetch       = veg_con[0].fetch;
	  }

	  // set prior moist and ice for subsidence calculations
#if EXCESS_ICE
          for ( dist = 0; dist < Ndist; dist++ ) {
            for(lidx=0;lidx<options.Nlayer;lidx++) {
	      moist_prior[dist][iveg][band][lidx] = cell[dist][iveg][band].layer[lidx].moist;
	      evap_prior[dist][iveg][band][lidx] = 0; //initialize
	    }
	  }
#endif //EXCESS_ICE

	  ErrorFlag = surface_fluxes(overstory, bare_albedo, height, ice0, moist, 
#if EXCESS_ICE
				     SubsidenceUpdate, evap_prior[DRY][iveg][band], evap_prior[WET][iveg][band],
#endif
				     prcp->mu[iveg], surf_atten, &(Melt[band*2]), &Le, 
				     cell[WET][iveg][0].aero_resist, &(cell[WET][iveg][0].aero_resist_used),
				     &(cell[DRY][iveg][band].baseflow), 
				     &(cell[WET][iveg][band].baseflow), displacement, 
				     gauge_correction, &(cell[DRY][iveg][band].inflow), 
				     &(cell[WET][iveg][band].inflow), &out_prec[band*2], 
				     &out_rain[band*2], &out_snow[band*2],
				     ref_height, roughness, 
				     &(cell[DRY][iveg][band].runoff), 
				     &(cell[WET][iveg][band].runoff), &snow_inflow[band], 
				     tmp_wind, veg_con[iveg].root, Nbands, Ndist, 
				     options.Nlayer, Nveg, band, dp, iveg, rec, veg_class, 
				     atmos, dmy, &(energy[iveg][band]), gp, 
				     cell[DRY][iveg][band].layer, 
				     cell[WET][iveg][band].layer, &(snow[iveg][band]), 
				     soil_con, dry_veg_var, wet_veg_var, 
				     lag_one, sigma_slope, fetch);
	  
	  if ( ErrorFlag == ERROR ) return ( ERROR );
	  
	  atmos->out_prec += out_prec[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_rain += out_rain[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_snow += out_snow[band*2] * Cv * soil_con->AreaFract[band];

          /********************************************************
            Compute soil wetness and root zone soil moisture
          ********************************************************/
          // Loop through distributed precipitation fractions
          for ( dist = 0; dist < 2; dist++ ) {
            cell[dist][iveg][band].rootmoist = 0;
            cell[dist][iveg][band].wetness = 0;
            for(lidx=0;lidx<options.Nlayer;lidx++) {
#if SPATIAL_FROST
              tmp_total_moist = cell[dist][iveg][band].layer[lidx].moist;
              for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
                tmp_total_moist += cell[dist][iveg][band].layer[lidx].ice[frost_area];
	      }
#else
              tmp_total_moist = cell[dist][iveg][band].layer[lidx].moist + cell[dist][iveg][band].layer[lidx].ice;
#endif // SPATIAL_FROST 
              if (veg_con->root[lidx] > 0) {
                cell[dist][iveg][band].rootmoist += tmp_total_moist;
              }
#if EXCESS_ICE
	      cell[dist][iveg][band].wetness += (tmp_total_moist - soil_con->Wpwp[lidx])/(soil_con->effective_porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
#else
	      cell[dist][iveg][band].wetness += (tmp_total_moist - soil_con->Wpwp[lidx])/(soil_con->porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
#endif
            }
            cell[dist][iveg][band].wetness /= options.Nlayer;
          }

	} /** End Loop Through Elevation Bands **/
      } /** End Full Energy Balance Model **/


      /****************************
         Controls Debugging Output
      ****************************/
#if LINK_DEBUG
      
      for(j = 0; j < Ndist; j++) {
	if(iveg < Nveg) 
	  tmp_veg[j] = veg_var[j][iveg];
	else 
	  tmp_veg[j] = NULL;
      }
      ptr_energy = energy[iveg];
      tmp_snow   = snow[iveg];
      for(j = 0; j < Ndist; j++) {
	if(j == 0) 
	  tmp_mu = prcp->mu[iveg]; 
	else 
	  tmp_mu = 1. - prcp->mu[iveg]; 
	/** for debugging water balance: [0] = vegetation, 
	    [1] = ground snow, [2..Nlayer+1] = soil layers **/
	if(debug.PRT_BALANCE) {
	  for(band = 0; band < Nbands; band++) {
	    if(soil_con->AreaFract[band] > 0) {
	      debug.inflow[j][band][options.Nlayer+2] 
		+= out_prec[j+band*2] * soil_con->Pfactor[band];
	      debug.inflow[j][band][0]  = 0.;
	      debug.inflow[j][band][1]  = 0.;
	      debug.outflow[j][band][0] = 0.;
	      debug.outflow[j][band][1] = 0.;
	      if(iveg < Nveg) {
		/** Vegetation Present **/
		debug.inflow[j][band][0] += out_prec[j+band*2] * soil_con->Pfactor[band]; 
		debug.outflow[j][band][0] 
		  += veg_var[j][iveg][band].throughfall;
	      }
	      if(j == 0)
		debug.inflow[j][band][1] += snow_inflow[band];
	      debug.outflow[j][band][1] += Melt[band*2+j];
	    }
	  }  /** End loop through elevation bands **/
	}
         
	if(iveg != Nveg) {
	  write_debug(atmos, soil_con, cell[j][iveg], ptr_energy, 
		      tmp_snow, tmp_veg[j], &(dmy[rec]), gp, out_short, 
		      tmp_mu, Nveg, iveg, rec, gridcell, j, NEWCELL); 
	}
	else {
	  write_debug(atmos, soil_con, cell[j][iveg], ptr_energy, tmp_snow, 
		      tmp_veg[j], &(dmy[rec]), gp, out_short, tmp_mu, Nveg, 
		      iveg, rec, gridcell, j, NEWCELL); 
	}
      }
#endif // LINK_DEBUG

    } /** end current vegetation type **/
  } /** end of vegetation loop **/

  /****************************
     Calculate Subsidence
  ****************************/
#if EXCESS_ICE
  if(options.LAKES) 
    MaxVeg = Nveg+1;
  else
    MaxVeg = Nveg;
  
  total_subsidence = 0;
  total_meltwater = 0; //for lake model only
  for(lidx=0;lidx<options.Nlayer;lidx++) {//soil layer
    subsidence[lidx] = 0;
    if(soil_con->effective_porosity[lidx]>soil_con->porosity[lidx]){
      
      /* find average ice/porosity fraction and sub-grid with greatest ice/porosity fraction */
      ave_ice = 0;
      max_ice_layer = 0;
      for(band = 0; band < Nbands; band++) {//band
	if(soil_con->AreaFract[band] > 0) {
	  for(iveg = 0; iveg <= MaxVeg; iveg++){ //iveg  
	    if ((iveg <  Nveg && veg_con[iveg].Cv  > 0.) || 
		(iveg == Nveg && veg_con[0].Cv_sum < 1.) ||
		(iveg == MaxVeg && MaxVeg > Nveg && lake_con->Cl[0] > 0. )) {    
	      if ( iveg < Nveg ) Cv = veg_con[iveg].Cv;
	      else if( iveg == Nveg) Cv = 1. - veg_con[0].Cv_sum;
	      else
		Cv = lake_con->Cl[0];
	      for ( dist = 0; dist < Ndist; dist++ ) {// wet/dry
		if(dist==0) 
		  tmp_mu = prcp->mu[iveg];
		else 
		  tmp_mu = 1. - prcp->mu[iveg];
#if SPATIAL_FROST
		tmp_ice = 0;
		for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {//frost area
		  tmp_ice += (cell[dist][iveg][band].layer[lidx].ice[frost_area]
			      * soil_con->frost_fract[frost_area]);
		  ice_layer = cell[dist][iveg][band].layer[lidx].ice[frost_area];
		  if(ice_layer>=max_ice_layer)
		    max_ice_layer = ice_layer;
		} // frost area
#else //SPATIAL_FROST
		tmp_ice = cell[dist][iveg][band].layer[lidx].ice;
		ice_layer = cell[dist][iveg][band].layer[lidx].ice;
		if(ice_layer>=max_ice_layer)
		  max_ice_layer = ice_layer;	
#endif //SPATIAL_FROST
		ave_ice += tmp_ice * Cv * tmp_mu * soil_con->AreaFract[band];
	      }// wet/dry
	    }
	  }//iveg
	}
      } //band
      ave_ice_fract = ave_ice/soil_con->max_moist[lidx];

      /*check to see if threshold is exceeded by average ice/porosity fraction*/
      if(ave_ice_fract <= ICE_AT_SUBSIDENCE) {
	SubsidenceUpdate = 1;

	/*calculate subsidence based on maximum ice volume*/
	/*make slightly larger than what can contain max ice volume*/
	tmp_depth_prior = soil_con->depth[lidx];
	tmp_depth = 1.01*max_ice_layer/1000. + soil_con->depth[lidx]*(1.0 - soil_con->effective_porosity[lidx]);
	if(tmp_depth <= soil_con->min_depth[lidx])
	  tmp_depth = soil_con->min_depth[lidx];

	soil_con->depth[lidx] = (float)(int)(tmp_depth * 1000 + 0.5) / 1000;	
	subsidence[lidx] = 1000.0*(tmp_depth_prior - soil_con->depth[lidx]); //in mm
	total_subsidence += (tmp_depth_prior - soil_con->depth[lidx]);
	if(subsidence[lidx] > 0 ){
#if VERBOSE
	  fprintf(stderr,"Subsidence of %.3f m in layer %d:\n",subsidence[lidx]/1000.,lidx+1);
	  fprintf(stderr,"\t\tOccurred for record=%d: year=%d, month=%d, day=%d, hour=%d\n",rec,dmy[rec].year,dmy[rec].month,dmy[rec].day, dmy[rec].hour);
	  fprintf(stderr,"\t\tDepth of soil layer decreased from %.3f m to %.3f m.\n",tmp_depth_prior,soil_con->depth[lidx]);
#endif	  

	  /*update soil_con properties*/
#if VERBOSE
	  fprintf(stderr,"\t\tEffective porosity decreased from %.2f to %.2f.\n",soil_con->effective_porosity[lidx],1.0-(1.0-soil_con->effective_porosity[lidx])*tmp_depth_prior/soil_con->depth[lidx]);
#endif
	  soil_con->effective_porosity[lidx]=1.0-(1.0-soil_con->effective_porosity[lidx])*tmp_depth_prior/soil_con->depth[lidx];
	  if(tmp_depth <= soil_con->min_depth[lidx])
	    soil_con->effective_porosity[lidx]=soil_con->porosity[lidx];
#if VERBOSE
	  fprintf(stderr,"\t\tBulk density increased from %.2f kg/m^3 to %.2f kg/m^3.\n",soil_con->bulk_density[lidx],(1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]);
#endif
	  soil_con->bulk_density[lidx] = (1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]; //adjust bulk density
	  total_meltwater += soil_con->max_moist[lidx] - soil_con->depth[lidx] * soil_con->effective_porosity[lidx] * 1000.; //for lake model
	  soil_con->max_moist[lidx] = soil_con->depth[lidx] * soil_con->effective_porosity[lidx] * 1000.;
	  
	}//subsidence occurs
      }//threshold exceeded
      
    }//excess ice exists
  } //loop for each soil layer
  if(total_subsidence>0){

    /********update remaining soil_con properties**********/
#if VERBOSE
    fprintf(stderr,"Damping depth decreased from %.2f m to %.2f m.\n",soil_con->dp,soil_con->dp-total_subsidence);
#endif
    soil_con->dp -= total_subsidence;  //adjust damping depth

#if VERBOSE
    fprintf(stderr,"More updated parameters in soil_con: max_infil, Wcr, and Wpwp\n");
#endif  
    /* update Maximum Infiltration for Upper Layers */
    if(options.Nlayer==2)
      soil_con->max_infil = (1.0+soil_con->b_infilt)*soil_con->max_moist[0];
    else
      soil_con->max_infil = (1.0+soil_con->b_infilt)*(soil_con->max_moist[0]+soil_con->max_moist[1]);
    
    /* Soil Layer Critical and Wilting Point Moisture Contents */
    for(lidx=0;lidx<options.Nlayer;lidx++) {//soil layer
      soil_con->Wcr[lidx]  = soil_con->Wcr_FRACT[lidx] * soil_con->max_moist[lidx];
      soil_con->Wpwp[lidx] = soil_con->Wpwp_FRACT[lidx] * soil_con->max_moist[lidx];
      if(soil_con->Wpwp[lidx] > soil_con->Wcr[lidx]) {
	sprintf(ErrStr,"Updated wilting point moisture (%f mm) is greater than updated critical point moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
		soil_con->Wpwp[lidx], soil_con->Wcr[lidx], lidx);
	nrerror(ErrStr);
      }
      if(soil_con->Wpwp[lidx] < soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000.) {
	sprintf(ErrStr,"Updated wilting point moisture (%f mm) is less than updated residual moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
		soil_con->Wpwp[lidx], soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000., lidx);
	nrerror(ErrStr);
      }
    }      

    /* If BASEFLOW = NIJSSEN2001 then convert ARNO baseflow
       parameters d1, d2, d3, and d4 to Ds, Dsmax, Ws, and c */
    if(options.BASEFLOW == NIJSSEN2001) {
      lidx = options.Nlayer-1;
      soil_con->Dsmax = soil_con->Dsmax_orig * 
	pow((double)(1./(soil_con->max_moist[lidx]-soil_con->Ws_orig)), -soil_con->c) +
	soil_con->Ds_orig * soil_con->max_moist[lidx];
      soil_con->Ds = soil_con->Ds_orig * soil_con->Ws_orig / soil_con->Dsmax_orig;
      soil_con->Ws = soil_con->Ws_orig/soil_con->max_moist[lidx];
#if VERBOSE
      fprintf(stderr,"More updated parameters in soil_con: Dsmax, Ds, Ws\n");
#endif    
    }
    
    /*********** update root fractions ***************/
    fprintf(stderr,"Updated parameter in veg_con: root\n");
    calc_root_fractions(veg_con, soil_con);

    /**********redistribute soil moisture (call runoff function)*************/
    for(iveg = 0; iveg <= Nveg; iveg++){
      if ((iveg <  Nveg && veg_con[iveg].Cv  > 0.) || 
	  (iveg == Nveg && veg_con[0].Cv_sum < 1.)) {
	for ( band = 0; band < Nbands; band++ ) {
	  if( soil_con->AreaFract[band] > 0 ) { 

	    //set inflow for runoff call 
	    ppt[WET]=cell[WET][iveg][band].inflow;
	    ppt[DRY]=cell[DRY][iveg][band].inflow;
	    
	    /* If subsidence occurs, recalculate runoff, baseflow, and soil moisture
	       using soil moisture values from previous time-step; i.e.
	       as if prior runoff call did not occur.*/
	    for ( dist = 0; dist < Ndist; dist++ ) {
	      for(lidx=0;lidx<options.Nlayer;lidx++) {			
		cell[dist][iveg][band].layer[lidx].moist = moist_prior[dist][iveg][band][lidx];
		cell[dist][iveg][band].layer[lidx].evap = evap_prior[dist][iveg][band][lidx];
	      }
	    }

	    ErrorFlag = runoff(cell[WET][iveg][band].layer, cell[DRY][iveg][band].layer, &(energy[iveg][band]), 
			       soil_con, &(cell[WET][iveg][band].runoff), &(cell[DRY][iveg][band].runoff), 
			       &(cell[WET][iveg][band].baseflow), &(cell[DRY][iveg][band].baseflow), ppt, 
			       SubsidenceUpdate,
#if SPATIAL_FROST
			       soil_con->frost_fract,
#endif // SPATIAL_FROST
			       prcp->mu[iveg], gp->dt, options.Nnode, band, rec, iveg);
	    if ( ErrorFlag == ERROR ) return ( ERROR );

	  }
	}//band
      }
    }//veg
  
    /**********interpolate nodal temperatures to new depths and recalculate thermal properties***********/
    ErrorFlag = update_thermal_nodes(prcp, Nveg, options.Nnode, Ndist, soil_con, veg_con);
    if ( ErrorFlag == ERROR ) return ( ERROR );

  }//subsidence occurs

  /********************************************
    Save subsidence for output
  ********************************************/
  for(lidx=0;lidx<options.Nlayer;lidx++)
    soil_con->subsidence[lidx] = subsidence[lidx];
  
#endif //EXCESS_ICE

  /****************************
     Run Lake Model           
  ****************************/

  /** Compute total runoff and baseflow for all vegetation types
      within each snowband. **/
  if ( options.LAKES && lake_con->Cl[0] > 0 ) {

    sum_runoff = sum_baseflow = 0;
	
    // Loop through snow elevation bands
    for ( band = 0; band < Nbands; band++ ) {
      if ( soil_con->AreaFract[band] > 0 ) {
	
	// Loop through all vegetation types
	for ( iveg = 0; iveg <= Nveg; iveg++ ) {
	  
	  /** Solve Veg Type only if Coverage Greater than 0% **/
	  if ((iveg <  Nveg && veg_con[iveg].Cv  > 0.) || 
	      (iveg == Nveg && veg_con[0].Cv_sum < 1.)) {

	    // Loop through distributed precipitation fractions
	    for ( dist = 0; dist < 2; dist++ ) {
	      
	      if ( dist == 0 ) 
		tmp_mu = prcp->mu[iveg]; 
	      else 
		tmp_mu = 1. - prcp->mu[iveg]; 
	      
	      if ( iveg < Nveg ) {
		Cv = veg_con[iveg].Cv;
		sum_runoff += ( cell[dist][iveg][band].runoff * tmp_mu 
				* Cv * soil_con->AreaFract[band] );
		sum_baseflow += ( cell[dist][iveg][band].baseflow * tmp_mu 
				  * Cv * soil_con->AreaFract[band] );
	      }
	      else {
		Cv = 1. - veg_con[0].Cv_sum;
		sum_runoff += ( cell[dist][iveg][band].runoff * tmp_mu 
				* Cv * soil_con->AreaFract[band] );
		sum_baseflow += ( cell[dist][iveg][band].baseflow * tmp_mu 
				  * Cv * soil_con->AreaFract[band] ); 
	      }
	    }
	  }
	}
      }
    }

    /** Run lake model **/
    lake_var->runoff_in   = sum_runoff;
    lake_var->baseflow_in = sum_baseflow;
    rainonly = calc_rainonly(atmos->air_temp[NR], atmos->prec[NR], 
			     gp->MAX_SNOW_TEMP, gp->MIN_RAIN_TEMP, 1);
    if ( (int)rainonly == ERROR ) {
      return( ERROR );
    }
    lake_prec = ( gauge_correction[SNOW] * (atmos->prec[NR] - rainonly) 
		  + gauge_correction[RAIN] * rainonly );
    // atmos->out_prec += lake_prec * lake_con->Cl[0];
    
    ErrorFlag = lakemain(atmos, *lake_con, gauge_correction[SNOW] * (atmos->prec[NR] - rainonly),
			 gauge_correction[RAIN] * rainonly, soil_con,
#if EXCESS_ICE
			 SubsidenceUpdate, total_meltwater,
#endif
			 (float)gp->dt, prcp, NR, rec, 
			 gp->wind_h, gp, dmy, Nveg+1, 0);
    if ( ErrorFlag == ERROR ) return ( ERROR );
    
#if LINK_DEBUG
    if ( debug.PRT_LAKE ) { 
      if ( rec == 0 ) {
	// print file header
	fprintf(debug.fg_lake,"Date,Rec,AeroResist,BasinflowIn,BaseflowOut");
	for ( i = 0; i < MAX_LAKE_NODES; i++ )
	  fprintf(debug.fg_lake, ",Density%i", i);
	fprintf(debug.fg_lake,",Evap,IceFract,IceHeight,LakeDepth,RunoffIn,RunoffOut,SurfaceArea,SnowDepth,SnowMelt");
	for ( i = 0; i < MAX_LAKE_NODES; i++ )
	  fprintf(debug.fg_lake, ",Area%i", i);
	fprintf(debug.fg_lake,",SWE");
	for ( i = 0; i < MAX_LAKE_NODES; i++ )
	  fprintf(debug.fg_lake, ",Temp%i", i);
	fprintf(debug.fg_lake,",IceTemp,TpIn,Volume,Nodes,MinMax");
	fprintf(debug.fg_lake,",AlbedoLake,AlbedoOver,AlbedoUnder,AtmosError,AtmosLatent,AtmosLatentSub,AtmosSensible,LongOverIn,LongUnderIn,LongUnderOut,NetLongAtmos,NetLongOver,NetLongUnder,NetShortAtmos,NetShortGrnd,NetShortOver,NetShortUnder");
	fprintf(debug.fg_lake,",ShortOverIn,ShortUnderIn,advection,deltaCC,deltaH,error,fusion,grnd_flux,latent,latent_sub,longwave,melt_energy,out_long_surface,refreeze_energy,sensible,shortwave");
	fprintf(debug.fg_lake,",Qnet,albedo,coldcontent,coverage,density,depth,mass_error,max_swq,melt,pack_temp,pack_water,store_coverage,store_swq,surf_temp,surf_water,swq,swq_slope,vapor_flux,last_snow,store_snow\n");
      }

      // print lake variables
      fprintf(debug.fg_lake, "%i/%i/%i %i:00:00,%i", dmy[rec].month, dmy[rec].day, dmy[rec].year, dmy[rec].hour, rec);
      fprintf(debug.fg_lake, ",%f", lake_var->aero_resist);
      fprintf(debug.fg_lake, ",%f", lake_var->baseflow_in);
      fprintf(debug.fg_lake, ",%f", lake_var->baseflow_out);
      for ( i = 0; i < MAX_LAKE_NODES; i++ )
	fprintf(debug.fg_lake, ",%f", lake_var->density[i]);
      fprintf(debug.fg_lake, ",%f", lake_var->evapw);
      fprintf(debug.fg_lake, ",%f", lake_var->fraci);
      fprintf(debug.fg_lake, ",%f", lake_var->hice);
      fprintf(debug.fg_lake, ",%f", lake_var->ldepth);
      fprintf(debug.fg_lake, ",%f", lake_var->runoff_in);
      fprintf(debug.fg_lake, ",%f", lake_var->runoff_out);
      fprintf(debug.fg_lake, ",%f", lake_var->sarea);
      fprintf(debug.fg_lake, ",%f", lake_var->sdepth);
      fprintf(debug.fg_lake, ",%f", lake_var->snowmlt);
      for ( i = 0; i < MAX_LAKE_NODES; i++ )
	fprintf(debug.fg_lake, ",%f", lake_var->surface[i]);
      fprintf(debug.fg_lake, ",%f", lake_var->swe);
      for ( i = 0; i < MAX_LAKE_NODES; i++ )
	fprintf(debug.fg_lake, ",%f", lake_var->temp[i]);
      fprintf(debug.fg_lake, ",%f", lake_var->tempi);
      fprintf(debug.fg_lake, ",%f", lake_var->tp_in);
      fprintf(debug.fg_lake, ",%f", lake_var->volume);
      fprintf(debug.fg_lake, ",%i", lake_var->activenod);
      fprintf(debug.fg_lake, ",%i", lake_var->mixmax);
      
      // print lake energy variables
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AlbedoLake);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AlbedoOver);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AlbedoUnder);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AtmosError);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AtmosLatent);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AtmosLatentSub);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].AtmosSensible);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].LongOverIn);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].LongUnderIn);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].LongUnderOut);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetLongAtmos);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetLongOver);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetLongUnder);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetShortAtmos);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetShortGrnd);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetShortOver);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].NetShortUnder);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].ShortOverIn);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].ShortUnderIn);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].advection);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].deltaCC);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].deltaH);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].error);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].fusion);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].grnd_flux);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].latent);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].latent_sub);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].longwave);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].melt_energy);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].out_long_surface);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].refreeze_energy);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].sensible);
      fprintf(debug.fg_lake, ",%f", energy[Nveg+1][0].shortwave);
      
      // print lake snow variables
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].Qnet);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].albedo);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].coldcontent);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].coverage);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].density);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].depth);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].mass_error);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].max_swq);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].melt);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].pack_temp);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].pack_water);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].store_coverage);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].store_swq);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].surf_temp);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].surf_water);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].swq);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].swq_slope);
      fprintf(debug.fg_lake, ",%f", snow[Nveg+1][0].vapor_flux);
      fprintf(debug.fg_lake, ",%i", snow[Nveg+1][0].last_snow);
      fprintf(debug.fg_lake, ",%i\n", snow[Nveg+1][0].store_snow);

    }
#endif // LINK_DEBUG

  }

  return (0);
}

