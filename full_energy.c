#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void full_energy(char                 NEWCELL,
		 int                  gridcell,
                 int                  rec,
                 atmos_data_struct   *atmos,
                 dist_prcp_struct    *prcp,
                 dmy_struct          *dmy,
                 global_param_struct *gp,
#if LAKE_MODEL
		 lake_con_struct     *lake_con,
#endif // LAKE_MODEL 
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

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  char                   overstory;
  char                   SOLVE_SURF_ENERGY;
  int                    i, j, k;
  int                    lidx;
  int                    Ndist;
  int                    dist;
  int                    iveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands;
  int                    hour;
  double                 out_prec[2*MAX_BANDS];
  double                 out_rain[2*MAX_BANDS];
  double                 out_snow[2*MAX_BANDS];
  double                 tmp_surf_temp;
  double                 last_T1;
  double                 out_short=0;
  double                 inshort;
  double                 inlong;
  double                 dp;
  double                 ice0;
  double                 moist;
  double                 surf_atten;
  double                 Tsurf;
  double                 Tgrnd;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 rainfall[2]; 
  double                 wind_h;
  double                 height;
  double                 displacement[3];
  double                 roughness[3];
  double                 ref_height[3];
  double                 Cv;
  double                 Le;
  double                 Evap;
  double                 Melt[2*MAX_BANDS];
  double                 bare_albedo;
  double                 snow_inflow[MAX_BANDS];
#if LAKE_MODEL
  double                 lake_prec;
  double                 rainonly;
  double                 sum_runoff;
  double                 sum_baseflow;
#endif // LAKE_MODEL
  double                 step_rad;
  double                 step_net_short;
  double                 tmp_aero_resist;
  double                 tmp_throughfall[2][MAX_BANDS];
  double                 tmp_wind[3];
  double                 tmp_melt[MAX_BANDS*2];
  double                 tmp_vapor_flux[MAX_BANDS];
  double                 tmp_canopy_vapor_flux[MAX_BANDS];
  double                 tmp_canopyevap[2][MAX_BANDS];
  double                 tmp_snow_energy;
  double                 tmp_Wdew[2];
  double                 tmp_mu;
  double                 tmp_layerevap[2][MAX_BANDS][MAX_LAYERS];
  double                 tmp_Tmin;
  double                 tmp_total_moist;
  double                 gauge_correction[2];
  float 	         lag_one;
  float 	         sigma_slope;
  float  	         fetch;
#if LAKE_MODEL
  lake_var_struct       *lake_var;
#endif /* LAKE_MODEL */
  veg_var_struct         tmp_veg_var[2];
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

  /* set local pointers */
  cell    = prcp->cell;
  energy  = prcp->energy;
#if LAKE_MODEL
  lake_var = &prcp->lake_var;
#endif /* LAKE_MODEL */
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
      CalcAerodynamic(overstory, height, veg_lib[veg_class].trunk_ratio, 
                      soil_con->snow_rough, soil_con->rough, 
		      veg_lib[veg_class].wind_atten,
		      gp->wind_h, cell[WET][iveg][0].aero_resist, tmp_wind, 
                      displacement, ref_height, roughness, 
                      Nveg, iveg);
  
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

	  surface_fluxes(overstory, bare_albedo, height, ice0, moist, 
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
              tmp_total_moist = cell[dist][iveg][band].layer[lidx].moist + cell[dist][iveg][band].layer[lidx].ice;
              if (veg_con->root[lidx] > 0) {
                cell[dist][iveg][band].rootmoist += tmp_total_moist;
              }
              cell[dist][iveg][band].wetness += (tmp_total_moist - soil_con->Wpwp[lidx])/(soil_con->porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
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

#if LAKE_MODEL
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
    lake_prec = ( gauge_correction[SNOW] * (atmos->prec[NR] - rainonly) 
		  + gauge_correction[RAIN] * rainonly );
    // atmos->out_prec += lake_prec * lake_con->Cl[0];

    lakemain(atmos, *lake_con, gauge_correction[SNOW] * (atmos->prec[NR] - rainonly),gauge_correction[RAIN] * rainonly, 
	     soil_con,(float)gp->dt, prcp, NR, rec, 
    	     gp->wind_h, gp, dmy, Nveg+1, 0);

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

#endif /* LAKE_MODEL */

}

