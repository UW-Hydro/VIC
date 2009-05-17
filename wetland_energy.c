#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

int wetland_energy(int                  rec,
		   atmos_data_struct   *atmos,
		   dist_prcp_struct    *prcp,
		   dmy_struct          *dmy,
		   global_param_struct *gp,
		   soil_con_struct     *soil_con,
#if EXCESS_ICE
		   int                  SubsidenceUpdate,
#endif
		   int                  iveg,
		   int                  band,
		   double               lake_frac,
		   lake_con_struct      lake_con,
		   veg_con_struct      *veg_con)

/**********************************************************************
	wetland_energy	Laura Bowling		May 12, 2002

  Modifications:
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    used in flux calculations.					TJB
  2006-Sep-23 Implemented flexible output configuration; atmos->out_rain
	      and atmos->out_snow must be tracked now.			TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-03 Changed routine type to int so that it can return numeric
	      ERROR values.						KAC via GCT
  2007-Aug-21 Added features for EXCESS_ICE option.			JCA
  2007-Nov-06 Modified so that wetland vegetation defaults to the first
	      vegetation type in the parameter file, rather than
	      hard-wired.						LCB via TJB
  2008-Mar-01 Reinserted logic for QUICK_FS in calls to
	      distribute_node_moisture_properties().			TJB
  2008-Jun-16 Now runs even when wetland fraction is 0.			LCB via TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2009-Mar-15 Added code to assign values to aero_resist_used.		TJB
  2009-May-17 Added asat to cell_data.					TJB
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
  int                    Ndist;
  int                    dist;
  int                    Nveg;
  int                    veg_class;
  int                    Nbands;
  int                    hour;
  int                    ErrorFlag;
  int                    lidx;
  double                 tmp_moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
#if EXCESS_ICE
  double                 ppt[2]; 
  double                 dummy[MAX_LAYERS]; //dummy array, not needed here
#endif
  double                 rad_atten;
  double                 LAI;
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
  double                 step_rad;
  double                 step_net_short;
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
  float                  root[MAX_LAYERS];
  double                 gauge_correction[2];
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
  float                  lag_one;
  float                  sigma_slope;
  float                  fetch;
  double MOIST1, MOIST2, EVAP, DSWE, ERR, WDEW1, WDEW2;

  /* set local pointers */
  cell    = prcp->cell;
  energy  = prcp->energy;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;

  /* set variables for distributed precipitation */
  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;

  /* No snow bands in lake model. */
  Nbands = 1;

  /* Set number of vegetation types */
  Nveg      =  iveg-1;

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

  MOIST1=0.0;
  for(i=0; i<options.Nlayer; i++) {
    MOIST1 += cell[WET][iveg][band].layer[i].moist;
  }
  DSWE = snow[iveg][band].swq;
  WDEW1 = veg_var[WET][iveg][band].Wdew;

//  /** Solve Wetland only if Lake Coverage Less than 100% **/
//  if (fabs(lake_frac - 1.0) > SMALL) {

    /* fraction of lake area. */
    Cv = 1. - lake_frac;

    /**************************************************
        Initialize Model Parameters
    **************************************************/
      
    /* Initialize energy balance variables */
    energy[iveg][band].shortwave = 0;
    energy[iveg][band].longwave  = 0.;

    /* Initialize snow variables */
    snow[iveg][band].vapor_flux        = 0.;
    snow[iveg][band].canopy_vapor_flux = 0.;
    snow_inflow[band]                  = 0.;
    Melt[band*2]                       = 0.;
  
    /* Initialize precipitation storage */
    for ( j = 0; j < 2; j++ ) {
      out_prec[j] = 0;
      out_rain[j] = 0;
      out_snow[j] = 0;
    }

    /** Define vegetation class number **/
    veg_class = veg_con[0].veg_class;
    wind_h = veg_lib[veg_class].wind_h;

    /** Compute Surface Attenuation due to Vegetation Coverage **/
    surf_atten = exp(-veg_lib[veg_class].rad_atten
                     * veg_lib[veg_class].LAI[dmy[rec].month-1]);

    /* Initialize soil thermal properties for the top two layers */
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      prepare_full_energy(iveg, Nveg, options.Nnode, prcp, 
			  soil_con, &moist, &ice0);
    }

    /** Compute Bare Soil (free of snow) Albedo **/
    bare_albedo = veg_lib[veg_class].albedo[dmy[rec].month-1];       
  
    /*************************************
	Compute the aerodynamic resistance 
    *************************************/

    /* Set surface descriptive variables */
    displacement[0] = veg_lib[veg_class].displacement[dmy[rec].month-1];
    roughness[0]    = veg_lib[veg_class].roughness[dmy[rec].month-1];
    overstory       = veg_lib[veg_class].overstory;
    if ( roughness[0] == 0 ) roughness[0] = soil_con->rough;

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
    if ( ( ErrorFlag = CalcAerodynamic(overstory, height, veg_lib[veg_class].trunk_ratio,
                    soil_con->snow_rough, soil_con->rough,
                    veg_lib[veg_class].wind_atten, gp->wind_h,
                    cell[WET][iveg][0].aero_resist, tmp_wind,
                    displacement, ref_height, roughness,
                    Nveg, iveg) ) == ERROR)
      return(ERROR);
    cell[WET][iveg][0].aero_resist_used[0] = cell[WET][iveg][0].aero_resist[0];
    cell[WET][iveg][0].aero_resist_used[1] = cell[WET][iveg][0].aero_resist[1];

    /******************************
        Solve ground surface fluxes 
    ******************************/
    
    wet_veg_var = &(veg_var[WET][iveg][band]);
    dry_veg_var = &(veg_var[DRY][iveg][band]);
    lag_one     = veg_con[0].lag_one;
    sigma_slope = veg_con[0].sigma_slope;
    fetch       = veg_con[0].fetch;

#if EXCESS_ICE
    if(SubsidenceUpdate == 1 ){ //if subsidence has occurred in this time-step,
                                //call runoff before surface_fluxes, then don't
                                //call runoff again in surface_fluxes.
      //set inflow for runoff call       
      ppt[WET]=cell[WET][iveg][band].inflow;
      ppt[DRY]=cell[DRY][iveg][band].inflow;
      ErrorFlag = runoff(cell[WET][iveg][band].layer, cell[DRY][iveg][band].layer, &(energy[iveg][band]), 
			 soil_con, &(cell[WET][iveg][band].runoff), &(cell[DRY][iveg][band].runoff), 
			 &(cell[WET][iveg][band].baseflow), &(cell[DRY][iveg][band].baseflow),
			 &(cell[WET][iveg][band].asat), &(cell[DRY][iveg][band].asat), ppt, 
			 SubsidenceUpdate,
#if SPATIAL_FROST
			 soil_con->frost_fract,
#endif // SPATIAL_FROST
			 prcp->mu[iveg], gp->dt, options.Nnode, band, rec, iveg);
      if ( ErrorFlag == ERROR ) return ( ERROR );

      SubsidenceUpdate = 2;
    }
#endif //EXCESS_ICE
    
    for ( lidx = 0; lidx < options.Nlayer; lidx++ ) 
      tmp_moist[iveg][band][lidx] = cell[0][iveg][band].layer[lidx].moist;
    ErrorFlag = distribute_node_moisture_properties(energy[iveg][band].moist,
						    energy[iveg][band].ice,
						    energy[iveg][band].kappa_node,
						    energy[iveg][band].Cs_node,
						    soil_con->Zsum_node,
						    energy[iveg][band].T,
						    soil_con->max_moist_node,
#if QUICK_FS
						    soil_con->ufwc_table_node,
#else
						    soil_con->expt_node,
						    soil_con->bubble_node,
#endif // QUICK_FS
#if EXCESS_ICE
						    soil_con->porosity_node,
						    soil_con->effective_porosity_node,
#endif // EXCESS_ICE
						    tmp_moist[iveg][band],
						    soil_con->depth,
						    soil_con->soil_density,
						    soil_con->bulk_density,
						    soil_con->quartz,
						    options.Nnode, options.Nlayer,
						    soil_con->FS_ACTIVE);
#if EXCESS_ICE
    if ( ErrorFlag == ERROR ) {//this means there has been a new wetland area
			       //since subsidence, so first redistribute soil moisture
      
      SubsidenceUpdate = 1;
      //set inflow for runoff call       
      ppt[WET]=cell[WET][iveg][band].inflow;
      ppt[DRY]=cell[DRY][iveg][band].inflow;
      ErrorFlag = runoff(cell[WET][iveg][band].layer, cell[DRY][iveg][band].layer, &(energy[iveg][band]), 
			 soil_con, &(cell[WET][iveg][band].runoff), &(cell[DRY][iveg][band].runoff), 
			 &(cell[WET][iveg][band].baseflow), &(cell[DRY][iveg][band].baseflow), ppt, 
			 SubsidenceUpdate,
#if SPATIAL_FROST
			 soil_con->frost_fract,
#endif // SPATIAL_FROST
			 prcp->mu[iveg], gp->dt, options.Nnode, band, rec, iveg);
      if ( ErrorFlag == ERROR ) return ( ERROR );
    
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) 
	tmp_moist[iveg][band][lidx] = cell[0][iveg][band].layer[lidx].moist;
// NOTE: this is the version of distribute_node_moisture_properties for EXCESS_ICE = TRUE
      ErrorFlag = distribute_node_moisture_properties(energy[iveg][band].moist,
						      energy[iveg][band].ice,
						      energy[iveg][band].kappa_node,
						      energy[iveg][band].Cs_node,
						      soil_con->Zsum_node,
						      energy[iveg][band].T,
						      soil_con->max_moist_node,
#if QUICK_FS
						      soil_con->ufwc_table_node,
#else
						      soil_con->expt_node,
						      soil_con->bubble_node,
#endif // QUICK_FS
						      soil_con->porosity_node,
						      soil_con->effective_porosity_node,
						      tmp_moist[iveg][band],
						      soil_con->depth,
						      soil_con->soil_density,
						      soil_con->bulk_density,
						      soil_con->quartz,
						      options.Nnode, options.Nlayer,
						      soil_con->FS_ACTIVE);
      if ( ErrorFlag == ERROR ) return ( ERROR );
      
      SubsidenceUpdate = 2;
    }
#else //!EXCESS_ICE
    if ( ErrorFlag == ERROR ) return ( ERROR );
#endif //EXCESS_ICE

    
#if EXCESS_ICE
    for ( lidx = 0; lidx < options.Nlayer; lidx++ ) 
      dummy[lidx]=0;
#endif    

    ErrorFlag = surface_fluxes(overstory, bare_albedo, height, ice0, moist, 
#if EXCESS_ICE
			       SubsidenceUpdate, dummy, dummy,
#endif
                               prcp->mu[iveg], surf_atten, &(Melt[band*2]), &Le,
                               cell[WET][iveg][0].aero_resist,cell[WET][iveg][0].aero_resist_used,
                               &(cell[DRY][iveg][band].baseflow),
                               &(cell[WET][iveg][band].baseflow),
                               &(cell[DRY][iveg][band].asat),
                               &(cell[WET][iveg][band].asat), displacement,
                               gauge_correction, &(cell[DRY][iveg][band].inflow),
                               &(cell[WET][iveg][band].inflow), &out_prec[band*2],
                               &out_rain[band*2], &out_snow[band*2],
                               ref_height, roughness,
                               &(cell[DRY][iveg][band].runoff),
                               &(cell[WET][iveg][band].runoff), &snow_inflow[band],
                               tmp_wind, veg_con[0].root, Nbands, Ndist,
                               options.Nlayer, Nveg, band, dp, iveg, rec, veg_class,
                               atmos, dmy, &(energy[iveg][band]), gp,
                               cell[DRY][iveg][band].layer,
                               cell[WET][iveg][band].layer, &(snow[iveg][band]),
                               soil_con, dry_veg_var, wet_veg_var, lag_one, sigma_slope, fetch);

    if ( ErrorFlag == ERROR )
      // Failed in surface_fluxes, return error flag
      return( ErrorFlag );
    
    EVAP = MOIST2=0.0;
    for(i=0; i<options.Nlayer; i++) {
      MOIST2 += cell[WET][iveg][band].layer[i].moist;
      EVAP += cell[WET][iveg][band].layer[i].evap;
    }
    EVAP += veg_var[WET][iveg][band].canopyevap;
    EVAP += snow[iveg][band].vapor_flux * 1000.;
    EVAP += snow[iveg][band].canopy_vapor_flux * 1000.;
    WDEW2 = veg_var[WET][iveg][band].Wdew;
    DSWE -= snow[iveg][band].swq;
    ERR = out_prec[band*2] - (cell[WET][iveg][band].baseflow+cell[WET][iveg][band].runoff+EVAP)-(-1000*DSWE+(MOIST2-MOIST1)+(WDEW2-WDEW1));
//    if(fabs(ERR) > SMALL) {
//      fprintf(stderr, "prec=%f, b=%f, R=%f ET=%f DeltaM=%f DSWE=%f Ddew=%f ERR=%e\n",out_prec[band*2],cell[WET][iveg][band].baseflow,cell[WET][iveg][band].runoff, EVAP, MOIST2-MOIST1, -1000*DSWE, WDEW2-WDEW1, ERR);
//    }

    atmos->out_prec += out_prec[band*2] * Cv * lake_con.Cl[0];
    atmos->out_rain += out_rain[band*2] * Cv * lake_con.Cl[0];
    atmos->out_snow += out_snow[band*2] * Cv * lake_con.Cl[0];

//  } /** end if fabs(lake_frac-1.0) > SMALL **/

  return (0);  

}
