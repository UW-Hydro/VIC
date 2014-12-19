#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  full_energy(int                  gridcell,
                 int                  rec,
                 atmos_data_struct   *atmos,
                 all_vars_struct     *all_vars,
                 dmy_struct          *dmy,
                 global_param_struct *gp,
		 lake_con_struct     *lake_con,
                 soil_con_struct     *soil_con,
                 veg_con_struct      *veg_con,
                 veg_hist_struct    **veg_hist)
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
	      of soil wetness and root zone soil moisture happens here.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.		GCT/KAC
  2007-May-01 Added case of SPATIAL_FROST = TRUE in modifications
	      from 2006-Sep-23.							GCT
  2007-Aug-10 Added features for EXCESS_ICE option.
	      Including calculating subsidence for each layer and
	      updating soil depth, effective porosity,
	      bulk density, and soil moisture and fluxes by calling
	      runoff function if subsidence occurs.				JCA
  2007-Sep-07 No longer resets ice content to previous time-step ice content if
	      subsidence has occurred.						JCA
  2007-Sep-19 Added MAX_SUBSIDENCE parameter to EXCESS_ICE option.		JCA
  2007-Sep-19 Fixed bug in subsidence calculation.				JCA
  2007-Nov-06 Added veg_con to parameter list of lakemain().  Replaced
	      lake.fraci with lake.areai.					LCB via TJB
  2008-Jan-23 Changed ice0 from a scalar to an array.  Previously,
	      when options.SNOW_BAND > 1, the value of ice0 computed
	      for earlier bands was always overwritten by the value
	      of ice0 computed for the final band (even if the final
	      band had 0 area).							JS via TJB
  2008-May-05 Changed moist from a scalar to an array (moist0).  Previously,
	      when options.SNOW_BAND > 1, the value of moist computed
	      for earlier bands was always overwritten by the value
	      of moist computed for the final band (even if the final
	      band had 0 area).							KAC via TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.					TJB
  2009-May-17 Added asat to cell_data.						TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.						TJB
  2009-Jun-09 Modified to compute aero_resist for all potential evap
	      landcover types.							TJB
  2009-Jun-09 Cell_data structure now only stores final aero_resist
	      values (called "aero_resist").  Preliminary uncorrected
	      aerodynamic resistances for current vegetation and various
	      reference land cover types for use in potential evap
	      calculations is stored in temporary array aero_resist.		TJB
  2009-Jun-26 Simplified argument list of runoff() by passing all cell_data
	      variables via a single reference to the cell data structure.	TJB
  2009-Jul-22 Fixed error in assignment of cell.aero_resist.			TJB
  2009-Jul-31 Wetland portion of lake/wetland tile is now processed in
	      full_energy() instead of wetland_energy().  Lake funcions are
	      now called directly from full_energy instead of lakemain().	TJB
  2009-Sep-28 Moved lake_snow and lake_energy into lake_var structure.
	      Removed calls to initialize_prcp and update_prcp.			TJB
  2009-Sep-30 Miscellaneous fixes for lake model.				TJB
  2009-Oct-05 Miscellaneous fixes for lake model, including updating/
	      rescaling of lake and wetland storages and fluxes to account
	      for changes in lake area.						TJB
  2009-Nov-09 Changed definition of lake->sarea to include ice extent; other
	      changes to handle case when lake fraction goes to 0.		LCB via TJB
  2010-Mar-31 Added runoff_in.							TJB
  2010-Sep-24 Changed atmos.runoff_in to atmos.channel_in.  Added
	      lake_var.channel_in to store it.					TJB
  2010-Nov-02 Changed units of lake_var moisture fluxes to volume (m3).		TJB
  2010-Nov-26 Changed argument list of water_balance().				TJB
  2011-May-31 Prepare_full_energy() is now always called.			TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2012-Jan-01 Modified condition for determining whether to simulate lakes
	      to check whether lake_idx >= 0.					TJB
  2012-Jan-16 Removed LINK_DEBUG code						BN
  2012-Aug-28 Added accumulation of rain and snow over lake to grid cell
	      totals.								TJB
  2013-Jul-25 Added photosynthesis terms.					TJB
  2013-Jul-25 Added soil carbon terms.						TJB
  2013-Dec-26 Removed EXCESS_ICE option.					TJB
  2013-Dec-27 Removed (unused) SPATIAL_FROST code.				TJB
  2014-Mar-28 Removed DIST_PRCP option.						TJB
  2014-Apr-25 Added non-climatological veg params.				TJB
  2014-Apr-25 Added partial vegcover fraction.					TJB

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
  char                   overstory;
  int                    i, j, p;
  int                    lidx;
  int                    iveg;
  int                    Nveg;
  int                    veg_class;
  int                    band;
  int                    Nbands;
  int                    ErrorFlag;
  double                 out_prec[2*MAX_BANDS];
  double                 out_rain[2*MAX_BANDS];
  double                 out_snow[2*MAX_BANDS];
  double                 out_short=0;
  double                 dp;
  double                 ice0[MAX_BANDS];
  double                 moist0[MAX_BANDS];
  double                 surf_atten;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 wind_h;
  double                 height;
  double                 displacement[3];
  double                 roughness[3];
  double                 ref_height[3];
  double               **aero_resist;
  double                 Cv;
  double                 Le;
  double                 Melt[2*MAX_BANDS];
  double                 bare_albedo;
  double                 snow_inflow[MAX_BANDS];
  double                 rainonly;
  double                 sum_runoff;
  double                 sum_baseflow;
  double                 tmp_wind[3];
  double                 gauge_correction[2];
  float 	         lag_one;
  float 	         sigma_slope;
  float  	         fetch;
  int                    pet_veg_class;
  double                 lakefrac;
  double                 fraci;
  double                 wetland_runoff;
  double                 wetland_baseflow;
  double                 oldsnow;
  double                 snowprec;
  double                 rainprec;
  int                    cidx;
  lake_var_struct       *lake_var;
  cell_data_struct     **cell;
  veg_var_struct       **veg_var;
  energy_bal_struct    **energy;
  snow_data_struct     **snow;

  /* Allocate aero_resist array */
  aero_resist = (double**)calloc(N_PET_TYPES+1,sizeof(double*));
  for (p=0; p<N_PET_TYPES+1; p++) {
    aero_resist[p] = (double*)calloc(3,sizeof(double));
  }

  /* set local pointers */
  cell    = all_vars->cell;
  energy  = all_vars->energy;
  lake_var = &all_vars->lake_var;
  snow    = all_vars->snow;
  veg_var = all_vars->veg_var;

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

  /* Assign current veg albedo and LAI */
  if (rec >= 0) {
    // Loop over vegetated tiles
    for(iveg = 0; iveg < Nveg; iveg++){
      veg_class = veg_con[iveg].veg_class;
      if (veg_hist[rec][iveg].vegcover[0] < MIN_VEGCOVER)
        veg_hist[rec][iveg].vegcover[0] = MIN_VEGCOVER;
      for ( band = 0; band < Nbands; band++ ) {
        veg_var[iveg][band].vegcover = veg_hist[rec][iveg].vegcover[0];
        veg_var[iveg][band].albedo = veg_hist[rec][iveg].albedo[0];
        veg_var[iveg][band].LAI = veg_hist[rec][iveg].LAI[0];
        // Convert LAI from global to local
        veg_var[iveg][band].LAI /= veg_var[iveg][band].vegcover;
        veg_var[iveg][band].Wdew /= veg_var[iveg][band].vegcover;
        veg_var[iveg][band].Wdmax = veg_var[iveg][band].LAI*LAI_WATER_FACTOR;
        snow[iveg][band].snow_canopy /= veg_var[iveg][band].vegcover;
      }
    }
  }

  /**************************************************
    Solve Energy and/or Water Balance for Each
    Vegetation Type
  **************************************************/
  for(iveg = 0; iveg <= Nveg; iveg++){

    /** Solve Veg Type only if Coverage Greater than 0% **/
    if (veg_con[iveg].Cv > 0.0) {
      Cv = veg_con[iveg].Cv;
      Nbands = options.SNOW_BAND;

      /** Lake-specific processing **/
      if (veg_con[iveg].LAKE) {

        /* Update areai to equal new ice area from previous time step. */
        lake_var->areai = lake_var->new_ice_area;

        /* Compute lake fraction and ice-covered fraction */
        if (lake_var->areai < 0) lake_var->areai = 0;
	if (lake_var->sarea > 0) {
	  fraci = lake_var->areai/lake_var->sarea;
	  if(fraci > 1.0) fraci = 1.0;
	}
	else
	  fraci = 0.0;
	lakefrac = lake_var->sarea/lake_con->basin[0];

        Nbands = 1;
        Cv *= (1-lakefrac);

        if (Cv == 0)
          continue;

      }

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
      veg_class = veg_con[iveg].veg_class;

      /** Initialize other veg vars **/
      if (iveg < Nveg) {
	for(band=0; band<Nbands; band++) {
          veg_var[iveg][band].rc = HUGE_RESIST;
        }
      }

      /** Assign wind_h **/
      /** Note: this is ignored below **/
      wind_h = veg_lib[veg_class].wind_h;

      /** Compute Surface Attenuation due to Vegetation Coverage **/
      surf_atten = (1-veg_var[iveg][0].vegcover)*1.0
                   + veg_var[iveg][0].vegcover
                   * exp(-veg_lib[veg_class].rad_atten * veg_var[iveg][0].LAI);

      /* Initialize soil thermal properties for the top two layers */
      prepare_full_energy(iveg, Nveg, options.Nnode, all_vars, soil_con, moist0, ice0);

      /** Compute Bare (free of snow) Albedo **/
      if (iveg!=Nveg){
        bare_albedo = veg_var[iveg][0].albedo;
      }
      else {
        bare_albedo = BARE_SOIL_ALBEDO;
      }

      /*************************************
	Compute the aerodynamic resistance 
	for current veg cover and various
	types of potential evap
      *************************************/

      /* Loop over types of potential evap, plus current veg */
      /* Current veg will be last */
      for (p=0; p<N_PET_TYPES+1; p++) {

        /* Initialize wind speeds */
        tmp_wind[0] = atmos->wind[NR];
        tmp_wind[1] = -999.;
        tmp_wind[2] = -999.;
 
        /* Set surface descriptive variables */
        if (p < N_PET_TYPES_NON_NAT) {
	  pet_veg_class = veg_lib[0].NVegLibTypes+p;
        }
        else {
	  pet_veg_class = veg_class;
        }
        displacement[0] = veg_lib[pet_veg_class].displacement[dmy[rec].month-1];
        roughness[0]    = veg_lib[pet_veg_class].roughness[dmy[rec].month-1];
        overstory       = veg_lib[pet_veg_class].overstory;
	if (p >= N_PET_TYPES_NON_NAT)
          if ( roughness[0] == 0 ) roughness[0] = soil_con->rough;

        /* Estimate vegetation height */
        height = calc_veg_height(displacement[0]);

        /* Estimate reference height */
        if(displacement[0] < wind_h) 
          ref_height[0] = wind_h;
        else 
          ref_height[0] = displacement[0] + wind_h + roughness[0];

        /* Compute aerodynamic resistance over various surface types */
        /* Do this not only for current veg but also all types of PET */
        ErrorFlag = CalcAerodynamic(overstory, height,
		        veg_lib[pet_veg_class].trunk_ratio, 
                        soil_con->snow_rough, soil_con->rough, 
		        veg_lib[pet_veg_class].wind_atten,
			aero_resist[p], tmp_wind,
		        displacement, ref_height,
		        roughness);
        if ( ErrorFlag == ERROR ) return ( ERROR );  

      }

      /* Initialize final aerodynamic resistance values */
      for ( band = 0; band < Nbands; band++ ) {
        if( soil_con->AreaFract[band] > 0 ) {
          cell[iveg][band].aero_resist[0] = aero_resist[N_PET_TYPES][0];
          cell[iveg][band].aero_resist[1] = aero_resist[N_PET_TYPES][1];
        }
      }

      /******************************
        Compute nitrogen scaling factors and initialize other veg vars
      ******************************/
      if (options.CARBON && iveg < Nveg) {
	for(band=0; band<Nbands; band++) {
          for (cidx=0; cidx<options.Ncanopy; cidx++) {
            veg_var[iveg][band].rsLayer[cidx] = HUGE_RESIST;
          }
          veg_var[iveg][band].aPAR = 0;
          if (dmy->hour == 0) {
            calc_Nscale_factors(veg_lib[veg_class].NscaleFlag,
                                veg_con[iveg].CanopLayerBnd,
                                veg_lib[veg_class].LAI[dmy[rec].month-1],
                                soil_con->lat,
                                soil_con->lng,
                                soil_con->time_zone_lng,
                                dmy[rec],
                                veg_var[iveg][band].NscaleFactor);
          }
          if (dmy[rec].month == 1 && dmy[rec].day == 1) {
            veg_var[iveg][band].AnnualNPPPrev = veg_var[iveg][band].AnnualNPP;
            veg_var[iveg][band].AnnualNPP = 0;
          }
        }
      }

      /******************************
        Solve ground surface fluxes 
      ******************************/
  
      for ( band = 0; band < Nbands; band++ ) {
	if( soil_con->AreaFract[band] > 0 ) {

	  lag_one     = veg_con[iveg].lag_one;
	  sigma_slope = veg_con[iveg].sigma_slope;
	  fetch       = veg_con[iveg].fetch;

	  /* Initialize pot_evap */
	  for (p=0; p<N_PET_TYPES; p++)
	    cell[iveg][band].pot_evap[p] = 0;

	  ErrorFlag = surface_fluxes(overstory, bare_albedo, height, ice0[band], moist0[band], 
				     surf_atten, &(Melt[band*2]), &Le, 
				     aero_resist,
				     displacement, gauge_correction,
				     &out_prec[band*2], 
				     &out_rain[band*2], &out_snow[band*2],
				     ref_height, roughness, 
				     &snow_inflow[band], 
				     tmp_wind, veg_con[iveg].root, Nbands, 
				     options.Nlayer, Nveg, band, dp, iveg, rec, veg_class, 
				     atmos, dmy, &(energy[iveg][band]), gp, 
				     &(cell[iveg][band]),
				     &(snow[iveg][band]), 
				     soil_con, &(veg_var[iveg][band]), 
				     lag_one, sigma_slope, fetch, veg_con[iveg].CanopLayerBnd);
	  
	  if ( ErrorFlag == ERROR ) return ( ERROR );
	  
	  atmos->out_prec += out_prec[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_rain += out_rain[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_snow += out_snow[band*2] * Cv * soil_con->AreaFract[band];

          /********************************************************
            Compute soil wetness and root zone soil moisture
          ********************************************************/
          cell[iveg][band].rootmoist = 0;
          cell[iveg][band].wetness = 0;
          for(lidx=0;lidx<options.Nlayer;lidx++) {
            if (veg_con->root[lidx] > 0) {
              cell[iveg][band].rootmoist += cell[iveg][band].layer[lidx].moist;
            }
	    cell[iveg][band].wetness += (cell[iveg][band].layer[lidx].moist - soil_con->Wpwp[lidx])/(soil_con->porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
          }
          cell[iveg][band].wetness /= options.Nlayer;

	} /** End non-zero area band **/
      } /** End Loop Through Elevation Bands **/

    } /** end non-zero area veg tile **/
  } /** end of vegetation loop **/

  /* Convert LAI back to global */
  if (rec >= 0) {
    for(iveg = 0; iveg < Nveg; iveg++){
      for ( band = 0; band < Nbands; band++ ) {
        veg_var[iveg][band].LAI *= veg_var[iveg][band].vegcover;
        veg_var[iveg][band].Wdmax *= veg_var[iveg][band].vegcover;
      }
    }
  }

  for (p=0; p<N_PET_TYPES+1; p++) {
    free((char *)aero_resist[p]);
  }
  free((char *)aero_resist);

  /****************************
     Run Lake Model           
  ****************************/

  /** Compute total runoff and baseflow for all vegetation types
      within each snowband. **/
  if ( options.LAKES && lake_con->lake_idx >= 0 ) {

    wetland_runoff = wetland_baseflow = 0;
    sum_runoff = sum_baseflow = 0;
	
    // Loop through all vegetation tiles
    for ( iveg = 0; iveg <= Nveg; iveg++ ) {
	  
      /** Solve Veg Tile only if Coverage Greater than 0% **/
      if (veg_con[iveg].Cv  > 0.) {

	Cv = veg_con[iveg].Cv;
        Nbands = options.SNOW_BAND;
        if (veg_con[iveg].LAKE) {
          Cv *= (1-lakefrac);
          Nbands = 1;
        }

        // Loop through snow elevation bands
        for ( band = 0; band < Nbands; band++ ) {
          if ( soil_con->AreaFract[band] > 0 ) {
	
            if (veg_con[iveg].LAKE) {
              wetland_runoff += ( cell[iveg][band].runoff
                          * Cv * soil_con->AreaFract[band] );
              wetland_baseflow += ( cell[iveg][band].baseflow
                          * Cv * soil_con->AreaFract[band] );
              cell[iveg][band].runoff = 0;
              cell[iveg][band].baseflow = 0;
            }
            else {
              sum_runoff += ( cell[iveg][band].runoff
                            * Cv * soil_con->AreaFract[band] );
              sum_baseflow += ( cell[iveg][band].baseflow
                              * Cv * soil_con->AreaFract[band] );
              cell[iveg][band].runoff *= (1-lake_con->rpercent);
              cell[iveg][band].baseflow *= (1-lake_con->rpercent);
            }

	  }
        }
      }
    }

    /** Run lake model **/
    iveg = lake_con->lake_idx;
    band = 0;
    lake_var->runoff_in   = (sum_runoff * lake_con->rpercent + wetland_runoff)*soil_con->cell_area*0.001; // m3
    lake_var->baseflow_in = (sum_baseflow * lake_con->rpercent + wetland_baseflow)*soil_con->cell_area*0.001; // m3
    lake_var->channel_in  = atmos->channel_in[NR]*soil_con->cell_area*0.001; // m3
    lake_var->prec        = atmos->prec[NR]*lake_var->sarea*0.001; // m3
    rainonly = calc_rainonly(atmos->air_temp[NR], atmos->prec[NR], 
			     gp->MAX_SNOW_TEMP, gp->MIN_RAIN_TEMP);
    if ( (int)rainonly == ERROR ) {
      return( ERROR );
    }

    /**********************************************************************
       Solve the energy budget for the lake.
     **********************************************************************/

    oldsnow = lake_var->snow.swq;
    snowprec = gauge_correction[SNOW] * (atmos->prec[NR] - rainonly);
    rainprec = gauge_correction[SNOW] * rainonly;
    atmos->out_prec += (snowprec + rainprec) * lake_con->Cl[0] * lakefrac;
    atmos->out_rain += rainprec * lake_con->Cl[0] * lakefrac;
    atmos->out_snow += snowprec * lake_con->Cl[0] * lakefrac;

    ErrorFlag = solve_lake(snowprec, rainprec, atmos->air_temp[NR],
                           atmos->wind[NR], atmos->vp[NR] / 1000.,
                           atmos->shortwave[NR], atmos->longwave[NR],
                           atmos->vpd[NR] / 1000.,
                           atmos->pressure[NR] / 1000.,
                           atmos->density[NR], lake_var, *lake_con,
                           *soil_con, gp->dt, rec, gp->wind_h, dmy[rec], fraci);
    if ( ErrorFlag == ERROR ) return (ERROR);

    /**********************************************************************
       Solve the water budget for the lake.
     **********************************************************************/

    ErrorFlag = water_balance(lake_var, *lake_con, gp->dt, all_vars, rec, iveg, band, lakefrac, *soil_con, *veg_con);
    if ( ErrorFlag == ERROR ) return (ERROR);

  } // end if (options.LAKES && lake_con->lake_idx >= 0)

  return (0);
}

