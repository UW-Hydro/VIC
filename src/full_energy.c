#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  full_energy(int                  gridcell,
                 int                  rec,
                 atmos_data_struct   *atmos,
                 all_vars_struct     *all_vars,
                 all_vars_struct     *all_vars_crop,
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
  int                    target_idx;
  int                    thresh_idx;
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
  double                 moistfract;
  double                 irr_sm_thresh;
  double                 irr_sm_target;
  double                 extract_water;
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
  int month2,layer;
  float fraction_irrveg;
  double irrig_old;
  double irrig_new;
  double noirrig_old;
  double noirrig_new;
  double soilmoist_irr_old;
  double soilmoist_noirr_old;
  double soilmoist_noirr_new;
  double soilmoist_irr_new;
  double irrig_est;
  double irrig;
  double old_crop_frac;
  double new_crop_frac;
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
  int                    cridx;
  double                 tmp_double;
  int                    index;
  int                    frost_area;
  lake_var_struct       *lake_var;
  cell_data_struct     **cell;
  veg_var_struct       **veg_var;
  energy_bal_struct    **energy;
  snow_data_struct     **snow;
  float moistlayer,moisttotal1,moisttotal2,moistnoncrop;
  int localextract,ivegdummy,recdummy; //ingjerd added. at the moment: water extracted, but not applied to fields, so doesn't really do any good!
  double Cvdummy,extractwaterdummy;


  recdummy=153;

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

  // Crop fraction
  if (options.CROPFRAC) { //true or false set in global file, default=false
    if (rec >= 0) {
      for(iveg = 0; iveg < Nveg; iveg++){
        if (veg_con[iveg].crop_frac_active) {
          veg_class = veg_con[iveg].veg_class;
          for ( band = 0; band < Nbands; band++ ) {
            if (rec == 0)
            old_crop_frac = veg_hist[rec][iveg].crop_frac[0];
            else
            old_crop_frac = veg_var[iveg][band].crop_frac;
            if(snow[iveg][band].swq<0.0001 && atmos->air_temp[NR]>7)  //ingjerd added. ok for all bands, veg types etc???? should be same air temp as below! (7 degrees)
            new_crop_frac = veg_hist[rec][iveg].crop_frac[0];
            else new_crop_frac = old_crop_frac;
            moistlayer=moisttotal1=moisttotal2=0.;
            if (new_crop_frac != old_crop_frac) {
              // The portion that grows needs to assimilate state variables
              // from the portion that shrinks; canopy storages need to be rescaled
              if (new_crop_frac > old_crop_frac) { // crop fraction is growing
                for(lidx=0;lidx<options.Nlayer;lidx++) {
                  moistlayer=all_vars_crop->cell[0][band].layer[lidx].moist*(1-old_crop_frac)+all_vars_crop->cell[1][band].layer[lidx].moist*old_crop_frac;
                  moisttotal1+=moistlayer;
                  all_vars_crop->cell[1][band].layer[lidx].moist = (all_vars_crop->cell[1][band].layer[lidx].moist*old_crop_frac+all_vars_crop->cell[0][band].layer[lidx].moist*(new_crop_frac-old_crop_frac))/new_crop_frac;
                  moistlayer=all_vars_crop->cell[0][band].layer[lidx].moist*(1-new_crop_frac)+all_vars_crop->cell[1][band].layer[lidx].moist*new_crop_frac;
                  moisttotal2+=moistlayer;
                  //all_vars_crop->snow[1][band].swq=(all_vars_crop->snow[1][band].swq*old_crop_frac+all_vars_crop->snow[0][band].swq*(new_crop_frac-old_crop_frac))/new_crop_frac;
                  //all_vars_crop->snow[0][band].swq=all_vars_crop->snow[0][band].swq;
                }
              }
              else { // fallow fraction is growing
                for(lidx=0;lidx<options.Nlayer;lidx++) {
                  all_vars_crop->cell[0][band].layer[lidx].moist = (all_vars_crop->cell[0][band].layer[lidx].moist*(1-old_crop_frac)+all_vars_crop->cell[1][band].layer[lidx].moist*(old_crop_frac-new_crop_frac))/(1-new_crop_frac);
                }
              }
              all_vars_crop->veg_var[1][band].Wdew *= old_crop_frac/new_crop_frac;
              all_vars_crop->snow[1][band].snow_canopy *= old_crop_frac/new_crop_frac;
            }
            veg_var[iveg][band].crop_frac = new_crop_frac;
          }
        }
      }
    }
  }

  /* Assign current veg albedo and LAI (standard vic veg ) */
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
        veg_var[iveg][band].Wdmax = veg_var[iveg][band].LAI*LAI_WATER_FACTOR;
        // Handle crop tiles (fallow and crop sub-tiles)
        if (options.CROPFRAC && veg_con[iveg].crop_frac_active) {
          for (cridx=veg_con[iveg].crop_frac_idx; cridx<veg_con[iveg].crop_frac_idx+2; cridx++) {
            if (cridx % 2 == 0) { //fallow part
              all_vars_crop->veg_var[cridx][band].vegcover = MIN_VEGCOVER;
              all_vars_crop->veg_var[cridx][band].albedo = BARE_SOIL_ALBEDO;
              all_vars_crop->veg_var[cridx][band].LAI = 0;
              all_vars_crop->veg_var[cridx][band].Wdew = 0;
              all_vars_crop->veg_var[cridx][band].Wdmax = 0;
              all_vars_crop->snow[cridx][band].snow_canopy = 0;
            }
            else { //crop part
              // Convert LAI from global to local
              //veg_var[iveg][band].vegcover=1;
              all_vars_crop->veg_var[cridx][band].vegcover = veg_var[iveg][band].vegcover/veg_var[iveg][band].crop_frac;
              if (all_vars_crop->veg_var[cridx][band].vegcover > 1) all_vars_crop->veg_var[cridx][band].vegcover = 1;
              if (all_vars_crop->veg_var[cridx][band].vegcover < MIN_VEGCOVER) all_vars_crop->veg_var[cridx][band].vegcover = MIN_VEGCOVER;
              all_vars_crop->veg_var[cridx][band].albedo = (veg_var[iveg][band].albedo - (1-veg_var[iveg][band].crop_frac)*BARE_SOIL_ALBEDO)/veg_var[iveg][band].crop_frac;
              if (all_vars_crop->veg_var[cridx][band].albedo < H2O_SURF_ALBEDO) all_vars_crop->veg_var[cridx][band].albedo = H2O_SURF_ALBEDO;
              all_vars_crop->veg_var[cridx][band].LAI = veg_var[iveg][band].LAI/(all_vars_crop->veg_var[cridx][band].vegcover*veg_var[iveg][band].crop_frac);
              all_vars_crop->veg_var[cridx][band].Wdew = veg_var[iveg][band].Wdew/(all_vars_crop->veg_var[cridx][band].vegcover*veg_var[iveg][band].crop_frac);
              all_vars_crop->veg_var[cridx][band].Wdmax = all_vars_crop->veg_var[cridx][band].LAI*LAI_WATER_FACTOR;
              all_vars_crop->snow[cridx][band].snow_canopy = snow[iveg][band].snow_canopy/(all_vars_crop->veg_var[cridx][band].vegcover*veg_var[iveg][band].crop_frac);
            }
          }
        }
        else {
          // Convert LAI from global to local
          veg_var[iveg][band].LAI /= veg_var[iveg][band].vegcover;
          veg_var[iveg][band].Wdew /= veg_var[iveg][band].vegcover;
          veg_var[iveg][band].Wdmax = veg_var[iveg][band].LAI*LAI_WATER_FACTOR;
          snow[iveg][band].snow_canopy /= veg_var[iveg][band].vegcover;
        }
      }
    }
  }

  /***************************************
    Irrigation demand.
    Based on upper layer soil moisture.
  ******************************************/
  irrig = 0;
  if (options.IRRIGATION) {
    if (rec >= 0) {
      for(iveg = 0; iveg < Nveg; iveg++){
        veg_class = veg_con[iveg].veg_class; //ingjerd added
        for ( band = 0; band < Nbands; band++ ) {
          cell[iveg][band].irr_extract=0; //ingjerd added, for initializing. needed, or maybe not?
          veg_var[iveg][band].irrig=0.; //ingjerd added for initialization. hm... not sure about this!
          if (options.CROPFRAC  && veg_con[iveg].crop_frac_active)
          all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irrig=0.; //ingjerd added for initialization. hm... not sure about this!
          if(veg_lib[veg_class].irr_active[dmy[rec].month-1]) {  // only apply irrigation to specified crops
            if(snow[iveg][band].swq<0.001 && atmos->air_temp[NR]>7) {
              thresh_idx = 0;
              target_idx = 0;
              if (options.CROPFRAC && veg_con[iveg].crop_frac_active) {
                moistfract=all_vars_crop->cell[veg_con[iveg].crop_frac_idx+1][band].layer[target_idx].moist;
              }
              else if (veg_con[iveg].crop_frac_active)
              moistfract=cell[iveg][band].layer[target_idx].moist;
              if (veg_lib[veg_class].irr_sm_thresh == IRR_SAT)
              irr_sm_thresh=soil_con->max_moist[thresh_idx];
              else if (veg_lib[veg_class].irr_sm_thresh == IRR_FC)
              irr_sm_thresh=soil_con->Wcr[thresh_idx]/0.7;
              else
              irr_sm_thresh=soil_con->Wcr[thresh_idx]; // critical point in most cases
              if (irr_sm_thresh > soil_con->max_moist[thresh_idx])
              irr_sm_thresh = soil_con->max_moist[thresh_idx];
              if (veg_lib[veg_class].irr_sm_target == IRR_SAT)
              irr_sm_target=soil_con->max_moist[target_idx];
              else
              irr_sm_target=soil_con->Wcr[target_idx]/0.7; // field capacity in most cases
              if (irr_sm_target > soil_con->max_moist[target_idx])
              irr_sm_target = soil_con->max_moist[target_idx];
              if (options.CROPFRAC && veg_con[iveg].crop_frac_active) {
                if (moistfract < irr_sm_thresh)
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irr_apply = TRUE;
                else if (moistfract >= irr_sm_target)
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irr_apply = FALSE;
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx][band].irr_apply = FALSE;
              }
              else if (veg_con[iveg].crop_frac_active) {
                if (moistfract < irr_sm_thresh)
                veg_var[iveg][band].irr_apply = TRUE;
                else if (moistfract >= irr_sm_target)
                veg_var[iveg][band].irr_apply = FALSE;
              }

              if(options.CROPFRAC && veg_con[iveg].crop_frac_active) {
                irrig_est = (soil_con->max_moist[0]-all_vars_crop->cell[veg_con[iveg].crop_frac_idx+1][band].layer[target_idx].moist);
                if (all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irr_apply && atmos->prec[NR]<irrig_est)
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irrig = irrig_est-atmos->prec[NR];
                else
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irrig = 0;
                all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx][band].irrig = 0;

                Cv = veg_con[iveg].Cv;
                Cv *= veg_var[iveg][0].crop_frac;
                irrig += all_vars_crop->veg_var[veg_con[iveg].crop_frac_idx+1][band].irrig*Cv*soil_con->AreaFract[band];
              }
              else if (veg_con[iveg].crop_frac_active) { // veg_lib cropping calendar (irr can happen this month)
                irrig_est = (soil_con->max_moist[0]-cell[iveg][band].layer[target_idx].moist);
                if (veg_var[iveg][band].irr_apply && atmos->prec[NR]<irrig_est)
                veg_var[iveg][band].irrig = irrig_est-atmos->prec[NR]; //mm, irrigated area
                else
                veg_var[iveg][band].irrig = 0;
                Cv = veg_con[iveg].Cv;
                irrig += veg_var[iveg][band].irrig*Cv*soil_con->AreaFract[band]; //mm (cell average)
              }
            }
          }
        }
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
      if(options.CROPFRAC && veg_con[iveg].crop_frac_active) {
      for (cridx=veg_con[iveg].crop_frac_idx; cridx<veg_con[iveg].crop_frac_idx+2; cridx++) {

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

     if(cridx % 2 == 0) Cv *= (1-veg_var[iveg][0].crop_frac);
     else Cv *= veg_var[iveg][0].crop_frac;

      /**************************************************
        Initialize Model Parameters
      **************************************************/

      for(band = 0; band < Nbands; band++) {
	if(soil_con->AreaFract[band] > 0) {

	  /* Initialize energy balance variables */
	  all_vars_crop->energy[cridx][band].shortwave = 0;
	  all_vars_crop->energy[cridx][band].longwave  = 0.;

	  /* Initialize snow variables */
	  all_vars_crop->snow[cridx][band].vapor_flux        = 0.;
	  all_vars_crop->snow[cridx][band].canopy_vapor_flux = 0.;
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
          all_vars_crop->veg_var[cridx][band].rc = HUGE_RESIST;
        }
      }

      /** Assign wind_h **/
      /** Note: this is ignored below **/
      if (cridx % 2 == 0)
        wind_h = 0.1;
      else
        wind_h = veg_lib[veg_class].wind_h;

      /** Compute Surface Attenuation due to Vegetation Coverage **/
      surf_atten = (1-all_vars_crop->veg_var[cridx][0].vegcover)*1.0
                   + all_vars_crop->veg_var[cridx][0].vegcover
                   * exp(-veg_lib[veg_class].rad_atten * all_vars_crop->veg_var[cridx][0].LAI);

      /* Initialize soil thermal properties for the top two layers */
      prepare_full_energy(cridx, Nveg, options.Nnode, all_vars_crop, soil_con, moist0, ice0);

      /** Compute Bare (free of snow) Albedo **/
      bare_albedo = all_vars_crop->veg_var[cridx][0].albedo;

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
        if (cridx % 2 == 0) {
          displacement[0] = 0.1;
          roughness[0]    = soil_con->rough;
          overstory       = 0;
        }
        else {
          displacement[0] = veg_lib[pet_veg_class].displacement[dmy[rec].month-1];
          roughness[0]    = veg_lib[pet_veg_class].roughness[dmy[rec].month-1];
          overstory       = veg_lib[pet_veg_class].overstory;
	  if (p >= N_PET_TYPES_NON_NAT)
            if ( roughness[0] == 0 ) roughness[0] = soil_con->rough;
        }

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
          all_vars_crop->cell[cridx][band].aero_resist[0] = aero_resist[N_PET_TYPES][0];
          all_vars_crop->cell[cridx][band].aero_resist[1] = aero_resist[N_PET_TYPES][1];
        }
      }

      /******************************
        Compute nitrogen scaling factors and initialize other veg vars
      ******************************/
      if (options.CARBON && iveg < Nveg) {
	for(band=0; band<Nbands; band++) {
          for (cidx=0; cidx<options.Ncanopy; cidx++) {
            all_vars_crop->veg_var[cridx][band].rsLayer[cidx] = HUGE_RESIST;
          }
          all_vars_crop->veg_var[cridx][band].aPAR = 0;
          if (dmy->hour == 0) {
            calc_Nscale_factors(veg_lib[veg_class].NscaleFlag,
                                veg_con[iveg].CanopLayerBnd,
                                veg_lib[veg_class].LAI[dmy[rec].month-1],
                                soil_con->lat,
                                soil_con->lng,
                                soil_con->time_zone_lng,
                                dmy[rec],
                                all_vars_crop->veg_var[cridx][band].NscaleFactor);
          }
          if (dmy[rec].month == 1 && dmy[rec].day == 1) {
            all_vars_crop->veg_var[cridx][band].AnnualNPPPrev = all_vars_crop->veg_var[cridx][band].AnnualNPP;
            all_vars_crop->veg_var[cridx][band].AnnualNPP = 0;
          }
        }
      }

      /******************************
        Solve ground surface fluxes (for crops = crop_active
      ******************************/

      for ( band = 0; band < Nbands; band++ ) {
	if( soil_con->AreaFract[band] > 0 ) {

	  lag_one     = veg_con[iveg].lag_one;
	  sigma_slope = veg_con[iveg].sigma_slope;
	  fetch       = veg_con[iveg].fetch;

	  /* Initialize pot_evap */
	  for (p=0; p<N_PET_TYPES; p++)
	    all_vars_crop->cell[cridx][band].pot_evap[p] = 0;

          // Limit irrigation to available water
          if(options.IRRIGATION && veg_lib[veg_class].irr_active[dmy[rec].month-1] && (!options.CROPFRAC || cridx % 2 == 1) && !options.IRR_FREE && atmos->air_temp[NR]>7) {

	    all_vars_crop->veg_var[cridx][band].irrig = 0.; //initialize

	    // Extract irrigation from local water as necessary. Local water takes priority over irr_run and irr_with
	    extract_water=0.; //cell average
	    localextract=1;
	    if (localextract==1 && irrig>0) {
	      // extract water from runoff, then baseflow, from each already simulated vegtiles, until reaching irrig.
	      // forutsetning: irrigated crop is last on vegparam list, and only one veg is irrigated.
	      for(ivegdummy = 0; ivegdummy < iveg; ivegdummy++) { //loop through previous ivegs
		cell[ivegdummy][band].irr_extract=0;
		if (veg_con[ivegdummy].Cv > 0.0) {
		  Cvdummy = veg_con[ivegdummy].Cv;
		  if (extract_water+cell[ivegdummy][band].runoff*Cvdummy*soil_con->AreaFract[band] <= irrig) {
		    extract_water += cell[ivegdummy][band].runoff*Cvdummy*soil_con->AreaFract[band];
		    cell[ivegdummy][band].irr_extract += cell[ivegdummy][band].runoff;
		    cell[ivegdummy][band].runoff = 0;
		  }
		  else if(extract_water<irrig) {
		    cell[ivegdummy][band].irr_extract += (irrig-extract_water)/(Cvdummy*soil_con->AreaFract[band]);
		    cell[ivegdummy][band].runoff -= (irrig-extract_water)/(Cvdummy*soil_con->AreaFract[band]);
		    extract_water += irrig-extract_water;
		  }
		  if (extract_water+cell[ivegdummy][band].baseflow*Cvdummy*soil_con->AreaFract[band] <= irrig) {
		    extract_water += cell[ivegdummy][band].baseflow*Cvdummy*soil_con->AreaFract[band];
		    cell[ivegdummy][band].irr_extract += cell[ivegdummy][band].baseflow;
		    cell[ivegdummy][band].baseflow = 0;
		  }
		  else if(extract_water<irrig) {
		    if(rec==recdummy)
		      fprintf(stderr,"ivegdummyD %d Cvdummy %f missing %.3f base init average %f\n",
		            ivegdummy,Cvdummy,irrig-extract_water,cell[ivegdummy][band].baseflow*Cvdummy*soil_con->AreaFract[band]);
		    cell[ivegdummy][band].irr_extract += (irrig-extract_water)/(Cvdummy*soil_con->AreaFract[band]);
		    cell[ivegdummy][band].baseflow -= (irrig-extract_water)/(Cvdummy*soil_con->AreaFract[band]);
		    extract_water += (irrig-extract_water);
		     if(rec==recdummy)
		       fprintf(stderr,"\t missing %.3f base final average %f\n\n",
			       irrig-extract_water,cell[ivegdummy][band].baseflow*Cvdummy*soil_con->AreaFract[band]);
		  }
		} //if current ivegdummy>0
	      } //loop through previous tiles
	    } //if localextract=1 and irrig>0

	    if(irrig>extract_water) { //still need more
	      extract_water+=atmos->irr_run[NR]+atmos->irr_with[NR];
	      if(extract_water>irrig)
		extract_water=irrig;
	    }

	    all_vars_crop->veg_var[cridx][band].irrig = extract_water/(Cv*soil_con->AreaFract[band]); //local cridx tile, irrig possible

	  } //options.irrigation limited

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
				     atmos, dmy, &(all_vars_crop->energy[cridx][band]), gp,
				     &(all_vars_crop->cell[cridx][band]),
				     &(all_vars_crop->snow[cridx][band]),
				     soil_con, &(all_vars_crop->veg_var[cridx][band]),
				     lag_one, sigma_slope, fetch, veg_con[iveg].CanopLayerBnd);

	  if ( ErrorFlag == ERROR ) return ( ERROR );

	  atmos->out_prec += out_prec[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_rain += out_rain[band*2] * Cv * soil_con->AreaFract[band];
	  atmos->out_snow += out_snow[band*2] * Cv * soil_con->AreaFract[band];

          /********************************************************
            Compute soil wetness and root zone soil moisture
          ********************************************************/
          all_vars_crop->cell[cridx][band].rootmoist = 0;
          all_vars_crop->cell[cridx][band].wetness = 0;
          for(lidx=0;lidx<options.Nlayer;lidx++) {
            if (veg_con->root[lidx] > 0) {
              all_vars_crop->cell[cridx][band].rootmoist += all_vars_crop->cell[cridx][band].layer[lidx].moist;
            }
	    all_vars_crop->cell[cridx][band].wetness += (all_vars_crop->cell[cridx][band].layer[lidx].moist - soil_con->Wpwp[lidx])/(soil_con->porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
          }
          all_vars_crop->cell[cridx][band].wetness /= options.Nlayer;

	} /** End non-zero area band **/

      } /** End Loop Through Elevation Bands **/

      } /** end loop over crop subtiles **/
      // Transfer states and fluxes from crop subtiles back to original data structures

      for ( band = 0; band < Nbands; band++ ) {
	if( soil_con->AreaFract[band] > 0 ) {

	  //ingjerd added the below sentence, for wb tracking with irrigation and crop_frac included . more is needed?*/
	  veg_var[iveg][band].irrig=all_vars_crop->veg_var[1][band].irrig*veg_var[iveg][band].crop_frac;

	  // Copy veg_var state data
	  veg_var[iveg][band].Wdew = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Wdew + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Wdew;
	  if (options.CARBON) {
	    for (cidx=0; cidx<options.Ncanopy; cidx++) {
	      veg_var[iveg][band].NscaleFactor[cidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].NscaleFactor[cidx] + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].NscaleFactor[cidx];
	      veg_var[iveg][band].aPARLayer[cidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].aPARLayer[cidx] + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].aPARLayer[cidx];
	      veg_var[iveg][band].CiLayer[cidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].CiLayer[cidx] + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].CiLayer[cidx];
	      veg_var[iveg][band].rsLayer[cidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].rsLayer[cidx] + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].rsLayer[cidx];
	    }
	  }
	  veg_var[iveg][band].Ci = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Ci + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Ci;
	  tmp_double = (1-veg_var[iveg][band].crop_frac)/all_vars_crop->veg_var[0][band].rc + veg_var[iveg][band].crop_frac/all_vars_crop->veg_var[1][band].rc;
	  veg_var[iveg][band].rc = 1/tmp_double;
	  veg_var[iveg][band].NPPfactor = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].NPPfactor + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].NPPfactor;
	  veg_var[iveg][band].AnnualNPP = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].AnnualNPP + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].AnnualNPP;
	  veg_var[iveg][band].AnnualNPPPrev = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].AnnualNPPPrev + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].AnnualNPPPrev;

	  // Copy veg_var flux data
	  veg_var[iveg][band].canopyevap = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].canopyevap + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].canopyevap;
	  veg_var[iveg][band].throughfall = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].throughfall + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].throughfall;
	  veg_var[iveg][band].aPAR = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].aPAR + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].aPAR;
	  veg_var[iveg][band].GPP = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].GPP + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].GPP;
	  veg_var[iveg][band].Rphoto = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Rphoto + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Rphoto;
	  veg_var[iveg][band].Rdark = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Rdark + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Rdark;
	  veg_var[iveg][band].Rmaint = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Rmaint + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Rmaint;
	  veg_var[iveg][band].Rgrowth = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Rgrowth + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Rgrowth;
	  veg_var[iveg][band].Raut = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Raut + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Raut;
	  veg_var[iveg][band].NPP = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].NPP + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].NPP;
	  veg_var[iveg][band].Litterfall = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->veg_var[0][band].Litterfall + veg_var[iveg][band].crop_frac*all_vars_crop->veg_var[1][band].Litterfall;

	  // Copy cell state data
	  for (cidx=0; cidx<2; cidx++) {
	    tmp_double = (1-veg_var[iveg][band].crop_frac)/all_vars_crop->cell[0][band].aero_resist[cidx] + veg_var[iveg][band].crop_frac/all_vars_crop->cell[1][band].aero_resist[cidx];
	    cell[iveg][band].aero_resist[cidx] = 1/tmp_double;
	  }

	  cell[iveg][band].asat = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].asat + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].asat;
	  cell[iveg][band].CLitter = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].CLitter + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].CLitter;
	  cell[iveg][band].CInter = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].CInter + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].CInter;
	  cell[iveg][band].CSlow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].CSlow + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].CSlow;

	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    cell[iveg][band].layer[lidx].bare_evap_frac = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].bare_evap_frac + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].bare_evap_frac;

	    cell[iveg][band].layer[lidx].Cs = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].Cs + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].Cs;

	    cell[iveg][band].layer[lidx].kappa = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].kappa + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].kappa;

	    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
	      cell[iveg][band].layer[lidx].ice[frost_area] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].ice[frost_area] + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].ice[frost_area];

	    cell[iveg][band].layer[lidx].moist = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].moist + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].moist;

	    cell[iveg][band].layer[lidx].phi = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].phi + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].phi;

	    cell[iveg][band].layer[lidx].zwt = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].zwt + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].zwt;

	    cell[iveg][band].layer[lidx].T = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].T + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].T;
	  }

	  // Copy cell flux data
	  cell[iveg][band].baseflow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].baseflow + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].baseflow;
	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    cell[iveg][band].layer[lidx].evap = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].layer[lidx].evap + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].layer[lidx].evap;
	  }
	  cell[iveg][band].inflow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].inflow + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].inflow;
	  for ( cidx = 0; cidx < N_PET_TYPES; cidx++ ) {
	    cell[iveg][band].pot_evap[cidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].pot_evap[cidx] + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].pot_evap[cidx];
	  }
	  cell[iveg][band].runoff = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].runoff + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].runoff;
	  cell[iveg][band].RhLitter = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].RhLitter + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].RhLitter;
	  cell[iveg][band].RhLitter2Atm = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].RhLitter2Atm + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].RhLitter2Atm;
	  cell[iveg][band].RhInter = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].RhInter + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].RhInter;
	  cell[iveg][band].RhSlow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].RhSlow + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].RhSlow;
	  cell[iveg][band].RhTot = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].RhTot + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].RhTot;
	  cell[iveg][band].rootmoist = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].rootmoist + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].rootmoist;
	  cell[iveg][band].wetness = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].wetness + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].wetness;
	  cell[iveg][band].zwt = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].zwt + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].zwt;
	  cell[iveg][band].zwt_lumped = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->cell[0][band].zwt_lumped + veg_var[iveg][band].crop_frac*all_vars_crop->cell[1][band].zwt_lumped;

	  // Copy snow state data
	  snow[iveg][band].albedo = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].albedo + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].albedo;
	  snow[iveg][band].canopy_albedo = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].canopy_albedo + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].canopy_albedo;
	  snow[iveg][band].coldcontent = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].coldcontent + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].coldcontent;
	  snow[iveg][band].coverage = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].coverage + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].coverage;
	  snow[iveg][band].density = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].density + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].density;
	  snow[iveg][band].depth = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].depth + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].depth;
	  snow[iveg][band].last_snow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].last_snow + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].last_snow;
	  snow[iveg][band].max_snow_depth = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].max_snow_depth + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].max_snow_depth;
	  snow[iveg][band].MELTING = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].MELTING + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].MELTING;
	  snow[iveg][band].pack_temp = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].pack_temp + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].pack_temp;
	  snow[iveg][band].pack_water = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].pack_water + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].pack_water;
	  snow[iveg][band].snow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].snow + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].snow;
	  snow[iveg][band].snow_canopy = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].snow_canopy + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].snow_canopy;
	  snow[iveg][band].store_coverage = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].store_coverage + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].store_coverage;
	  snow[iveg][band].store_snow = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].store_snow + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].store_snow;
	  snow[iveg][band].store_swq = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].store_swq + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].store_swq;
	  snow[iveg][band].surf_temp = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].surf_temp + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].surf_temp;
	  snow[iveg][band].surf_temp_fbcount = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].surf_temp_fbcount + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].surf_temp_fbcount;
	  snow[iveg][band].surf_temp_fbflag = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].surf_temp_fbflag + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].surf_temp_fbflag;
	  snow[iveg][band].surf_water = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].surf_water + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].surf_water;
	  snow[iveg][band].swq = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].swq + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].swq;
	  snow[iveg][band].snow_distrib_slope = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].snow_distrib_slope + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].snow_distrib_slope;
	  snow[iveg][band].tmp_int_storage = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].tmp_int_storage + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].tmp_int_storage;

	  // Copy snow flux data
	  snow[iveg][band].blowing_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].blowing_flux + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].blowing_flux;
	  snow[iveg][band].canopy_vapor_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].canopy_vapor_flux + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].canopy_vapor_flux;
	  snow[iveg][band].mass_error = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].mass_error + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].mass_error;
	  snow[iveg][band].melt = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].melt + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].melt;
	  snow[iveg][band].Qnet = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].Qnet + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].Qnet;
	  snow[iveg][band].surface_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].surface_flux + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].surface_flux;
	  snow[iveg][band].transport = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].transport + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].transport;
	  snow[iveg][band].vapor_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->snow[0][band].vapor_flux + veg_var[iveg][band].crop_frac*all_vars_crop->snow[1][band].vapor_flux;

	  // Copy energy state data
	  energy[iveg][band].AlbedoLake = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AlbedoLake + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AlbedoLake;
	  energy[iveg][band].AlbedoOver = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AlbedoOver + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AlbedoOver;
	  energy[iveg][band].AlbedoUnder = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AlbedoUnder + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AlbedoUnder;
	  energy[iveg][band].frozen = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].frozen + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].frozen;
	  energy[iveg][band].Nfrost = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Nfrost + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Nfrost;
	  energy[iveg][band].Nthaw = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Nthaw + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Nthaw;
	  energy[iveg][band].T1_index = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].T1_index + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].T1_index;
	  energy[iveg][band].Tcanopy = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tcanopy + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tcanopy;
	  energy[iveg][band].Tcanopy_fbcount = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tcanopy_fbcount + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tcanopy_fbcount;
	  energy[iveg][band].Tcanopy_fbflag = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tcanopy_fbflag + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tcanopy_fbflag;
	  energy[iveg][band].Tfoliage = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tfoliage + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tfoliage;
	  energy[iveg][band].Tfoliage_fbcount = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tfoliage_fbcount + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tfoliage_fbcount;
	  energy[iveg][band].Tfoliage_fbflag = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tfoliage_fbflag + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tfoliage_fbflag;
	  energy[iveg][band].Tsurf = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tsurf + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tsurf;
	  energy[iveg][band].Tsurf_fbcount = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tsurf_fbcount + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tsurf_fbcount;
	  energy[iveg][band].Tsurf_fbflag = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Tsurf_fbflag + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Tsurf_fbflag;
	  energy[iveg][band].unfrozen = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].unfrozen + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].unfrozen;
	  for (lidx=0; lidx<2; lidx++) {
	    energy[iveg][band].Cs[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Cs[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Cs[lidx];
	    energy[iveg][band].kappa[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].kappa[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].kappa[lidx];
	  }
	  for (index=0; index<options.Nnode; index++) {
	    energy[iveg][band].Cs[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].Cs[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].Cs[lidx];
	    energy[iveg][band].fdepth[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].fdepth[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].fdepth[lidx];
	    energy[iveg][band].ice[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].ice[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].ice[lidx];
	    energy[iveg][band].kappa_node[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].kappa_node[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].kappa_node[lidx];
	    energy[iveg][band].moist[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].moist[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].moist[lidx];
	    energy[iveg][band].T[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].T[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].T[lidx];
	    energy[iveg][band].T_fbcount[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].T_fbcount[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].T_fbcount[lidx];
	    energy[iveg][band].T_fbflag[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].T_fbflag[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].T_fbflag[lidx];
	    energy[iveg][band].tdepth[lidx] = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].tdepth[lidx] + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].tdepth[lidx];
	  }

	  // Copy energy flux data
	  energy[iveg][band].advected_sensible = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].advected_sensible + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].advected_sensible;
	  energy[iveg][band].advection = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].advection + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].advection;
	  energy[iveg][band].AtmosError = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AtmosError + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AtmosError;
	  energy[iveg][band].AtmosLatent = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AtmosLatent + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AtmosLatent;
	  energy[iveg][band].AtmosLatentSub = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AtmosLatentSub + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AtmosLatentSub;
	  energy[iveg][band].AtmosSensible = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].AtmosSensible + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].AtmosSensible;
	  energy[iveg][band].canopy_advection = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].canopy_advection + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].canopy_advection;
	  energy[iveg][band].canopy_latent = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].canopy_latent + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].canopy_latent;
	  energy[iveg][band].canopy_latent_sub = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].canopy_latent_sub + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].canopy_latent_sub;
	  energy[iveg][band].canopy_refreeze = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].canopy_refreeze + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].canopy_refreeze;
	  energy[iveg][band].canopy_sensible = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].canopy_sensible + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].canopy_sensible;
	  energy[iveg][band].deltaCC = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].deltaCC + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].deltaCC;
	  energy[iveg][band].deltaH = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].deltaH + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].deltaH;
	  energy[iveg][band].error = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].error + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].error;
	  energy[iveg][band].fusion = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].fusion + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].fusion;
	  energy[iveg][band].grnd_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].grnd_flux + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].grnd_flux;
	  energy[iveg][band].latent = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].latent + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].latent;
	  energy[iveg][band].latent_sub = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].latent_sub + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].latent_sub;
	  energy[iveg][band].longwave = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].longwave + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].longwave;
	  energy[iveg][band].LongOverIn = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].LongOverIn + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].LongOverIn;
	  energy[iveg][band].LongUnderIn = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].LongUnderIn + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].LongUnderIn;
	  energy[iveg][band].LongUnderOut = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].LongUnderOut + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].LongUnderOut;
	  energy[iveg][band].melt_energy = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].melt_energy + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].melt_energy;
	  energy[iveg][band].NetLongAtmos = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetLongAtmos + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetLongAtmos;
	  energy[iveg][band].NetLongOver = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetLongOver + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetLongOver;
	  energy[iveg][band].NetLongUnder = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetLongUnder + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetLongUnder;
	  energy[iveg][band].NetShortAtmos = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetShortAtmos + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetShortAtmos;
	  energy[iveg][band].NetShortGrnd = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetShortGrnd + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetShortGrnd;
	  energy[iveg][band].NetShortOver = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetShortOver + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetShortOver;
	  energy[iveg][band].NetShortUnder = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].NetShortUnder + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].NetShortUnder;
	  energy[iveg][band].out_long_canopy = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].out_long_canopy + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].out_long_canopy;
	  energy[iveg][band].out_long_surface = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].out_long_surface + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].out_long_surface;
	  energy[iveg][band].refreeze_energy = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].refreeze_energy + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].refreeze_energy;
	  energy[iveg][band].sensible = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].sensible + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].sensible;
	  energy[iveg][band].shortwave = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].shortwave + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].shortwave;
	  energy[iveg][band].ShortOverIn = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].ShortOverIn + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].ShortOverIn;
	  energy[iveg][band].ShortUnderIn = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].ShortUnderIn + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].ShortUnderIn;
	  energy[iveg][band].snow_flux = (1-veg_var[iveg][band].crop_frac)*all_vars_crop->energy[0][band].snow_flux + veg_var[iveg][band].crop_frac*all_vars_crop->energy[1][band].snow_flux;


        }
      }

      } //end irractive for veg type?

      else {

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

          // Limit irrigation to available water
          if(options.IRRIGATION && veg_lib[veg_class].irr_active[dmy[rec].month-1] && !options.IRR_FREE) {
 	    if(irrig>(atmos->irr_run[NR]+atmos->irr_with[NR]) && atmos->air_temp[NR]>7)
	    veg_var[iveg][band].irrig *= (atmos->irr_run[NR]+atmos->irr_with[NR])/irrig;
          }

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

     } /** end if crop_frac_active **/

    } /** end non-zero area veg tile **/
  } /** end of vegetation loop **/
  /* Convert LAI back to global */
  if (rec >= 0) {
    for(iveg = 0; iveg < Nveg; iveg++){
      for ( band = 0; band < Nbands; band++ ) {
        // Handle non-crop tiles
        if (!options.CROPFRAC || !veg_con[iveg].crop_frac_active) {
          veg_var[iveg][band].LAI *= veg_var[iveg][band].vegcover;
          veg_var[iveg][band].Wdmax *= veg_var[iveg][band].vegcover;
        }
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
