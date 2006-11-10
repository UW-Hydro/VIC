#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
              soil_con_struct   *soil_con,
	      veg_con_struct    *veg_con,
              lake_con_struct   *lake_con,
              out_data_file_struct   *out_data_files,
              out_data_struct   *out_data,
              save_data_struct  *save_data,
	      dmy_struct        *dmy,
              int                rec,
	      int                Nnodes)
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
  6-8-2000 modified to handle spatially distribute frozen soil     KAC
  10-6-2000 modified to handle partial snow cover                  KAC
  02-27-01 modified to output lake model variables                 KAC
  11-18-02 updated output of lake variables to reflect algorithm
           changes.  Also added output variables for blowing snow
           algorithm.                                              LCB
  03-12-03 modified to add additional energy balance variable storage
           when output of snow bands is selected.                  KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  30-Oct-03 Snow_flux was incorrectly set to Tcanopy.  Fixed.   TJB
  25-Aug-04 Sub_snow was incorrectly set to blowing_flux.  Now it is
            set to vapor_flux.                                  TJB
  28-Sep-04 Now out_data->aero_resist stores the aerodynamic resistance
            used in flux calculations.                          TJB
  2005-Mar-24 Modified to compute ALMA output variables.        TJB
  2005-Apr-23 Now aero_cond is aggregated instead of aero_resist.       TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT options; uses the
	      new save_data structure; implemented aggregation.  TJB
  2006-Oct-10 Shortened the names of output variables whose names were
	      too long; fixed typos in others; created new OUT_IN_LONG
	      variable.  TJB
  2006-Nov-07 Added OUT_SOIL_TNODE.  TJB
  2006-Nov-07 Assigned value to overstory.  TJB
  2006-Nov-07 Removed LAKE_MODEL option. TJB

**********************************************************************/
{
  extern global_param_struct global_param;
  extern veg_lib_struct  *veg_lib;
  extern option_struct    options;
#if LINK_DEBUG
  extern debug_struct     debug;
#endif

  int                     veg;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  int                     overstory;
#if SPATIAL_FROST
  int                     frost_area;
#endif
  char              *AboveTreeLine;
  double            *AreaFract;
  double            *depth;
  double            *dz;
#if SPATIAL_FROST
  double            *frost_fract;
  double             frost_slope;
#endif // SPATIAL_FROST
  double             dp;
  int                skipyear;
  double                  Cv;
  double                  Clake;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
  double                  tmp_ice;
  double                  tmp_fract;
  double                  cv_baresoil;
  double                  cv_veg;
  double                  cv_snow;
  double                  rad_temp;
  double                  surf_temp;
  double                  inflow;
  double                  outflow;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];
  int                     n;
  int                     v;
  int                     i;
  int                     dt_sec;
  int                     out_dt_sec;
  int                     out_step_ratio;
  static int              step_count;

  cell_data_struct     ***cell;
  energy_bal_struct     **energy;
  lake_var_struct         lake_var;
  snow_data_struct      **snow;
  veg_var_struct       ***veg_var;

  AboveTreeLine = soil_con->AboveTreeLine;
  AreaFract = soil_con->AreaFract;
  depth = soil_con->depth;
  dz = soil_con->dz_node;
#if SPATIAL_FROST
  frost_fract = soil_con->frost_fract;
  frost_slope = soil_con->frost_slope;
#endif // SPATIAL_FROST
  dp = soil_con->dp;
  skipyear = global_param.skipyear;
  dt_sec = global_param.dt*SECPHOUR;
  out_dt_sec = global_param.out_dt*SECPHOUR;
  out_step_ratio = (int)(out_dt_sec/dt_sec);
  step_count++;

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

  cv_baresoil = 0;
  cv_veg = 0;
  cv_snow = 0;

  // Initialize output data to zero
  zero_output_list(out_data);

  // Set output versions of input forcings
  out_data[OUT_AIR_TEMP].data[0]  = atmos->air_temp[NR];
  out_data[OUT_DENSITY].data[0]   = atmos->density[NR];
  out_data[OUT_LONGWAVE].data[0]  = atmos->longwave[NR];
  out_data[OUT_PREC].data[0]      = atmos->out_prec;
  out_data[OUT_PRESSURE].data[0]  = atmos->pressure[NR];
  out_data[OUT_QAIR].data[0]      = EPS * atmos->vp[NR]/atmos->pressure[NR];
  out_data[OUT_RAINF].data[0]     = atmos->out_rain;
  out_data[OUT_REL_HUMID].data[0] = 100.*atmos->vp[NR]/(atmos->vp[NR]+atmos->vpd[NR]);
  out_data[OUT_SHORTWAVE].data[0] = atmos->shortwave[NR];
  out_data[OUT_SNOWF].data[0]     = atmos->out_snow;
  out_data[OUT_VP].data[0]        = atmos->vp[NR];
  out_data[OUT_WIND].data[0]      = atmos->wind[NR];
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp->cell;
  energy  = prcp->energy;
  lake_var = prcp->lake_var;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;
 
  /****************************************
    Store Output for all Vegetation Types (except lakes)
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {

    if ( veg < veg_con[0].vegetat_type_num ) 
      Cv = veg_con[veg].Cv;
    else
      Cv = (1.0 - veg_con[0].Cv_sum);

    if ( Cv > 0 ) {

      overstory = veg_lib[veg_con[veg].veg_class].overstory;

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
	  if(AreaFract[band] > 0. && ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !overstory)))) {

	    tmp_evap = 0.0;
	    for(index=0;index<options.Nlayer;index++)
	      tmp_evap += cell[dist][veg][band].layer[index].evap;
	    if ( veg < veg_con[0].vegetat_type_num )
	      out_data[OUT_TRANSP_VEG].data[0] += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    else 
	      out_data[OUT_EVAP_BARE].data[0] += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];

	    tmp_evap += snow[veg][band].vapor_flux * 1000.;
	    out_data[OUT_SUB_SNOW].data[0] += snow[veg][band].vapor_flux * 1000. 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    out_data[OUT_SUB_SURFACE].data[0] += snow[veg][band].surface_flux * 1000. 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    out_data[OUT_SUB_BLOWING].data[0] += snow[veg][band].blowing_flux * 1000. 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    if ( veg <= veg_con[0].vegetat_type_num ) {
	      tmp_evap += snow[veg][band].canopy_vapor_flux * 1000.;
	      out_data[OUT_SUB_CANOP].data[0] += snow[veg][band].canopy_vapor_flux 
		* 1000. * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    }
	    if ( veg < veg_con[0].vegetat_type_num ) {
	      tmp_evap += veg_var[dist][veg][band].canopyevap;
	      out_data[OUT_EVAP_CANOP].data[0] += veg_var[dist][veg][band].canopyevap 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    }
	    out_data[OUT_EVAP].data[0] += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 

	    /** record runoff **/
	    out_data[OUT_RUNOFF].data[0]   += cell[dist][veg][band].runoff 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** record baseflow **/
	    out_data[OUT_BASEFLOW].data[0] += cell[dist][veg][band].baseflow 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 

	    /** record inflow **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data[OUT_INFLOW].data[0] += (cell[dist][veg][band].inflow 
				   + veg_var[dist][veg][band].canopyevap) 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    else 
	      out_data[OUT_INFLOW].data[0] += (cell[dist][veg][band].inflow) 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** record canopy interception **/
	    if ( veg < veg_con[0].vegetat_type_num ) 
	      out_data[OUT_WDEW].data[0] += veg_var[dist][veg][band].Wdew 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	  
            /** record aerodynamic conductance **/
            if (cell[WET][veg][0].aero_resist_used > SMALL) {
              out_data[OUT_AERO_COND].data[0] += (1/cell[WET][veg][0].aero_resist_used)
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }
            else {
              out_data[OUT_AERO_COND].data[0] = HUGE_RESIST;
              out_data[OUT_AERO_RESIST].data[0] = cell[WET][veg][0].aero_resist_used;
            }
	    
	    /** record layer moistures **/
	    for(index=0;index<options.Nlayer;index++) {
	      tmp_moist = cell[dist][veg][band].layer[index].moist;
#if SPATIAL_FROST
              tmp_ice = 0;
              for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
                tmp_ice  += (cell[dist][veg][band].layer[index].ice[frost_area]
                             * frost_fract[frost_area]);
#else
              tmp_ice   = cell[dist][veg][band].layer[index].ice;
#endif
	      tmp_moist -= tmp_ice;
	      if(options.MOISTFRACT) {
		tmp_moist /= depth[index] * 1000.;
		tmp_ice /= depth[index] * 1000.;
	      }
	      out_data[OUT_SOIL_LIQ].data[index] += tmp_moist
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      out_data[OUT_SOIL_ICE].data[index] += tmp_ice
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }
	    out_data[OUT_SOIL_WET].data[0] += cell[dist][veg][band].wetness
              * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    out_data[OUT_ROOTMOIST].data[0] += cell[dist][veg][band].rootmoist
              * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];

	    /** record layer temperatures **/
	    for(index=0;index<options.Nlayer;index++) {
	      out_data[OUT_SOIL_TEMP].data[index] += cell[dist][veg][band].layer[index].T
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }

	    /** record thermal node temperatures **/
	    for(index=0;index<options.Nnode;index++) {
	      out_data[OUT_SOIL_TNODE].data[index] += energy[veg][band].T[index]
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }

	  }
	}
      }

      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0. && ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !overstory)))) {

	  /**********************************
	    Record Frozen Soil Variables
	  **********************************/

	  /** record freezing and thawing front depths **/
	  if(options.FROZEN_SOIL) {
	    for(index = 0; index < MAX_FRONTS; index++) {
	      if(energy[veg][band].fdepth[index] != MISSING)
		out_data[OUT_FDEPTH].data[index] += energy[veg][band].fdepth[index] 
		  * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	      if(energy[veg][band].tdepth[index] != MISSING)
		out_data[OUT_TDEPTH].data[index] += energy[veg][band].tdepth[index] 
		  * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	    }
	  }

#if SPATIAL_FROST
          tmp_fract = 0;
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
            if ( cell[0][veg][band].layer[0].ice[frost_area] )
              tmp_fract  += frost_fract[frost_area];
#else
          if ( cell[0][veg][band].layer[0].ice > 0 )
            tmp_fract   = 1.;
#endif
          out_data[OUT_SURF_FROST_FRAC].data[0] += tmp_fract * Cv * AreaFract[band] * TreeAdjustFactor[band];

          tmp_fract = 0;
#if SPATIAL_FROST
          if ( (energy[veg][band].T[0] + frost_slope / 2.) > 0 ) {
            if ( (energy[veg][band].T[0] - frost_slope / 2.) <= 0 )
              tmp_fract
                += linear_interp( 0, (energy[veg][band].T[0]
                                      + frost_slope / 2.),
                                  (energy[veg][band].T[0]
                                   - frost_slope / 2.), 1, 0)
                * Cv * AreaFract[band] * TreeAdjustFactor[band];
          }
          else
            tmp_fract += 1 * Cv * AreaFract[band] * TreeAdjustFactor[band];
#else
          if ( energy[veg][band].T[0] <= 0 )
            tmp_fract = 1 * Cv * AreaFract[band] * TreeAdjustFactor[band];
#endif

	  /**********************************
            Record Energy Balance Variables
	  **********************************/

	  /** record surface radiative temperature **/
          if ( overstory && snow[veg][band].snow )
            rad_temp = energy[veg][band].Tcanopy + KELVIN;
          else
            rad_temp = energy[veg][band].Tsurf + KELVIN;

          /** record landcover temperature **/
          if(veg == veg_con[0].vegetat_type_num) {
            // landcover is bare soil
            out_data[OUT_BARESOILT].data[0] += (rad_temp-KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
            cv_baresoil += AreaFract[band] * TreeAdjustFactor[band];
          }
          else {
            // landcover is vegetation
            if ( overstory && !snow[veg][band].snow )
              // here, rad_temp will be wrong since it will pick the understory temperature
              out_data[OUT_VEGT].data[0] += energy[veg][band].Tfoliage * Cv * AreaFract[band] * TreeAdjustFactor[band];
            else
              out_data[OUT_VEGT].data[0] += (rad_temp-KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
            cv_veg += Cv * AreaFract[band] * TreeAdjustFactor[band];
          }

	  /** record mean surface temperature [C]  **/
	  out_data[OUT_SURF_TEMP].data[0] += (energy[veg][band].T[0] + energy[veg][band].T[1])/2. * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record net shortwave radiation **/
	  out_data[OUT_NET_SHORT].data[0] += energy[veg][band].NetShortAtmos
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record net longwave radiation **/
	  out_data[OUT_NET_LONG].data[0]  += energy[veg][band].NetLongAtmos
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
          /** record incoming longwave radiation at ground surface (under veg) **/
          if ( snow[veg][band].snow && overstory )
            out_data[OUT_IN_LONG].data[0] += energy[veg][band].LongOverIn
              * Cv * AreaFract[band] * TreeAdjustFactor[band];
          else
            out_data[OUT_IN_LONG].data[0] += energy[veg][band].LongUnderIn
              * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record albedo **/
          if ( snow[veg][band].snow && overstory )
            out_data[OUT_ALBEDO].data[0]    += energy[veg][band].AlbedoOver
              * Cv * AreaFract[band] * TreeAdjustFactor[band];
          else
            out_data[OUT_ALBEDO].data[0]    += energy[veg][band].AlbedoUnder
              * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record latent heat flux **/
	  out_data[OUT_LATENT].data[0]    -= energy[veg][band].AtmosLatent
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

          /** record latent heat flux from sublimation **/
          out_data[OUT_LATENT_SUB].data[0]
            -= energy[veg][band].AtmosLatentSub * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record sensible heat flux **/
	  out_data[OUT_SENSIBLE].data[0]  -= energy[veg][band].AtmosSensible
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record ground heat flux (+ heat storage) **/
	  out_data[OUT_GRND_FLUX].data[0] -= energy[veg][band].grnd_flux
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record heat storage **/
	  out_data[OUT_DELTAH].data[0]    -= energy[veg][band].deltaH
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record heat of fusion **/
	  out_data[OUT_FUSION].data[0]    -= energy[veg][band].fusion
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record energy balance error **/
	  out_data[OUT_ENERGY_ERROR].data[0] += energy[veg][band].error
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record radiative effective temperature [K], 
	      emissivities set = 1.0  **/
	  out_data[OUT_RAD_TEMP].data[0] += ((rad_temp) * (rad_temp) 
				 * (rad_temp) * (rad_temp)) 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /*****************************
	    Record Snow Pack Variables 
	  *****************************/
	  
	  /** record snow water equivalence **/
	  out_data[OUT_SWE].data[0]
	    += snow[veg][band].swq * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack depth **/
	  out_data[OUT_SNOW_DEPTH].data[0]
	    += snow[veg][band].depth * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
	  
          /** record snowpack albedo, temperature **/
          if (snow[veg][band].swq> 0.0) {
            out_data[OUT_SALBEDO].data[0]
              += snow[veg][band].albedo * Cv * AreaFract[band] * TreeAdjustFactor[band];
            out_data[OUT_SNOW_SURF_TEMP].data[0]
              += snow[veg][band].surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
            out_data[OUT_SNOW_PACK_TEMP].data[0]
              += snow[veg][band].pack_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
            cv_snow += Cv * AreaFract[band] * TreeAdjustFactor[band];
          }

	  /** record canopy intercepted snow **/
	  if ( veg < veg_con[0].vegetat_type_num )
	      out_data[OUT_SNOW_CANOPY].data[0]
	      += (snow[veg][band].snow_canopy) 
	      * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack melt **/
	  out_data[OUT_SNOW_MELT].data[0]
	    += snow[veg][band].melt * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow cover fraction **/
          if ( snow[veg][band].snow && overstory ) {
            if ( snow[veg][band].snow_canopy > 0 )
              out_data[OUT_SNOW_COVER].data[0] += 1. * Cv * AreaFract[band] * TreeAdjustFactor[band];
          }
          else
            out_data[OUT_SNOW_COVER].data[0]
              += snow[veg][band].coverage * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack cold content **/
	  out_data[OUT_DELTACC].data[0]
	    += energy[veg][band].deltaCC * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack advection **/
          if ( snow[veg][band].snow && overstory )
            out_data[OUT_ADVECTION].data[0]
              += energy[veg][band].canopy_advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  out_data[OUT_ADVECTION].data[0]
	    += energy[veg][band].advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow energy flux **/
	  out_data[OUT_SNOW_FLUX].data[0]
	    += energy[veg][band].snow_flux * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record refreeze energy **/
          if ( snow[veg][band].snow && overstory ) {
	    out_data[OUT_RFRZ_ENERGY].data[0]
	      += energy[veg][band].canopy_refreeze 
	      * Cv * AreaFract[band] * TreeAdjustFactor[band];
          }
	  out_data[OUT_RFRZ_ENERGY].data[0]
	    += energy[veg][band].refreeze_energy 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

          /** record melt energy **/
          out_data[OUT_MELT_ENERGY].data[0]
            += energy[veg][band].melt_energy
            * Cv * AreaFract[band] * TreeAdjustFactor[band];

          /** record advected sensible heat energy **/
          if ( !overstory )
            out_data[OUT_ADV_SENS].data[0]
              -= energy[veg][band].advected_sensible
              * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** if snow elevation bands are to be printed separately **/
	  if(options.PRT_SNOW_BAND) {
	    
	    /** record band snow water equivalent **/
	    out_data[OUT_SWE_BAND].data[band]
	      += snow[veg][band].swq * Cv * 1000.;
	    
	    /** record band snowpack depth **/
	    out_data[OUT_SNOW_DEPTH_BAND].data[band]
	      += snow[veg][band].depth * Cv * 100.;
	    
	    /** record band canopy intercepted snow **/
	    if ( veg < veg_con[0].vegetat_type_num )
		out_data[OUT_SNOW_CANOPY_BAND].data[band]
		+= (snow[veg][band].snow_canopy) * Cv * 1000.;
	    
	    /** record band snow melt **/
	    out_data[OUT_SNOW_MELT_BAND].data[band]
	      += snow[veg][band].melt * Cv;
	    
	    /** record band snow coverage **/
	    out_data[OUT_SNOW_COVER_BAND].data[band]
	      += snow[veg][band].coverage * Cv;
	    
	    /** record band cold content **/
	    out_data[OUT_DELTACC_BAND].data[band]
	      += energy[veg][band].deltaCC * Cv;
	    
	    /** record band advection **/
	    out_data[OUT_ADVECTION_BAND].data[band]
	      += energy[veg][band].advection * Cv;
	    
	    /** record band snow flux **/
	    out_data[OUT_SNOW_FLUX_BAND].data[band]
	      += energy[veg][band].snow_flux * Cv;
	    
	    /** record band refreeze energy **/
	    out_data[OUT_RFRZ_ENERGY_BAND].data[band]
	      += energy[veg][band].refreeze_energy * Cv;
	    
            /** record band melt energy **/
            out_data[OUT_MELT_ENERGY_BAND].data[band]
              += energy[veg][band].melt_energy * Cv;

            /** record band advected sensble heat **/
            out_data[OUT_ADV_SENS_BAND].data[band]
              -= energy[veg][band].advected_sensible * Cv;

            /** record surface layer temperature **/
            out_data[OUT_SNOW_SURFT_BAND].data[band]
              += snow[veg][band].surf_temp * Cv;

            /** record pack layer temperature **/
            out_data[OUT_SNOW_PACKT_BAND].data[band]
              += snow[veg][band].pack_temp * Cv;

            /** record latent heat of sublimation **/
            out_data[OUT_LATENT_SUB_BAND].data[band]
              += energy[veg][band].latent_sub * Cv;

	    /** record band net downwards shortwave radiation **/
	    out_data[OUT_NET_SHORT_BAND].data[band]
	      += energy[veg][band].NetShortAtmos * Cv;

	    /** record band net downwards longwave radiation **/
	    out_data[OUT_NET_LONG_BAND].data[band]
	      += energy[veg][band].NetLongAtmos * Cv;

	    /** record band albedo **/
            if ( snow[veg][band].snow && overstory )
              out_data[OUT_ALBEDO_BAND].data[band]
	        += energy[veg][band].AlbedoOver * Cv;
            else
              out_data[OUT_ALBEDO_BAND].data[band]
	        += energy[veg][band].AlbedoUnder * Cv;

	    /** record band net latent heat flux **/
	    out_data[OUT_LATENT_BAND].data[band]
	      -= energy[veg][band].latent * Cv;

	    /** record band net sensible heat flux **/
	    out_data[OUT_SENSIBLE_BAND].data[band]
	      -= energy[veg][band].sensible * Cv;

	    /** record band net ground heat flux **/
	    out_data[OUT_GRND_FLUX_BAND].data[band]
	      -= energy[veg][band].grnd_flux * Cv;

	  }
	}
      }

    }
  }
 
  /********************
    Lake Model Output
  ********************/

  if ( options.LAKES ) {
    if (lake_con->Cl[0] > 0 ) {

      /** If lake fraction exists store energy and water fluxes */

      // Fraction of wetland with open water
      Clake = lake_var.sarea/lake_con->basin[0];

      // Fraction of grid cell that is lakes and wetlands.
      Cv = lake_con->Cl[0];
      mu = 1;
      band = 0;
      veg = veg_con[0].vegetat_type_num + 1;

      tmp_evap = lake_var.evapw * Clake * Cv;

      for ( index = 0; index < options.Nlayer; index++ )
        tmp_evap += cell[0][veg][band].layer[index].evap*(1.-Clake)*Cv;
      tmp_evap += veg_var[0][veg][band].canopyevap*(1.-Clake)*Cv;

      tmp_evap += snow[veg][band].vapor_flux * 1000. * Cv;

      out_data[OUT_EVAP].data[0] += tmp_evap * AreaFract[band] * TreeAdjustFactor[band];

      out_data[OUT_SUB_SNOW].data[0] += snow[veg][band].vapor_flux * 1000. * Cv * AreaFract[band] * TreeAdjustFactor[band];
      out_data[OUT_SUB_SURFACE].data[0] += snow[veg][band].surface_flux * 1000. * Cv * AreaFract[band] * TreeAdjustFactor[band];
      out_data[OUT_SUB_BLOWING].data[0] += snow[veg][band].blowing_flux * 1000. * Cv * AreaFract[band] * TreeAdjustFactor[band];
      out_data[OUT_EVAP_LAKE].data[0] = lake_var.evapw * Clake * Cv * AreaFract[band] * TreeAdjustFactor[band]; // mm over gridcell
      out_data[OUT_SNOW_MELT].data[0] += snow[veg][band].melt * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record runoff **/
      out_data[OUT_RUNOFF].data[0]   = lake_var.runoff_out;

      /** record baseflow **/
      out_data[OUT_BASEFLOW].data[0] = lake_var.baseflow_out;

      /** record freezing and thawing front depths **/
      if(options.FROZEN_SOIL) {
        if(Clake != 1.) {
          for(index = 0; index < MAX_FRONTS; index++) {
            if(energy[veg][band].fdepth[index] != MISSING)
              out_data[OUT_FDEPTH].data[index] += energy[veg][band].fdepth[index]
                * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
            if(energy[veg][band].tdepth[index] != MISSING)
              out_data[OUT_TDEPTH].data[index] += energy[veg][band].tdepth[index]
                * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];
          }
        }
      }

      /** record aerodynamic conductivity **/
      if (lake_var.aero_resist_used > SMALL) {
        out_data[OUT_AERO_COND].data[0] += (1/lake_var.aero_resist_used) * Clake * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
      }
      else {
        out_data[OUT_AERO_COND].data[0] = HUGE_RESIST;
        out_data[OUT_AERO_RESIST].data[0] = lake_var.aero_resist_used;
      }
      if (cell[WET][veg][0].aero_resist_used > SMALL) {
        out_data[OUT_AERO_COND].data[0] += (1/cell[WET][veg][0].aero_resist_used)
              * Cv * mu * (1.-Clake) * AreaFract[band] * TreeAdjustFactor[band];
      }
      else {
        out_data[OUT_AERO_COND].data[0] = HUGE_RESIST;
        out_data[OUT_AERO_RESIST].data[0] = cell[WET][veg][0].aero_resist_used;
      }

      /** record lake moistures **/
      out_data[OUT_LAKE_MOIST].data[0] = (lake_var.volume / lake_con->basin[0]) * 1000. * Clake * Cv; // mm over gridcell
      out_data[OUT_LAKE_ICE].data[0]   = ( ice_density * lake_var.fraci
                               * lake_var.hice / RHO_W );

      for ( index = 0; index < options.Nlayer; index++ ) {
        tmp_moist = cell[0][veg][band].layer[index].moist;
#if SPATIAL_FROST
        tmp_ice = 0;
        for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          tmp_ice  += (cell[0][veg][band].layer[index].ice[frost_area]
                       * frost_fract[frost_area]);
#else
        tmp_ice   = cell[0][veg][band].layer[index].ice;
#endif
        tmp_moist -= tmp_ice;
        if(options.MOISTFRACT) {
          tmp_moist /= depth[index] * 1000.;
          tmp_ice /= depth[index] * 1000.;
        }
        out_data[OUT_SOIL_LIQ].data[index] += tmp_moist * Cv * AreaFract[band] * TreeAdjustFactor[band];
        out_data[OUT_SOIL_ICE].data[index] += tmp_ice * Cv * AreaFract[band] * TreeAdjustFactor[band];
      }
      out_data[OUT_SOIL_WET].data[0] += cell[0][veg][band].wetness * Cv * AreaFract[band] * TreeAdjustFactor[band];
      out_data[OUT_ROOTMOIST].data[0] += cell[0][veg][band].rootmoist * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record layer temperatures **/
      for(index=0;index<options.Nlayer;index++) {
        out_data[OUT_SOIL_TEMP].data[index] += cell[0][veg][band].layer[index].T * Cv * AreaFract[band] * TreeAdjustFactor[band];
      }

      /** record thermal node temperatures **/
      for(index=0;index<options.Nnode;index++) {
        out_data[OUT_SOIL_TNODE].data[index] += energy[veg][band].T[index] * Cv * AreaFract[band] * TreeAdjustFactor[band];
      }

      /***************************************
        Record Lake Energy Balance Variables
      ***************************************/

      /** record surface radiative temperature **/
      rad_temp = energy[veg][band].Tsurf + KELVIN;

      /** record lake surface temperature **/
      surf_temp = lake_var.temp[0];

      /** record landcover temperature **/
      if(veg == veg_con[0].vegetat_type_num) {
        // landcover is bare soil
        out_data[OUT_BARESOILT].data[0] += (rad_temp-KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
        cv_baresoil += AreaFract[band] * TreeAdjustFactor[band];
      }
      else {
        // landcover is vegetation
        out_data[OUT_VEGT].data[0] += (rad_temp-KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
        cv_veg += Cv * AreaFract[band] * TreeAdjustFactor[band];
      }

      /** record net shortwave radiation **/
      out_data[OUT_NET_SHORT].data[0] += energy[veg][band].NetShortAtmos * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record net longwave radiation **/
      out_data[OUT_NET_LONG].data[0] += energy[veg][band].NetLongAtmos * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record incoming longwave radiation at ground surface (under veg) **/
      out_data[OUT_IN_LONG].data[0] += (energy[veg][band].NetLongAtmos + STEFAN_B
                              * (rad_temp) * (rad_temp) * (rad_temp) * (rad_temp))
                              * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record albedo **/
      out_data[OUT_ALBEDO].data[0] += (energy[veg][band].AlbedoLake*Clake
                              + energy[veg][band].AlbedoUnder*(1.-Clake))* Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record latent heat flux **/
      out_data[OUT_LATENT].data[0] -= energy[veg][band].AtmosLatent * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record latent heat flux **/
      out_data[OUT_LATENT_SUB].data[0] -= energy[veg][band].AtmosLatentSub * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record sensible heat flux **/
      out_data[OUT_SENSIBLE].data[0] -= ( energy[veg][band].AtmosSensible
                               + energy[veg][band].snow_flux ) * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record ground heat flux (+ heat storage) **/
      out_data[OUT_GRND_FLUX].data[0] -= (energy[veg][band].grnd_flux) * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record heat storage **/
      out_data[OUT_DELTAH].data[0] -= energy[veg][band].deltaH * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record energy balance error **/
      out_data[OUT_ENERGY_ERROR].data[0] += energy[veg][band].error * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record radiative effective temperature [K],
          emissivities set = 1.0  **/
      out_data[OUT_RAD_TEMP].data[0] += ((rad_temp) * (rad_temp)
                             * (rad_temp) * (rad_temp)) * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record mean surface temperature [C]  **/
      out_data[OUT_SURF_TEMP].data[0] += surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /**********************************
        Record Lake Snow Pack Variables
      **********************************/

      /** record snow water equivalence **/
      out_data[OUT_SWE].data[0] += snow[veg][band].swq * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];

      /** record snowpack depth **/
      out_data[OUT_SNOW_DEPTH].data[0] += snow[veg][band].depth * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];

      /** record snow cover fraction **/
      out_data[OUT_SNOW_COVER].data[0] += snow[veg][band].coverage * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record snowpack albedo, temperature **/
      if (snow[veg][band].swq> 0.0) {
        out_data[OUT_SALBEDO].data[0]
          += snow[veg][band].albedo * Cv * AreaFract[band] * TreeAdjustFactor[band];
        out_data[OUT_SNOW_SURF_TEMP].data[0]
          += snow[veg][band].surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
        out_data[OUT_SNOW_PACK_TEMP].data[0]
          += snow[veg][band].pack_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
        cv_snow += Cv * AreaFract[band] * TreeAdjustFactor[band];
      }

      /** record snowpack cold content **/
      out_data[OUT_DELTACC].data[0] += energy[veg][band].deltaCC * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record snowpack advection **/
      out_data[OUT_ADVECTION].data[0] += energy[veg][band].advection * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record snow energy flux **/
      out_data[OUT_SNOW_FLUX].data[0] += energy[veg][band].snow_flux * Cv * AreaFract[band] * TreeAdjustFactor[band];

      /** record refreeze energy **/
      out_data[OUT_RFRZ_ENERGY].data[0] += energy[veg][band].refreeze_energy * Cv * AreaFract[band] * TreeAdjustFactor[band];

      // Store Lake Specific Variables
      out_data[OUT_LAKE_ICE_TEMP].data[0]   = lake_var.tempi;
      out_data[OUT_LAKE_ICE_HEIGHT].data[0] = lake_var.hice;
      out_data[OUT_LAKE_ICE_FRACT].data[0]  = lake_var.fraci*Cv;
      out_data[OUT_LAKE_DEPTH].data[0]      = lake_var.ldepth;
      //out_data[OUT_LAKE_SURF_AREA].data[0]  = Clake;
      out_data[OUT_LAKE_SURF_AREA].data[0]  = lake_var.sarea;
      out_data[OUT_LAKE_VOLUME].data[0]     = lake_var.volume;
      out_data[OUT_LAKE_SURF_TEMP].data[0]  = lake_var.temp[0];
      out_data[OUT_SURFSTOR].data[0]        = (lake_var.volume/lake_var.sarea) * 1000. * Clake * Cv * AreaFract[band] * TreeAdjustFactor[band];

    }
  }

  /*****************************************
    Finish aggregation of special-case variables
   *****************************************/
  // Normalize quantities that aren't present over entire grid cell
  if (cv_baresoil > 0) {
    out_data[OUT_BARESOILT].data[0] /= cv_baresoil;
  }
  if (cv_veg > 0) {
    out_data[OUT_VEGT].data[0] /= cv_veg;
  }
  if (cv_snow > 0) {
    out_data[OUT_SALBEDO].data[0] /= cv_snow;
    out_data[OUT_SNOW_SURF_TEMP].data[0] /= cv_snow;
    out_data[OUT_SNOW_PACK_TEMP].data[0] /= cv_snow;
  }

  // Radiative temperature
  out_data[OUT_RAD_TEMP].data[0] = pow(out_data[OUT_RAD_TEMP].data[0],0.25);

  // Aerodynamic conductance and resistance
  if (out_data[OUT_AERO_COND].data[0] < HUGE_RESIST) {
    out_data[OUT_AERO_RESIST].data[0] = 1 / out_data[OUT_AERO_COND].data[0];
  }

  /*****************************************
    Compute derived variables
   *****************************************/
  // Water balance terms
  out_data[OUT_DELSOILMOIST].data[0] = 0;
  for (index=0; index<options.Nlayer; index++) {
    out_data[OUT_SOIL_MOIST].data[index] = out_data[OUT_SOIL_LIQ].data[index]+out_data[OUT_SOIL_ICE].data[index];
    out_data[OUT_DELSOILMOIST].data[0] += out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMLIQFRAC].data[index] = out_data[OUT_SOIL_LIQ].data[index]/out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMFROZFRAC].data[index] = 1 - out_data[OUT_SMLIQFRAC].data[index];
  }
  out_data[OUT_DELSOILMOIST].data[0] -= save_data->total_soil_moist;
  out_data[OUT_DELSWE].data[0] = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] - save_data->swe;
  out_data[OUT_DELINTERCEPT].data[0] = out_data[OUT_WDEW].data[0] - save_data->wdew;

  // Energy terms
  out_data[OUT_REFREEZE].data[0] = (out_data[OUT_RFRZ_ENERGY].data[0]/Lf)*dt_sec;
  out_data[OUT_R_NET].data[0] = out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0];

  // Save current moisture state for use in next time step
  save_data->total_soil_moist = 0;
  for (index=0; index<options.Nlayer; index++) {
    save_data->total_soil_moist += out_data[OUT_SOIL_MOIST].data[index];
  }
  save_data->surfstor = out_data[OUT_SURFSTOR].data[0];
  save_data->swe = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0];
  save_data->wdew = out_data[OUT_WDEW].data[0];

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data[OUT_PREC].data[0];
  outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0];
  storage = 0.;
  for(index=0;index<options.Nlayer;index++)
    if(options.MOISTFRACT)
      storage += (out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index]) 
	* depth[index] * 1000;
    else
      storage += out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index];
  storage += out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] + out_data[OUT_WDEW].data[0] + out_data[OUT_SURFSTOR].data[0];
  calc_water_balance_error(rec,inflow,outflow,storage);

  /********************
    Check Energy Balance 
    ********************/
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0],
			      out_data[OUT_LATENT].data[0]+out_data[OUT_LATENT_SUB].data[0],
			      out_data[OUT_SENSIBLE].data[0]+out_data[OUT_ADV_SENS].data[0],
			      out_data[OUT_GRND_FLUX].data[0]+out_data[OUT_DELTAH].data[0]+out_data[OUT_FUSION].data[0],
			      out_data[OUT_ADVECTION].data[0] - out_data[OUT_DELTACC].data[0]
			      - out_data[OUT_SNOW_FLUX].data[0] + out_data[OUT_RFRZ_ENERGY].data[0]);

  /********************
    Temporal Aggregation 
    ********************/
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    if (out_data[v].aggtype == AGG_TYPE_END) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_SUM) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_AVG) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i]/out_step_ratio;
      }
    }
  }

  /********************
    Output procedure
    (only execute when we've completed an output interval)
    ********************/
  if (step_count == out_step_ratio) {

    /***********************************************
      Change of units for ALMA-compliant output
    ***********************************************/
    if (options.ALMA_OUTPUT) {
      out_data[OUT_BASEFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_BARE].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_LAKE].aggdata[0] /= out_dt_sec;
      out_data[OUT_INFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_PREC].aggdata[0] /= out_dt_sec;
      out_data[OUT_RAINF].aggdata[0] /= out_dt_sec;
      out_data[OUT_REFREEZE].aggdata[0] /= out_dt_sec;
      out_data[OUT_RUNOFF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOW_MELT].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOWF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_BLOWING].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] += out_data[OUT_SUB_CANOP].aggdata[0];
      out_data[OUT_SUB_SURFACE].aggdata[0] /= out_dt_sec;
      out_data[OUT_TRANSP_VEG].aggdata[0] /= out_dt_sec;
      out_data[OUT_BARESOILT].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_PACK_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_SURF_TEMP].aggdata[0] += KELVIN;
      for (index=0; index<options.Nlayer; index++) {
        out_data[OUT_SOIL_TEMP].aggdata[index] += KELVIN;
      }
      for (index=0; index<options.Nnode; index++) {
        out_data[OUT_SOIL_TNODE].aggdata[index] += KELVIN;
      }
      out_data[OUT_SURF_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_VEGT].aggdata[0] += KELVIN;
      out_data[OUT_FDEPTH].aggdata[0] /= 100;
      out_data[OUT_TDEPTH].aggdata[0] /= 100;
      out_data[OUT_DELTACC].aggdata[0] *= out_dt_sec;
      out_data[OUT_DELTAH].aggdata[0] *= out_dt_sec;
      out_data[OUT_AIR_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_PRESSURE].aggdata[0] *= 1000;
      out_data[OUT_VP].aggdata[0] *= 1000;
    }

    /*************
      Write Data
    *************/
    if(rec >= skipyear) {
      if (options.BINARY_OUTPUT) {
        for (v=0; v<N_OUTVAR_TYPES; v++) {
          for (i=0; i<out_data[v].nelem; i++) {
            out_data[v].aggdata[i] *= out_data[v].mult;
          }
        }
      }
      write_data(out_data_files, out_data, dmy, global_param.out_dt);
    }

    // Reset the step count
    step_count = 0;

    // Reset the agg data
    for (v=0; v<N_OUTVAR_TYPES; v++) {
      for (i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = 0;
      }
    }

  } // End of output procedure

}
