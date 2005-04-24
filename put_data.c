#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void put_data(soil_con_struct	*soil_con,
	      global_param_struct *global_param,
	      int                Nnodes,
              int                rec,
	      atmos_data_struct *atmos,
	      dist_prcp_struct  *prcp,
	      dmy_struct        *dmy,
#if LAKE_MODEL
	      lake_con_struct   *lake_con,
#endif // LAKE_MODEL
              outfiles_struct   *outfiles,
	      veg_con_struct    *veg_con)
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
  30-Oct-03 Snow_flux was incorrectly set to Tcanopy.  Fixed.	TJB
  25-Aug-04 Sub_snow was incorrectly set to blowing_flux.  Now it is
	    set to vapor_flux.					TJB
  28-Sep-04 Now out_data->aero_resist stores the aerodynamic resistance
	    used in flux calculations.				TJB
  2005-Mar-24 Modified to compute ALMA output variables.	TJB
  2005-Apr-23 Now aero_cond is aggregated instead of aero_resist.	TJB

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
  static double           snow_depth;
#endif

  out_data_struct        *out_data;
  out_data_alma_struct   *out_data_alma;

  int                     veg;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  int n;
  int                     overstory;
#if SPATIAL_FROST
  int                     frost_area;
#endif
  char              *AboveTreeLine;
  double            *porosity;
  double            *Wpwp;
  double            *AreaFract;
  double            *depth;
  double            *dz;
#if SPATIAL_FROST
  double            *frost_fract;
  double             frost_slope;
#endif // SPATIAL_FROST
  double             dp;
  int                dt;
  int                skipyear;
  double             MIN_RAIN_TEMP;
  double             MAX_SNOW_TEMP;
  double                  Cv;
  double                  Clake;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
  double tmp_total_moist;
  double tmp_root_moist;
  double tmp_total_soil;
  double tmp_temp;
  double cv_baresoil;
  double cv_foliage;
  double cv_snow;
  double                  tmp_ice;
  double                  tmp_fract;
  double                  rad_temp;
  double                  surf_temp;
  double                  inflow;
  double                  outflow;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];

  cell_data_struct     ***cell;
  energy_bal_struct     **energy;
#if LAKE_MODEL
  lake_var_struct         lake_var;
#endif // LAKE_MODEL
  snow_data_struct      **snow;
  veg_var_struct       ***veg_var;

  AboveTreeLine = soil_con->AboveTreeLine;
  porosity = soil_con->porosity;
  Wpwp = soil_con->Wpwp;
  AreaFract = soil_con->AreaFract;
  depth = soil_con->depth;
  dz = soil_con->dz_node;
#if SPATIAL_FROST
  frost_fract = soil_con->frost_fract;
  frost_slope = soil_con->frost_slope;
#endif // SPATIAL_FROST
  dp = soil_con->dp;
  dt = global_param->dt;
  skipyear = global_param->skipyear;
  MIN_RAIN_TEMP = global_param->MIN_RAIN_TEMP;
  MAX_SNOW_TEMP = global_param->MAX_SNOW_TEMP;

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
  out_data->rel_humid      = ( 100. * atmos->vp[NR] 
			       / (atmos->vp[NR] + atmos->vpd[NR]) );
  if(out_data->air_temp > (MAX_SNOW_TEMP+MIN_RAIN_TEMP)/2) {
    out_data->rain = out_data->prec / (dt*3600);
    out_data->snow = 0.;
  }
  else {
    out_data->rain = 0.;
    out_data->snow = out_data->prec / (dt*3600);
  }
  out_data->surf_temp      = 0.;
  out_data->snow_albedo = SPVAL;
  for(band=0;band<Nbands;band++) {
    out_data->snow_surf_temp[band] = SPVAL;
    out_data->snow_pack_temp[band] = SPVAL;
  }
  out_data->baresoilt = SPVAL;
  out_data->tfoliage = SPVAL;
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp->cell;
  energy  = prcp->energy;
#if LAKE_MODEL
  lake_var = prcp->lake_var;
#endif // LAKE_MODEL
  snow    = prcp->snow;
  veg_var = prcp->veg_var;
  
  /****************************************
    Store Output for all Vegetation Types (except lakes)
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {
    
    if ( veg < veg_con[0].vegetat_type_num ) {
      Cv = veg_con[veg].Cv;
      overstory = veg_lib[veg_con[veg].veg_class].overstory;
    }
    else {
      Cv = (1.0 - veg_con[0].Cv_sum);
      overstory = FALSE;
    }

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
	    for ( index = 0; index < options.Nlayer; index++ )
	      tmp_evap += cell[dist][veg][band].layer[index].evap;
	    if ( veg < veg_con[0].vegetat_type_num )
	      out_data->evap_veg += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    else 
	      out_data->evap_bare += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];

	    tmp_evap += snow[veg][band].vapor_flux * 1000.;
	    out_data->sub_total += snow[veg][band].vapor_flux * 1000.
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    out_data->sub_surface += snow[veg][band].surface_flux * 1000.
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    out_data->sub_blowing += snow[veg][band].blowing_flux * 1000.
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    out_data->sub_snow[veg] += snow[veg][band].vapor_flux * 1000.; 
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
	  
	    /** record aerodynamic conductance **/
            if (cell[WET][veg][0].aero_resist_used > SMALL) {
	      out_data->aero_cond += (1/cell[WET][veg][0].aero_resist_used)
	        * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }
            else {
              out_data->aero_cond = HUGE_RESIST;
              out_data->aero_resist = cell[WET][veg][0].aero_resist_used;
            }
	    
	    /** record layer moistures **/
	    tmp_total_moist = tmp_total_soil = tmp_root_moist = 0.0;
	    for ( index = 0; index < options.Nlayer; index++ ) {
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
	      tmp_moist *= Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      tmp_ice *= Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      tmp_total_moist += tmp_moist + tmp_ice - Wpwp[index];
	      tmp_total_soil += depth[index] * porosity[index] * 1000. - Wpwp[index];
              if (tmp_total_moist - tmp_ice > 0.0)
	        tmp_root_moist += tmp_total_moist - tmp_ice;
	      out_data->moist[index] += tmp_moist;
	      out_data->ice[index]   += tmp_ice;
	    }
	    out_data->soil_wetness += tmp_total_moist/tmp_total_soil;
	    out_data->rootmoist += tmp_root_moist;
#endif

	    /** record layer temperatureses **/
	    for ( index = 0; index < options.Nlayer; index++ ) {
	      tmp_temp = cell[dist][veg][band].layer[index].T;

	      tmp_temp *= Cv * AreaFract[band];
	      out_data->layer_temp[index] += tmp_temp;
	    }

	  }
	}
      }

      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0.&& ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !veg_lib[veg_con[veg].veg_class].overstory)))) {

#if !OPTIMIZE

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
	  
#if SPATIAL_FROST
	  tmp_fract = 0;
	  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	    if ( cell[0][veg][band].layer[0].ice[frost_area] )
	      tmp_fract  += frost_fract[frost_area]; 
#else
  	  if ( cell[0][veg][band].layer[0].ice > 0 )
	    tmp_fract   = 1.;
#endif
	  out_data->surf_frost_fract += tmp_fract * Cv * AreaFract[band] * TreeAdjustFactor[band];

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
          /** note: these are undefined if the landcover class is not present **/
	  if(veg == veg_con[0].vegetat_type_num) {
            // landcover is bare soil
            if (out_data->baresoilt > SPVAL/2) {
              // first assignment of baresoilt
	      out_data->baresoilt = rad_temp * AreaFract[band] * TreeAdjustFactor[band];
            }
            else {
              // not first assignment of baresoilt
	      out_data->baresoilt += rad_temp * AreaFract[band] * TreeAdjustFactor[band];
            }
            cv_baresoil += AreaFract[band] * TreeAdjustFactor[band];
          }
          else {
            // landcover is vegetation
            if (overstory && !snow[veg][band].snow) {
              // here, rad_temp will be wrong since it will pick the understory temperature
              // so we'll use Tfoliage instead
              if (out_data->tfoliage > SPVAL/2) {
                // first assignment of tfoliage
	        out_data->tfoliage = (energy[veg][band].Tfoliage + KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
              }
              else {
                // not first assignment of tfoliage
	        out_data->tfoliage += (energy[veg][band].Tfoliage + KELVIN) * Cv * AreaFract[band] * TreeAdjustFactor[band];
              }
            }
            else {
              // use rad_temp in all other cases
              if (out_data->tfoliage > SPVAL/2) {
                // first assignment of tfoliage
	        out_data->tfoliage = rad_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
              }
              else {
                // not first assignment of tfoliage
	        out_data->tfoliage += rad_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
              }
            }
            cv_foliage += Cv * AreaFract[band] * TreeAdjustFactor[band];
          }

	  /** record soil surface temperature **/
	  surf_temp = energy[veg][band].T[0];
	  
	  /** record net shortwave radiation **/
	  out_data->net_short += energy[veg][band].NetShortAtmos
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record net longwave radiation **/
	  out_data->net_long  += energy[veg][band].NetLongAtmos
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record incoming longwave radiation **/
	  if ( snow[veg][band].snow && overstory )
	    out_data->in_long  += (energy[veg][band].LongOverIn 
				   * Cv * AreaFract[band] * TreeAdjustFactor[band]);
	  else
	    out_data->in_long  += (energy[veg][band].LongUnderIn 
				   * Cv * AreaFract[band] * TreeAdjustFactor[band]);
	  
	  /** record albedo **/
	  if ( snow[veg][band].snow && overstory )
	    out_data->albedo    += energy[veg][band].AlbedoOver
	      * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  else
	    out_data->albedo    += energy[veg][band].AlbedoUnder
	      * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record latent heat flux **/
	  out_data->latent    -= energy[veg][band].AtmosLatent
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record latent heat flux from sublimation **/
	  out_data->latent_sub[0] 
	    -= energy[veg][band].AtmosLatentSub * Cv * AreaFract[band] * TreeAdjustFactor[band];
	    
	  /** record sensible heat flux **/
	  out_data->sensible  -= energy[veg][band].AtmosSensible
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record ground heat flux (+ heat storage) **/
	  out_data->grnd_flux -= (energy[veg][band].grnd_flux)
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record heat storage **/
	  out_data->deltaH    -= energy[veg][band].deltaH
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record heat of fusion **/
	  out_data->fusion    -= energy[veg][band].fusion
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record energy balance error **/
	  out_data->energy_error += energy[veg][band].error
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record radiative effective temperature [K], 
	      emissivities set = 1.0  **/
	  out_data->rad_temp += ((rad_temp) * (rad_temp) 
				 * (rad_temp) * (rad_temp)) 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  out_data->AvgSurfT += rad_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
  
	  /** record mean surface temperature [C]  **/
	  /*out_data->surf_temp += surf_temp * Cv * AreaFract[band];*/
	  out_data->surf_temp += (energy[veg][band].T[0] 
				  + energy[veg][band].T[1]) / 2 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /*****************************
	    Record Snow Pack Variables 
	  *****************************/
	  
	  /** record snow water equivalence **/
	  out_data->swq[0]         
	    += snow[veg][band].swq * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
#endif // !OPTIMIZE
	  
	  /** record snowpack depth **/
	  out_data->snow_depth[0]  
	    += snow[veg][band].depth * Cv * 100. * AreaFract[band] * TreeAdjustFactor[band];

          /** record snowpack albedo, temperature **/
          /** Note: these quantities are undefined if swq is 0 **/
          if (out_data->swq[0] > 0.0) {
            if (out_data->snow_albedo >= SPVAL/2) {
	      out_data->snow_albedo
	        = snow[veg][band].albedo * Cv * AreaFract[band] * TreeAdjustFactor[band];
	      out_data->snow_surf_temp[0]   
	        = snow[veg][band].surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
	      out_data->snow_pack_temp[0]   
	        = snow[veg][band].pack_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
            }
            else {
	      out_data->snow_albedo
  	        += snow[veg][band].albedo * Cv * AreaFract[band] * TreeAdjustFactor[band];
  	      out_data->snow_surf_temp[0]   
  	        += snow[veg][band].surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
  	      out_data->snow_pack_temp[0]   
	        += snow[veg][band].pack_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
            }
            cv_snow += Cv;
          }

#if !OPTIMIZE

	  /** record canopy intercepted snow **/
	  if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	    out_data->swq[0] 
#else
	      out_data->snow_canopy[0] 
#endif
	      += (snow[veg][band].snow_canopy) 
	      * Cv * 1000. * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack melt **/
	  out_data->melt[0]  
	    += snow[veg][band].melt * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow cover fraction **/
	  if ( snow[veg][band].snow && overstory ) {
	    if ( snow[veg][band].snow_canopy > 0 )
	      out_data->coverage[0] += 1. * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  }
	  else
	    out_data->coverage[0]    
	      += snow[veg][band].coverage * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack cold content **/
	  out_data->deltaCC[0]           
	    += energy[veg][band].deltaCC * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack advection **/
	  if ( snow[veg][band].snow && overstory )
	    out_data->advection[0]         
	      += energy[veg][band].canopy_advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  out_data->advection[0]         
	    += energy[veg][band].advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow energy flux **/
	  out_data->snow_flux[0]         
 	    += energy[veg][band].snow_flux * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  /*  += energy[veg][band].Tcanopy * Cv * AreaFract[band] * TreeAdjustFactor[band]; */
	  
	  /** record refreeze energy **/
	  if ( snow[veg][band].snow && overstory ) {
	    out_data->refreeze_energy[0]   
	      += energy[veg][band].canopy_refreeze 
	      * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  }
	  out_data->refreeze_energy[0]   
	    += energy[veg][band].refreeze_energy 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record melt energy **/
	  out_data->melt_energy[0]   
	    += energy[veg][band].melt_energy 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record advected sensible heat energy **/
	  if ( !overstory )
	    out_data->advected_sensible[0]   
	      -= energy[veg][band].advected_sensible 
	      * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** if snow elevation bands are to be printed separately **/
	  if(options.PRT_SNOW_BAND) {
	    
	    /** record band snow water equivalent **/
	    out_data->swq[band+1]         
	      += snow[veg][band].swq * Cv  * 1000. * TreeAdjustFactor[band];
	    
	    /** record band snowpack depth **/
	    out_data->snow_depth[band+1]  
	      += snow[veg][band].depth * Cv * 100. * TreeAdjustFactor[band];
	    
	    /** record band canopy intercepted snow **/
	    if ( veg < veg_con[0].vegetat_type_num )
#if LDAS_OUTPUT
	      out_data->swq[band+1]
#else
		out_data->snow_canopy[band+1]  
#endif
		+= (snow[veg][band].snow_canopy) * Cv * 1000. * TreeAdjustFactor[band];
	    
	    /** record band snowpack melt **/
	    out_data->melt[band+1]  
	      += snow[veg][band].melt * Cv * TreeAdjustFactor[band];
	    
	    /** record band snow coverage **/
	    out_data->coverage[band+1]    
	      += snow[veg][band].coverage * Cv * TreeAdjustFactor[band];
	    
	    /** record band cold content **/
	    out_data->deltaCC[band+1]           
	      += energy[veg][band].deltaCC * Cv * TreeAdjustFactor[band];
	    
	    /** record band advection **/
	    out_data->advection[band+1]         
	      += energy[veg][band].advection * Cv * TreeAdjustFactor[band];
	    
	    /** record band snow flux **/
	    out_data->snow_flux[band+1]         
	      += energy[veg][band].snow_flux * Cv * TreeAdjustFactor[band];
	    
	    /** record band refreeze energy **/
	    out_data->refreeze_energy[band+1]   
	      += energy[veg][band].refreeze_energy * Cv * TreeAdjustFactor[band];
	    
	    /** record band melt energy **/
	    out_data->melt_energy[band+1]   
	      += energy[veg][band].melt_energy * Cv * TreeAdjustFactor[band];
	    
	    /** record band advected sensble heat **/
	    out_data->advected_sensible[band+1]   
	      -= energy[veg][band].advected_sensible * Cv * TreeAdjustFactor[band];
	    
	    /** record surface layer temperature **/
	    out_data->snow_surf_temp[band+1]   
	      += snow[veg][band].surf_temp * Cv * TreeAdjustFactor[band];
	    
	    /** record pack layer temperature **/
	    out_data->snow_pack_temp[band+1]   
	      += snow[veg][band].pack_temp * Cv * TreeAdjustFactor[band];
	    
	    /** record latent heat of sublimation **/
	    out_data->latent_sub[band+1]   
	      += energy[veg][band].latent_sub * Cv * TreeAdjustFactor[band];
	    
	  }

#endif /* not OPTIMIZE */

	}
      }
    }
  }

  // normalize quantities that may not be present over all of grid cell
  if (out_data->baresoilt < SPVAL/2) {
    out_data->baresoilt /= cv_baresoil;
  }
  if (out_data->tfoliage < SPVAL/2) {
    out_data->tfoliage /= cv_foliage;
  }
  if (out_data->snow_albedo < SPVAL/2) {
    out_data->snow_albedo /= cv_snow;
    out_data->snow_surf_temp[0] /= cv_snow;
    out_data->snow_pack_temp[0] /= cv_snow;
  }
  
#if LAKE_MODEL

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
      
      out_data->evap += tmp_evap; 

      out_data->sub_total += snow[veg][band].vapor_flux * 1000. * Cv; 
      out_data->sub_surface += snow[veg][band].surface_flux * 1000. * Cv; 
      out_data->sub_blowing += snow[veg][band].blowing_flux * 1000. * Cv; 
      out_data->sub_snow[veg] += snow[veg][band].vapor_flux * 1000.; 
      /*    out_data->evap_lake = lake_var.evapw * Clake * Cv; */  /* mm over gridcell */
      out_data->evap_lake = lake_var.evapw;
      out_data->melt[0] += snow[veg][band].melt * Cv;

      /** record runoff **/
      out_data->runoff   = lake_var.runoff_out;
	    
      /** record baseflow **/
      out_data->baseflow = lake_var.baseflow_out; 
	    
      /** record freezing and thawing front depths **/
      if(options.FROZEN_SOIL) {
	if(Clake != 1.) {
	  for(index = 0; index < MAX_FRONTS; index++) {
	    if(energy[veg][band].fdepth[index] != MISSING)
	      out_data->fdepth[index] += energy[veg][band].fdepth[index] 
		* Cv * 100. * AreaFract[band];
	    if(energy[veg][band].tdepth[index] != MISSING)
	      out_data->tdepth[index] += energy[veg][band].tdepth[index] 
		* Cv * 100. * AreaFract[band];
	  }
	}
      }

#if ! OPTIMIZE
	    
      /** record aerodynamic conductivity **/
      if (lake_var.aero_resist_used > SMALL) {
        out_data->aero_cond += (1/lake_var.aero_resist_used) * Clake * Cv * mu;
      }
      else {
        out_data->aero_cond = HUGE_RESIST;
        out_data->aero_resist = lake_var.aero_resist_used;
      }
      if (cell[WET][veg][0].aero_resist_used > SMALL) {
        out_data->aero_cond += (1/cell[WET][veg][0].aero_resist_used)
	      * Cv * mu * (1.-Clake);
      }
      else {
        out_data->aero_cond = HUGE_RESIST;
        out_data->aero_resist = cell[WET][veg][0].aero_resist_used;
      }

      /** record lake moistures **/
      out_data->lake_moist = lake_con->Cl[0] * 1000. * lake_var.volume / lake_con->basin[0]; // mm
      out_data->lake_ice   = ( ice_density * lake_var.fraci 
			       * lake_var.hice / RHO_W );

      tmp_total_moist = tmp_total_soil = tmp_root_moist = 0.0;
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
	tmp_moist *= Cv;
	tmp_ice *= Cv;
	tmp_total_moist += tmp_moist + tmp_ice - Wpwp[index];
	tmp_total_soil += depth[index] * porosity[index] * 1000. - Wpwp[index];
        if (tmp_total_moist - tmp_ice > 0.0)
	  tmp_root_moist += tmp_total_moist - tmp_ice;
	out_data->moist[index] += tmp_moist;
	out_data->ice[index]   += tmp_ice;
      }
      out_data->soil_wetness += tmp_total_moist/tmp_total_soil;
      out_data->rootmoist += tmp_root_moist;

      /** record layer temperatures **/
      for ( index = 0; index < options.Nlayer; index++ ) {
        tmp_temp = cell[0][veg][band].layer[index].T;
        tmp_temp *= Cv;
        out_data->layer_temp[index] += tmp_temp;
      }

      /***************************************
        Record Lake Energy Balance Variables
      ***************************************/

      /** record surface radiative temperature **/
      rad_temp = energy[veg][band].Tsurf + KELVIN;
      
      /** record lake surface temperature **/
      surf_temp = lake_var.temp[0];
      
      /** record net shortwave radiation **/
      out_data->net_short += energy[veg][band].NetShortAtmos * Cv;
      
      /** record net longwave radiation **/
      out_data->net_long  += energy[veg][band].NetLongAtmos * Cv;
      
      /** record incoming longwave radiation **/
      out_data->in_long  += ((energy[veg][band].NetLongAtmos + STEFAN_B 
			      * (rad_temp) * (rad_temp)
			      * (rad_temp) * (rad_temp)) * Cv);
      
      /** record albedo **/
      out_data->albedo    += (energy[veg][band].AlbedoLake*Clake  +
			      energy[veg][band].AlbedoUnder*(1.-Clake))* Cv;
      
      /** record latent heat flux **/
      out_data->latent        -= energy[veg][band].AtmosLatent * Cv;
      
      /** record latent heat flux **/
      out_data->latent_sub[0] -= energy[veg][band].AtmosLatentSub * Cv;
      
      /** record sensible heat flux **/
      out_data->sensible  -= ( energy[veg][band].AtmosSensible 
			       + energy[veg][band].snow_flux ) * Cv;
      
      /** record ground heat flux (+ heat storage) **/
      out_data->grnd_flux -= (energy[veg][band].grnd_flux) * Cv;
      
      /** record heat storage **/
      out_data->deltaH    -= energy[veg][band].deltaH * Cv;
      
      /** record energy balance error **/
      out_data->energy_error += energy[veg][band].error * Cv;
      
      /** record radiative effective temperature [K], 
	  emissivities set = 1.0  **/
      out_data->rad_temp += ((rad_temp) * (rad_temp) 
			     * (rad_temp) * (rad_temp)) * Cv;
      
      out_data->AvgSurfT += rad_temp*Cv;

      /** record mean surface temperature [C]  **/
      out_data->surf_temp += surf_temp * Cv;
      
      /**********************************
	Record Lake Snow Pack Variables 
      **********************************/
	  
      /** record snow water equivalence **/
      out_data->swq[0] += snow[veg][band].swq * Cv * 1000.;
      
      /** record snowpack depth **/
      out_data->snow_depth[0] += snow[veg][band].depth * Cv * 100.;
      
      /** record snow cover fraction **/
      out_data->coverage[0] += snow[veg][band].coverage * Cv;
      
      /** record snowpack cold content **/
      out_data->deltaCC[0] += energy[veg][band].deltaCC * Cv;
	  
      /** record snowpack advection **/
      out_data->advection[0] += energy[veg][band].advection * Cv;
	  
      /** record snow energy flux **/
      out_data->snow_flux[0] += energy[veg][band].snow_flux * Cv;
	  
      /** record refreeze energy **/
      out_data->refreeze_energy[0] += energy[veg][band].refreeze_energy * Cv;

      // Store Lake Specific Variables 
      out_data->lake_ice_temp   = lake_var.tempi;
      out_data->lake_ice_height = lake_var.hice;
      out_data->lake_ice_fract  = lake_var.fraci*Cv;
      out_data->lake_depth      = lake_var.ldepth;
      //     out_data->lake_surf_area  = Clake;
      out_data->lake_surf_area  = lake_var.sarea;
      out_data->lake_volume     = lake_var.volume;
      out_data->lake_surf_temp  = lake_var.temp[0];
      out_data->surfstor        = (lake_var.volume/lake_var.surface[0])*Cv*1000.;
      

#endif /* not OPTIMIZE */
	  
    }
  }
#endif /* LAKE_MODEL */
  
#if !OPTIMIZE

  /** record aerodynamic resistance **/
  if (out_data->aero_cond < HUGE_RESIST) {
    out_data->aero_resist = 1 / out_data->aero_cond;
  }

  /** record radiative temperature **/
  out_data->rad_temp = pow(out_data->rad_temp,0.25);

  /** record net radiation **/
  out_data->r_net    = out_data->net_short + out_data->net_long;

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data->prec;
  outflow = out_data->evap + out_data->runoff + out_data->baseflow;
  storage = 0.;
  for ( index = 0; index < options.Nlayer; index++ )
    if ( options.MOISTFRACT )
      storage += (out_data->moist[index] + out_data->ice[index]) 
	* depth[index] * 1000;
    else
      storage += out_data->moist[index] + out_data->ice[index];
  storage += out_data->swq[0] + out_data->snow_canopy[0] + out_data->Wdew;
#if LAKE_MODEL
  if ( options.LAKES ) 
    storage += ( out_data->lake_moist );
#endif // LAKE_MODEL
  calc_water_balance_error(rec,inflow,outflow,storage);
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data->net_short + out_data->net_long,
			      out_data->latent 
			      + out_data->latent_sub[0], 
			      out_data->sensible 
			      + out_data->advected_sensible[0],
			      out_data->grnd_flux + out_data->deltaH 
			      + out_data->fusion, 
			      out_data->advection[0] - out_data->deltaCC[0] 
			      + out_data->refreeze_energy[0] );

  /*************
    Write Data
  *************/

  if(rec >= skipyear) {
    if (options.ALMA_OUTPUT) {
      /* Convert to ALMA standard */
      out_data_alma = (out_data_alma_struct *) calloc(1,sizeof(out_data_alma_struct)); 
      conv_results_vic2alma(out_data, dt, depth, out_data_alma, rec);
      /* write ALMA variables */
      write_data_alma(out_data_alma, outfiles, dmy);
      // Free memory
      free((char *)out_data_alma); 
    }
    else {
      write_data(out_data, outfiles, dmy, dt);
    }
  }

#else /* not OPTIMIZE */

  if ( rec == 0 ) prtdt = 0;
  if ( prtdt == 0 ) {
    runoff = out_data->runoff;
    baseflow = out_data->baseflow;
    snow_depth = out_data->snow_depth[0];
    prtdt ++;
  }
  else {
    runoff += out_data->runoff;
    baseflow += out_data->baseflow;
    snow_depth += out_data->snow_depth[0];
    prtdt ++;
  }
  if ( prtdt == 24 / dt ) {
    out_data->runoff = runoff;
    out_data->baseflow = baseflow;
    out_data->snow_depth[0] = snow_depth / (double)(prtdt);
    prtdt = 0;
    if(rec >= skipyear)
      write_data(out_data, outfiles, dmy, dt);
  } 

#endif /* not OPTIMIZE */

  free((char *)out_data); 

}
