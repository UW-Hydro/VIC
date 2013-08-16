#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: put_data.c,v 4.2.2.9 2006/09/26 20:24:26 vicadmin Exp $";

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
	      veg_con_struct    *veg_con,
              out_data_file_struct   *out_data_files,
              out_data_struct   *out_data,
              save_data_struct  *save_data,
              double            *depth,
	      double            *dz,
	      double             dp,
	      double            *AreaFract,
	      char              *AboveTreeLine,
	      dmy_struct        *dmy,
              int                rec,
	      int                dt,
	      int                Nnodes,
	      int                skipyear,
	      double             elevation,
	      double            *Wcr,
	      double            *Wpwp,
	      double            *max_moist,
	      int month)
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
  03-12-03 modified to add additional energy balance variable storage 
           when output of snow bands is selected.                  KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  28-Sep-04 Replaced aero_resist[0] with aero_resist_used, the aerodynamic
	    resistance that was actually used in flux calculations.	TJB
  2005-11-21 (Port from 4.1.0) New aero-cond is aggregated instead of 
            aero-resist. GCT
  2006-Sep-11 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT options. TJB
  2006-Sep-14 Implemented ALMA-compliant input and output; uses the
	      new save_data structure; tracks more variables.  TJB
  2006-Sep-18 Implemented aggregation of output variables.  TJB
  2006-Sep-23 Changed MOL_WT_RATIO to EPS.  TJB

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
  double                  Cv;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_potevap_lake; //ingjerd dec2008
  double                  tmp_potevap; //ingjerd oct2009
  double                  tmp_moist;
  double                  tmp_root_moist; //ingjerd dec2008
  double                  tmp_ice;
  double                  cv_baresoil;
  double                  cv_veg;
  double                  cv_snow;
  double                  rad_temp;
  double                  surf_temp;
  double                  inflow;
  double                  outflow;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];
  double                  h; //ingjerd dec 2008
  double dummy;
  double dummy_moist;
  double dummy_depth1;
  double dummy_depth2;

  int                     v;
  int                     i;
  int                     dt_sec;
  int                     out_dt_sec;
  int                     out_step_ratio;
  int veg_class; //ingjerd jul 2009
  static int              step_count;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

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
  out_data[OUT_ORIGPREC].data[0]   = atmos->orig_prec[NR];
  //if( atmos->orig_prec[NR]>0) printf("put_data outorigprec %f orig_prec %f\n",  out_data[OUT_ORIGPREC].data[0],atmos->orig_prec[NR]);
  out_data[OUT_PRESSURE].data[0]  = atmos->pressure[NR];
  out_data[OUT_QAIR].data[0]      = EPS * atmos->vp[NR]/atmos->pressure[NR];
  out_data[OUT_RAINF].data[0]     = atmos->out_rain;
  out_data[OUT_REL_HUMID].data[0] = 100.*atmos->vp[NR]/(atmos->vp[NR]+atmos->vpd[NR]);
  out_data[OUT_SHORTWAVE].data[0] = atmos->shortwave[NR];
  out_data[OUT_SNOWF].data[0]     = atmos->out_snow;
  out_data[OUT_VP].data[0]        = atmos->vp[NR];
  out_data[OUT_WIND].data[0]      = atmos->wind[NR];
  out_data[OUT_VPD].data[0]       = atmos->vpd[NR]; //ingjerd feb 2009
  out_data[OUT_SLOPE].data[0]     = svp_slope(atmos->air_temp[NR])*1000.; //ingjerd feb 2009
  h = 287/9.81 * ((atmos->air_temp[NR] + 273.15) + 0.5 * (double)elevation * LAPSE_PM);
  out_data[OUT_GAMMA].data[0]     =  1628.6 * (PS_PM * exp(-(double)elevation/h))/
      (2501000 - 2361 * atmos->air_temp[NR]); //ingjerd feb 2009
 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
  /****************************************
    Store Output for all Vegetation Types
  ****************************************/
  for ( veg = 0 ; veg <= veg_con[0].vegetat_type_num ; veg++) {
    
    if ( veg < veg_con[0].vegetat_type_num ) 
      Cv = veg_con[veg].Cv;
    else
      Cv = (1.0 - veg_con[0].Cv_sum);

    /*ingjerd. if irrigated vegetation exist, adjust fractions according to month.
      do only if irrigation is taken into account. no irrigation: two classes with half the area. */
    if(veg_con->irrveg==1 && options.IRRIGATION) { 
      veg_class=veg_con[veg].veg_class;
      if(veg==veg_con[0].vegetat_type_num-1) Cv=Cv*2*veg_lib[veg_class].irrpercent[month]; //irrigated fraction this month
      if(veg==veg_con[0].vegetat_type_num-2) Cv=Cv*2*(1-veg_lib[veg_class].irrpercent[month]); //nonirrigated fraction this month
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
	    if ( veg <= veg_con[0].vegetat_type_num ) {
	      tmp_evap += snow[veg][band].canopy_vapor_flux * 1000.;
	      out_data[OUT_SUB_CANOP].data[0] += snow[veg][band].canopy_vapor_flux 
		* 1000. * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    }
	    if ( veg < veg_con[0].vegetat_type_num ) {
	      tmp_evap += veg_var[dist][veg][band].canopyevap;
	      out_data[OUT_EVAP_CANOP].data[0] += veg_var[dist][veg][band].canopyevap 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	      //** ingjerd record potevap, jun 2009 **/
	      out_data[OUT_POTEVAP_LAKE].data[0] +=  veg_var[dist][veg][band].potevap_lake 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      //  printf("put_data NR=%d veg=%d potevap_lake=%f %f %f %f\n",NR,veg,out_data[OUT_POTEVAP_LAKE].data[0],veg_var[dist][veg][band].potevap,Cv,mu);
	      out_data[OUT_POTEVAP].data[0] +=  veg_var[dist][veg][band].potevap 
		* Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	      //printf("put_data NR=%d veg=%d %f %f %f %f\n",NR,veg,out_data[OUT_POTEVAP].data[0],veg_var[dist][veg][band].potevap,Cv,mu);
	    }
	    out_data[OUT_EVAP].data[0] += tmp_evap * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 

	    /** record runoff **/
	    out_data[OUT_RUNOFF].data[0]   += cell[dist][veg][band].runoff 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
	    
	    /** record baseflow **/
	    out_data[OUT_BASEFLOW].data[0] += cell[dist][veg][band].baseflow 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 

	    /** record water extracted locally. ingjerd feb 2009, 
                thereafter subtract this from runoff/baseflow. 
		Works only for one snowband and distprec=FALSE! 
                **/
	    out_data[OUT_EXTRACT_WATER].data[0] += cell[dist][veg][band].extract_water 
	      * Cv * mu * AreaFract[band] * TreeAdjustFactor[band]; 
	    //if(rec==155) printf("put_data A rec=%d veg=%d Cv=%f prec=%f origprec=%f extract_water=%f runoff=%f base=%f\n",
	    //     rec,veg,Cv,out_data[OUT_PREC].data[0],out_data[OUT_ORIGPREC].data[0],out_data[OUT_EXTRACT_WATER].data[0],out_data[OUT_RUNOFF].data[0],out_data[OUT_BASEFLOW].data[0]);
	    if (veg == veg_con[0].vegetat_type_num-1 && out_data[OUT_EXTRACT_WATER].data[0]>0 ) { 
	      if(out_data[OUT_RUNOFF].data[0]>=out_data[OUT_EXTRACT_WATER].data[0]) {
		out_data[OUT_RUNOFF].data[0]-=out_data[OUT_EXTRACT_WATER].data[0];
	      }
	      else {
		dummy=out_data[OUT_EXTRACT_WATER].data[0]-out_data[OUT_RUNOFF].data[0];
		out_data[OUT_RUNOFF].data[0]=0.;
		if(out_data[OUT_BASEFLOW].data[0]<dummy) 
		  out_data[OUT_BASEFLOW].data[0]=0.;
		else
		  out_data[OUT_BASEFLOW].data[0]-=dummy;
	      }
	      if(out_data[OUT_BASEFLOW].data[0]<0) printf("put_data extractwater uffda! %f\n",out_data[OUT_BASEFLOW].data[0]);
	    }
	    //if(rec==155) printf("put_data B rec%d prec=%f extract=%f run=%f base=%f\n",
	    //     rec,out_data[OUT_PREC].data[0],out_data[OUT_EXTRACT_WATER].data[0],out_data[OUT_RUNOFF].data[0],out_data[OUT_BASEFLOW].data[0]);		    
	  

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
	      tmp_ice   = cell[dist][veg][band].layer[index].ice;
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

            /* ingjerd added rootstress statement */
	    out_data[OUT_ROOT_STRESS].data[0] += cell[dist][veg][band].rootstress
              * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];

	    /** record layer temperatures **/
	    for(index=0;index<options.Nlayer;index++) {
	      out_data[OUT_SOIL_TEMP].data[index] += cell[dist][veg][band].layer[index].T
                * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
            }

	  }
	}
      }

      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0. && ( veg == veg_con[0].vegetat_type_num || ( !AboveTreeLine[band] || (AboveTreeLine[band] && !veg_lib[veg_con[veg].veg_class].overstory)))) {

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
	  
	  /**********************************
            Record Energy Balance Variables
	  **********************************/

	  /** record surface radiative temperature **/
	  if(snow[veg][band].swq>0) 
	    rad_temp = snow[veg][band].surf_temp + KELVIN;
	  else
	    rad_temp = energy[veg][band].T[0] + KELVIN;
	  
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

	  /** record soil surface temperature **/
	  surf_temp = energy[veg][band].T[0];
	  // MODIFIED FOR USE WITH ROSEMOUNT SIMULATIONS
	  //surf_temp = (energy[veg][band].T[0] + energy[veg][band].T[1])/2.;
	  
	  /** record net shortwave radiation **/
	  out_data[OUT_NET_SHORT].data[0] += energy[veg][band].shortwave
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record net longwave radiation **/
	  out_data[OUT_NET_LONG].data[0]  += energy[veg][band].longwave
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record incoming longwave radiation **/
	  out_data[OUT_LONGWAVE].data[0]  += ((energy[veg][band].longwave + STEFAN_B 
				  * (rad_temp) * (rad_temp)
				  * (rad_temp) * (rad_temp))
				 * Cv * AreaFract[band] * TreeAdjustFactor[band]);
	  //printf("put_data rad_temp=%f net=%f lw=%f\n",rad_temp,energy[veg][band].longwave,out_data[OUT_LONGWAVE].data[0]);

	  /** record albedo **/
	  out_data[OUT_ALBEDO].data[0]    += energy[veg][band].albedo
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record latent heat flux **/
	  out_data[OUT_LATENT].data[0]    -= energy[veg][band].latent
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record sensible heat flux **/
	  out_data[OUT_SENSIBLE].data[0]  -= energy[veg][band].sensible
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record ground heat flux (+ heat storage) **/
	  out_data[OUT_GRND_FLUX].data[0] -= (energy[veg][band].grnd_flux
				  + energy[veg][band].deltaH)
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];

	  /** record heat storage **/
	  out_data[OUT_DELTAH].data[0]    -= energy[veg][band].deltaH
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record energy balance error **/
	  out_data[OUT_ENERGY_ERROR].data[0] += energy[veg][band].error
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record radiative effective temperature [K], 
	      emissivities set = 1.0  **/
	  out_data[OUT_RAD_TEMP].data[0] += ((rad_temp) * (rad_temp) 
				 * (rad_temp) * (rad_temp)) 
	    * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record mean surface temperature [C]  **/
	  out_data[OUT_SURF_TEMP].data[0] += surf_temp * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
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
	  
	  /** record snow cover fraction **/
	  out_data[OUT_SNOW_COVER].data[0]
	    += snow[veg][band].coverage * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack cold content **/
	  out_data[OUT_DELTACC].data[0]
	    += energy[veg][band].deltaCC * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snowpack advection **/
	  out_data[OUT_ADVECTION].data[0]
	    += energy[veg][band].advection * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record snow energy flux **/
	  out_data[OUT_SNOW_FLUX].data[0]
	    += energy[veg][band].snow_flux * Cv * AreaFract[band] * TreeAdjustFactor[band];
	  
	  /** record refreeze energy **/
	  out_data[OUT_REFREEZE_ENERGY].data[0]
	    += energy[veg][band].refreeze_energy 
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
	    out_data[OUT_REFREEZE_ENERGY_BAND].data[band]
	      += energy[veg][band].refreeze_energy * Cv;
	    
	    /** record band net downwards shortwave radiation **/
	    out_data[OUT_NET_SHORT_BAND].data[band]
	      += energy[veg][band].shortwave * Cv;

	    /** record band net downwards longwave radiation **/
	    out_data[OUT_NET_LONG_BAND].data[band]
	      += energy[veg][band].longwave * Cv;

	    /** record band albedo **/
	    out_data[OUT_ALBEDO_BAND].data[band]
	      += energy[veg][band].albedo * Cv;

	    /** record band net latent heat flux **/
	    out_data[OUT_LATENT_BAND].data[band]
	      -= energy[veg][band].latent * Cv;

	    /** record band net sensible heat flux **/
	    out_data[OUT_SENSIBLE_BAND].data[band]
	      -= energy[veg][band].sensible * Cv;

	    /** record band net ground heat flux **/
	    out_data[OUT_GRND_FLUX_BAND].data[band]
	      -= (energy[veg][band].grnd_flux
	         + energy[veg][band].deltaH) * Cv;

	  }
	}
      }

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
  out_data[OUT_SNOW_MELT].data[0] = out_data[OUT_SNOWF].data[0]-out_data[OUT_DELSWE].data[0]-out_data[OUT_SUB_SNOW].data[0]-out_data[OUT_SUB_CANOP].data[0];
  if (out_data[OUT_SNOW_MELT].data[0] < 0) out_data[OUT_SNOW_MELT].data[0] = 0;

  //ingjerd added 
  out_data[OUT_ESS_SNOW].data[0] = out_data[OUT_SUB_SNOW].data[0] + out_data[OUT_SUB_CANOP].data[0];  

  out_data[OUT_SOIL_MALL].data[0] = 0;
  for (index=0; index<options.Nlayer; index++)
  out_data[OUT_SOIL_MALL].data[0] += out_data[OUT_SOIL_MOIST].data[index];

  // Energy terms
  out_data[OUT_REFREEZE].data[0] = (out_data[OUT_REFREEZE_ENERGY].data[0]/Lf)*dt_sec;
  out_data[OUT_MELT_ENERGY].data[0] = out_data[OUT_SNOW_MELT].data[0]*Lf/dt_sec;
  out_data[OUT_R_NET].data[0] = out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0];

  // Save current moisture state for use in next time step
  save_data->total_soil_moist = 0;
  for (index=0; index<options.Nlayer; index++) {
    save_data->total_soil_moist += out_data[OUT_SOIL_MOIST].data[index];
  }
  save_data->swe = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0];
  save_data->wdew = out_data[OUT_WDEW].data[0];

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data[OUT_PREC].data[0]; 
  outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0] + out_data[OUT_EXTRACT_WATER].data[0]; //ingjerd added extract_water feb 2009
  //outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0];
  //if(rec==154 || rec==155 || rec==156) printf("put_data C rec=%d prec=%f extract=%f evap=%f run=%f base=%f in-out: %f snow:%f\n",rec,out_data[OUT_PREC].data[0],out_data[OUT_EXTRACT_WATER].data[0],out_data[OUT_EVAP].data[0],out_data[OUT_RUNOFF].data[0],out_data[OUT_BASEFLOW].data[0],out_data[OUT_PREC].data[0]-out_data[OUT_EXTRACT_WATER].data[0]-out_data[OUT_EVAP].data[0]-out_data[OUT_RUNOFF].data[0]-out_data[OUT_BASEFLOW].data[0],out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0]);
  storage = 0.;
  for(index=0;index<options.Nlayer;index++)
    if(options.MOISTFRACT)
      storage += (out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index]) 
	* depth[index] * 1000;
    else
      storage += out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index];
  storage += out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] + out_data[OUT_WDEW].data[0];

  //if(rec<35) 
  //  printf("put_data D rec=%d storage=%f\n",rec,storage);

  calc_water_balance_error(rec,inflow,outflow,storage);

  /********************
    Check Energy Balance 
    ********************/
  if(options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0],
			      out_data[OUT_LATENT].data[0], out_data[OUT_SENSIBLE].data[0],
			      out_data[OUT_GRND_FLUX].data[0], out_data[OUT_ADVECTION].data[0]
			      - out_data[OUT_DELTACC].data[0] - out_data[OUT_SNOW_FLUX].data[0]
			      + out_data[OUT_REFREEZE_ENERGY].data[0]);

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
      out_data[OUT_TRANSP_VEG].aggdata[0] /= out_dt_sec;
      out_data[OUT_INFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_PREC].aggdata[0] /= out_dt_sec;
      out_data[OUT_RAINF].aggdata[0] /= out_dt_sec;
      out_data[OUT_REFREEZE].aggdata[0] /= out_dt_sec;
      out_data[OUT_RUNOFF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOW_MELT].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOWF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] += out_data[OUT_SUB_CANOP].aggdata[0];
      out_data[OUT_BARESOILT].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_PACK_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_SURF_TEMP].aggdata[0] += KELVIN;
      for (index=0; index<options.Nlayer; index++) {
        out_data[OUT_SOIL_TEMP].aggdata[index] += KELVIN;
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
