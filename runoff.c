#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";
int  runoff(layer_data_struct *layer_wet,
	    layer_data_struct *layer_dry,
            energy_bal_struct *energy,
            soil_con_struct   *soil_con,
            double            *runoff_wet, 
	    double            *runoff_dry, 
	    double            *baseflow_wet,
	    double            *baseflow_dry,
	    double            *ppt, 
#if EXCESS_ICE
	    int                SubsidenceUpdate,
#endif // EXCESS_ICE
#if SPATIAL_FROST
	    double            *frost_fract,
#endif // SPATIAL_FROST
	    double             mu,
	    int                dt,
            int                Nnodes,
	    int                band,
	    int                rec,
	    int                iveg)
/**********************************************************************
	runoff.c	Keith Cherkauer		May 18, 1996

  This subroutine calculates infiltration and runoff from the surface,
  gravity driven drainage between all soil layers, and generates 
  baseflow from the bottom layer..
  
  sublayer indecies are always [layer number][sublayer number]
  [layer number] is the current VIC model moisture layer
  [sublayer number] is the current sublayer number where: 
         0 = thawed sublayer, 1 = frozen sublayer, and 2 = unfrozen sublayer.
	 when the model is run withoputfrozen soils, the sublayer number
	 is always = 2 (unfrozen).

  UNITS:	Ksat (mm/day)
		Q12  (mm/time step)
		moist (mm)
		inflow (mm)
                runoff (mm)

  Variables:
	ppt	incoming precipitation and snow melt
	mu	fraction of area that receives precipitation
	inflow	incoming water corrected for fractional area of precip (mu)

  MODIFICATIONS:
  5/22/96 Routine modified to account for spatially varying
	  precipitation, and it's effects on runoff.	KAC
  11/96	  Code modified to account for extra model layers
  	  needed for frozen soils modeling.		KAC
  1/9/97  Infiltration and other rate parameters modified
	  for time scales of less than 1 day.		KAC
  4-1-98  Soil moisture transport is now done on an hourly time
          step, irregardless to the model time step, to prevent
          numerical stabilities in the solution	Dag and KAC
  01-24-00 simplified handling of soil moisture for the
           frozen soil algorithm.  all option selection
	   now use the same soil moisture transport method   KAC
  6-8-2000 modified to handle spatially distributed soil frost  KAC
  06-07-03 modified so that infiltration is computed using only the
           top two soil moisture layers, rather than all but the
           bottom most layer.  This preserves the functionality
           of the original model design, but is more realistic for
           handling multiple soil moisture layers
  06-Sep-03   Changed calculation of dt_baseflow to go to zero when
              soil liquid moisture <= residual moisture.  Changed
              block that handles case of total soil moisture < residual
              moisture to not allow dt_baseflow to go negative.		TJB
  17-May-04   Changed block that handles baseflow when soil moisture
	      drops below residual moisture.  Now, the block is only
	      entered if baseflow > 0 and soil moisture < residual,
	      and the amount of water taken out of baseflow and given
	      to the soil cannot exceed baseflow.  In addition, error
	      messages are no longer printed, since it isn't an error
	      to be in that block.					TJB
  2007-Apr-04 Modified to return Error status from 
              distribute_node_moisture_properties			GCT/KAC
  2007-Apr-24 Passes soil_con->Zsum_node to distribute_node_moisture_properties.  JCA
  2007-Jun-13 Fixed bug arising from earlier fix to dt_baseflow
	      calculation.  Earlier fix took residual moisture
	      into account in the linear part of the baseflow eqn,
	      but not in the non-linear part.  Now we take residual
	      moisture into account correctly throughout the whole
	      equation.  Also re-wrote equation in simpler form.	TJB
  2007-Aug-15 Changed SPATIAL_FROST if statement to enclose the correct
              end-bracket for the frost_area loop.			JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
              Including adding SubsidenceUpdate flag for parts
              of the routine that will be used if redistributing
              soil moisture after subsidence.
  2007-Sep-18 Modified to correctly handle evaporation from spatially
	      distributed soil frost.  Original version could produce
	      negative soil moisture in fractions with high ice content
	      since only total evaporation was checked versus total
	      liquid water content, not versus available liquid water
	      in each frost subsection.					KAC via TJB
  2007-Sep-20 Removed logic that reset resid_moist[i].  Previously,
	      resid_moist[i] was reset to 0 for i > 0 when
	      resid_moist[0] == 0.  Such resetting of soil properties
	      was deemed unnecessary and confusing, since VIC would end
	      up using different residual moisture values than those
	      specified by the user.  If a user truly wants to specify
	      residual moisture in all layers to be 0, the user should
	      set these explicitly in the soil parameter file.  Also
	      fixed typo in fprintf() on line 289.			TJB
**********************************************************************/
{  
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct  debug;
#endif // LINK_DEBUG

  char               ErrStr[MAXSTRING];
  int                firstlayer, lindex, sub;
  int                i;
  int                last_layer[MAX_LAYERS*3];
  int                last_sub[MAX_LAYERS*3];
  int                froz_solid[MAX_LAYERS*3];
  int                last_index;
  int                last_cnt;
  int                tmp_index;
  int                time_step;
  int                Ndist;
  int                dist;
  int                storesub;
  int                tmpsub;
  int                tmplayer;
  int                frost_area;
  int                ErrorFlag;
  double             ex, A, i_0, basis, frac;
  double             inflow;
  double             last_moist;
  double             resid_moist[MAX_LAYERS];
  double             org_moist[MAX_LAYERS];
  double             moist[MAX_LAYERS];
  double             ice[MAX_LAYERS];
  double             max_moist[MAX_LAYERS];
  double             max_infil;
  double             Ksat[MAX_LAYERS];
  double             Q12[MAX_LAYERS-1];
  double            *kappa;
  double            *Cs;
  double            *M;
  double             Dsmax;
  double             top_moist;
  double             top_max_moist;
  double             tmp_inflow;
  double             tmp_moist;
  double             dt_inflow, dt_outflow;
  double             dt_runoff;
  double            *runoff_ptr;
  double            *baseflow_ptr;
  double             runoff[FROST_SUBAREAS];
  double             baseflow[FROST_SUBAREAS];
  double             actual_frost_fract[FROST_SUBAREAS];
  double             tmp_mu;
  double             dt_baseflow;
  double             rel_moist;
  double             evap[MAX_LAYERS][FROST_SUBAREAS];
  double             sum_moist;
  double             evap_percent;
  double             evap_sum;
#if LOW_RES_MOIST
  double             b[MAX_LAYERS];
  double             matric[MAX_LAYERS];
  double             avg_matric;
  double             spatial_fract;
#endif // LOW_RES_MOIST
#if EXCESS_ICE
  double             excess_water;
  double             net_excess_water;
  double             moist_prior;
  double             total_evap;
#endif //EXCESS_ICE
  layer_data_struct *layer;
  layer_data_struct  tmp_layer;

  /** Set Residual Moisture **/
  for ( i = 0; i < options.Nlayer; i++ ) 
    resid_moist[i] = soil_con->resid_moist[i] * soil_con->depth[i] * 1000.;

  /** Initialize Other Parameters **/
  if ( options.DIST_PRCP ) Ndist = 2;
  else Ndist = 1;
  tmp_mu = mu;

  /** Allocate and Set Values for Soil Sublayers **/
  for ( dist = 0; dist < Ndist; dist++ ) {
    
    /* Loop through wet and dry cell fractions */

    if(dist>0) {
      layer        = layer_dry;
      runoff_ptr   = runoff_dry;
      baseflow_ptr = baseflow_dry;
      mu           = (1. - mu);
    }
    else {
      layer        = layer_wet;
      runoff_ptr   = runoff_wet;
      baseflow_ptr = baseflow_wet;
    }

    *runoff_ptr   = 0;
    *baseflow_ptr = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
      baseflow[frost_area] = 0;
      
    if(mu>0.) {
	
#if SPATIAL_FROST
      for ( lindex = 0; lindex < options.Nlayer; lindex++ ) {
	evap[lindex][0] = layer[lindex].evap/(double)dt;
	org_moist[lindex] = layer[lindex].moist;
	layer[lindex].moist = 0;
        if ( evap[lindex][0] != 0 ) {
          // if there is evaporation
          sum_moist = 0;
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
            // compute total available soil moisture between frost sub areas.
            sum_moist += (org_moist[lindex] - layer[lindex].ice[frost_area]) * frost_fract[frost_area];
          // compute fraction of avaiable soil moisture that is evaporated
          evap_percent = evap[lindex][0] / sum_moist;
          // distribute evaporation between frost sub areas by percentage
          evap_sum = evap[lindex][0];
          for ( frost_area = FROST_SUBAREAS - 1; frost_area >= 0; frost_area-- ) {
            evap[lindex][frost_area] = ( org_moist[lindex]
                                         - layer[lindex].ice[frost_area] )
              * evap_percent;
            evap_sum -= evap[lindex][frost_area] * frost_fract[frost_area];
          }
          if ( evap_sum > SMALL || evap_sum < -SMALL ) {
            fprintf(stderr,"Evap_sum = %f\n", evap_sum);
          }
        }
        else {
          for ( frost_area = FROST_SUBAREAS - 1; frost_area > 0; frost_area-- )
            evap[lindex][frost_area] = evap[lindex][0];
        }
      }
      
      for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
	
#else
      // store current evaporation
      for ( lindex = 0; lindex < options.Nlayer; lindex++ )
        evap[lindex][0] = layer[lindex].evap/(double)dt;

      frost_area = 0;

#endif // SPATIAL_FROST
      
	/** ppt = amount of liquid water coming to the surface **/
	inflow = ppt[dist];
	
	/**************************************************
	  Initialize Variables
	**************************************************/
	for ( lindex = 0; lindex < options.Nlayer; lindex++ ) {
	  Ksat[lindex]         = soil_con->Ksat[lindex] / 24.;
#if LOW_RES_MOIST
	  b[lindex]            = (soil_con->expt[lindex] - 3.) / 2.;
#endif // LOW_RES_MOIST
	  
	  /** Set Layer Unfrozen Moisture Content **/
#if SPATIAL_FROST
	  moist[lindex] = org_moist[lindex] - layer[lindex].ice[frost_area];
#else
	  moist[lindex] = layer[lindex].moist - layer[lindex].ice;
#endif // SPATIAL_FROST
	  if(moist[lindex]<0) {
	    if(fabs(moist[lindex]) < 1e-6)
	      moist[lindex] = 0;
#if SPATIAL_FROST
	    else if (layer[lindex].moist < layer[lindex].ice[frost_area])
#else
	    else if (layer[lindex].moist < layer[lindex].ice)
#endif // SPATIAL_FROST
	      moist[lindex] = 0;
	    else {
	      fprintf(stderr, "ERROR in runoff(): Layer %d has negative soil moisture, %f\n",
                      lindex, moist[lindex]);
	      return(ERROR);
	    }
	  }
	  
	  /** Set Layer Ice Content **/
#if SPATIAL_FROST
	  ice[lindex]       = layer[lindex].ice[frost_area];
#else
	  ice[lindex]       = layer[lindex].ice;
#endif // SPATIAL_FROST
	  
	  /** Set Layer Maximum Moisture Content **/
	  max_moist[lindex] = soil_con->max_moist[lindex];

	}//each layer
	
	/******************************************************
          In case of subsidence, check if total soil column 
          moisture exceeds maximum capacity, and run simple
          scenario if true.
	******************************************************/
#if EXCESS_ICE
	if(SubsidenceUpdate == 1){ 
	  excess_water = 0;
	  net_excess_water = 0;
	  for ( lindex = 0; lindex < options.Nlayer; lindex++ ) {
	    net_excess_water += (moist[lindex]+ice[lindex] - max_moist[lindex]);
	    if( (moist[lindex]+ice[lindex]) > max_moist[lindex])
	      excess_water += (moist[lindex]+ice[lindex] - max_moist[lindex]);	
	  }	
	}

	if(SubsidenceUpdate == 1 && net_excess_water >= 0){//run simple scenario
	  /* set all layers to saturation*/
	  for ( lindex = 0; lindex < options.Nlayer; lindex++ ){
	    moist[lindex] = max_moist[lindex] - ice[lindex];
	    if(moist[lindex]<0){
	      fprintf(stderr, "ERROR in runoff(): Layer %d has negative soil moisture, %f\n",
                      lindex, moist[lindex]);
	      return(ERROR);
	    }
	  }
	  
	  /*estimate baseflow contribution, same method as below*/
	  lindex = options.Nlayer-1;
	  Dsmax = soil_con->Dsmax / 24.;  
	  for (time_step = 0; time_step < dt; time_step++) {
	    /** Compute relative moisture **/
	    rel_moist = (moist[lindex]-resid_moist[lindex])
	      / (soil_con->max_moist[lindex]-resid_moist[lindex]);
	    /** Compute baseflow as function of relative moisture **/
	    frac = Dsmax * soil_con->Ds / soil_con->Ws;
	    dt_baseflow = frac * rel_moist;
	    if (rel_moist > soil_con->Ws) {
	      frac = (rel_moist - soil_con->Ws) / (1 - soil_con->Ws);
	      dt_baseflow += Dsmax * (1 - soil_con->Ds / soil_con->Ws)
		* pow(frac,soil_con->c);
	    }	    
	    if(dt_baseflow < 0) 
	      dt_baseflow = 0;
	    baseflow[frost_area] += dt_baseflow;
	  }
	  
	  /*calculate total evap*/
	  total_evap = 0;
	  for ( lindex = 0; lindex < options.Nlayer; lindex++ ) 
	    total_evap += evap[lindex][frost_area]*(double)dt;

	  /* estimate runoff as sum of excess water */
	  runoff[frost_area] = net_excess_water + inflow - baseflow[frost_area] - total_evap;
	  if(runoff[frost_area] < 0) {
	    baseflow[frost_area] += runoff[frost_area];  
	    runoff[frost_area] = 0;
	  }

	}//end simple scenario

	else {
	  /******************************************************
           For now, do a crude redistribution of soil moisture, so
           that moist does not exceed max_moist for any layer. 
           Then continue with usual runoff routine.
           Eventually, may want to make this more sophisticated.
           (Note: This case is rare compared to case above.)
	  ******************************************************/
	  if(SubsidenceUpdate == 1 && excess_water > 0){
	    //fill from bottom up with excess water only
	    for(lindex=(options.Nlayer-1);lindex>=0;lindex--) {
	      if(max_moist[lindex] > (moist[lindex] + ice[lindex])) {//if not a subsidence layer
		if((max_moist[lindex] - (moist[lindex] + ice[lindex])) <= excess_water){//can't take all excess
		  if(excess_water > 0){
		    moist_prior = moist[lindex];
		    moist[lindex] = max_moist[lindex] - ice[lindex];//set to saturation
		    excess_water -= (moist[lindex]-moist_prior);
		  }
		}
		else {//can take all excess
		  if(excess_water > 0){
		    moist_prior = moist[lindex];
		    moist[lindex] += excess_water;//take-up all excess
		    excess_water -= (moist[lindex]-moist_prior);
		  }
		}
	      }
	      else //if a subsidence layer
		moist[lindex] = max_moist[lindex] - ice[lindex];//set to saturation  
	    }
	    
	  }

#endif //EXCESS_ICE	
	  
	  /******************************************************
          Runoff Based on Soil Moisture Level of Upper Layers
	  ******************************************************/
	  
	  top_moist = 0.;
	  top_max_moist=0.;
	  for(lindex=0;lindex<2;lindex++) {
	    top_moist += (moist[lindex] + ice[lindex]);
	    top_max_moist += max_moist[lindex];
	  }
	  if(top_moist>top_max_moist) top_moist = top_max_moist;
	  
	  /**************************************************
          Calculate Runoff from Surface
	  **************************************************/
	  
	  /** Runoff Calculations for Top Layer Only **/
	  /** A and i_0 as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
	  
	  max_infil = (1.0+soil_con->b_infilt) * top_max_moist;
	  
	  ex        = soil_con->b_infilt / (1.0 + soil_con->b_infilt);
	  A         = 1.0 - pow((1.0 - top_moist / top_max_moist),ex);
	  i_0       = max_infil * (1.0 - pow((1.0 - A),(1.0 / soil_con->b_infilt))); 
	  /* Maximum Inflow */
	  
	  /** equation (3a) Wood et al. **/
	  
	  if (inflow == 0.0) runoff[frost_area] = 0.0;
	  else if (max_infil == 0.0) runoff[frost_area] = inflow;
	  else if ((i_0 + inflow) > max_infil) 
	    runoff[frost_area] = inflow - top_max_moist + top_moist;
	  
	  /** equation (3b) Wood et al. (wrong in paper) **/
	  else {
	    basis = 1.0 - (i_0 + inflow) / max_infil;
	    runoff[frost_area] = (inflow - top_max_moist + top_moist 
				  + top_max_moist
				  * pow(basis,1.0*(1.0+soil_con->b_infilt)));
	  }
	  if ( runoff[frost_area] < 0. ) runoff[frost_area] = 0.;
	  
	  /**************************************************
	  Compute Flow Between Soil Layers (using an hourly time step)
	  **************************************************/
	  
	  dt_inflow  =  inflow / (double) dt;
	  dt_outflow =  0.0;
	  
	  for (time_step = 0; time_step < dt; time_step++) {
	    inflow   = dt_inflow;
	    last_cnt = 0;
	    
#if LOW_RES_MOIST
	    for( lindex = 0; lindex < options.Nlayer; lindex++ ) {
	      if( (tmp_moist = moist[lindex] - evap[lindex][frost_area]) 
		  < resid_moist[lindex] )
		tmp_moist = resid_moist[lindex];
	      if(tmp_moist > resid_moist[lindex])
		matric[lindex] = soil_con->bubble[lindex] 
		  * pow( (tmp_moist - resid_moist[lindex]) 
			 / (soil_con->max_moist[lindex] - resid_moist[lindex]), 
			 -b[lindex]);
	      else
		matric[lindex] = HUGE_RESIST;
	    }
#endif // LOW_RES_MOIST
	    
	    /*************************************
            Compute Drainage between Sublayers 
	    *************************************/
	    
	    for( lindex = 0; lindex < options.Nlayer-1; lindex++ ) {
	      
	      /** Brooks & Corey relation for hydraulic conductivity **/
	      
	      if((tmp_moist = moist[lindex] - evap[lindex][frost_area]) 
		 < resid_moist[lindex])
		tmp_moist = resid_moist[lindex];
	      
	      if(moist[lindex] > resid_moist[lindex]) {
#if LOW_RES_MOIST
		avg_matric = pow( 10, (soil_con->depth[lindex+1] 
				       * log10(fabs(matric[lindex]))
				       + soil_con->depth[lindex]
				       * log10(fabs(matric[lindex+1])))
				  / (soil_con->depth[lindex] 
				     + soil_con->depth[lindex+1]) );
		tmp_moist = resid_moist[lindex]
		  + ( soil_con->max_moist[lindex] - resid_moist[lindex] )
		  * pow( ( avg_matric / soil_con->bubble[lindex] ), -1/b[lindex] );
#endif // LOW_RES_MOIST
		Q12[lindex] 
		  = Ksat[lindex] * pow(((tmp_moist - resid_moist[lindex])
					/ (soil_con->max_moist[lindex]
					   - resid_moist[lindex])),
				       soil_con->expt[lindex]); 
	      }
	      else Q12[lindex] = 0.;
	      last_layer[last_cnt] = lindex;
	    }
	    
	    /**************************************************
            Solve for Current Soil Layer Moisture, and
            Check Versus Maximum and Minimum Moisture
            Contents.  
	    **************************************************/
	    
	    firstlayer = TRUE;
	    last_index = 0;
	    for ( lindex = 0; lindex < options.Nlayer - 1; lindex++ ) {
	      
	      if ( lindex == 0 ) dt_runoff = runoff[frost_area] / (double) dt;
	      else dt_runoff = 0;
	      
	      /* Store moistures for water balance debugging */
#if LINK_DEBUG
	      if ( debug.PRT_BALANCE ) {
		if ( time_step == 0 ) {
		  if ( firstlayer )
#if SPATIAL_FROST
		    debug.inflow[dist][band][lindex+2] 
		      += (inflow - dt_runoff) * frost_fract[frost_area];
#else
		  debug.inflow[dist][band][lindex+2] = inflow - dt_runoff;
#endif // SPATIAL_FROST
		  else {
#if SPATIAL_FROST
		    debug.inflow[dist][band][lindex+2] 
		      += inflow * frost_fract[frost_area];
		    debug.outflow[dist][band][lindex+1] 
		      += inflow * frost_fract[frost_area];
#else
		    debug.inflow[dist][band][lindex+2] = inflow;
		    debug.outflow[dist][band][lindex+1] = inflow;
#endif // SPATIAL_FROST
		  }
		}
		else {
		  if ( firstlayer )
		    debug.inflow[dist][band][lindex+2]  += inflow - dt_runoff;
		  else {
		    debug.inflow[dist][band][lindex+2]  += inflow;
		    debug.outflow[dist][band][lindex+1] += inflow;
		  }
		}
	      }
#endif // LINK_DEBUG
	      
	      /* transport moisture for all sublayers **/
#if LINK_DEBUG
	      if(debug.DEBUG || debug.PRT_BALANCE) 
		last_moist = moist[lindex];
#endif // LINK_DEBUG
	      
	      tmp_inflow = 0.;
	      
	      /** Update soil layer moisture content **/
	      moist[lindex] = moist[lindex] + (inflow - dt_runoff) 
		- (Q12[lindex] + evap[lindex][frost_area]);
	      
	      /** Verify that soil layer moisture is less than maximum **/
	      if((moist[lindex]+ice[lindex]) > max_moist[lindex]) {
		tmp_inflow = (moist[lindex]+ice[lindex]) - max_moist[lindex];
		moist[lindex] = max_moist[lindex] - ice[lindex];
		
		if(lindex==0) {
		  Q12[lindex] += tmp_inflow;
		  tmp_inflow = 0;
		}
		else {
		  tmplayer = lindex;
		  while(tmp_inflow > 0) {
		    tmplayer--;
#if LINK_DEBUG
		    if(debug.PRT_BALANCE) {
		      /** Update debugging storage terms **/
		      debug.inflow[dist][band][lindex+2]  -= tmp_inflow;
		      debug.outflow[dist][band][lindex+1] -= tmp_inflow;
		    }
#endif // LINK_DEBUG
		    if ( tmplayer < 0 ) {
		      /** If top layer saturated, add to runoff **/
		      runoff[frost_area] += tmp_inflow;
		      tmp_inflow = 0;
		    }
		    else {
		      /** else add excess soil moisture to next higher layer **/
		      moist[tmplayer] += tmp_inflow;
		      if((moist[tmplayer]+ice[tmplayer]) > max_moist[tmplayer]) {
			tmp_inflow = ((moist[tmplayer] + ice[tmplayer])
				      - max_moist[tmplayer]);
			moist[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
		      }
		      else tmp_inflow=0;
		    }
		  }
		} /** end trapped excess moisture **/
	      } /** end check if excess moisture in top layer **/
	      
	      firstlayer=FALSE;
	      
	      /** verify that current layer moisture is greater than minimum **/
	      if ((moist[lindex]+ice[lindex]) < resid_moist[lindex]) {
		/** moisture cannot fall below residual moisture content **/
		Q12[lindex] += moist[lindex] - resid_moist[lindex];
		moist[lindex] = resid_moist[lindex];
	      }
	      
	      inflow = (Q12[lindex]+tmp_inflow);
	      Q12[lindex] += tmp_inflow;
	      
	      last_index++;
	      
	    } /* end loop through soil layers */
	    
	    /**************************************************
	    Compute Baseflow
	    **************************************************/
	    
	    /** ARNO model for the bottom soil layer (based on bottom
		soil layer moisture from previous time step) **/
	    
	    lindex = options.Nlayer-1;
	    Dsmax = soil_con->Dsmax / 24.;
	    
#if LINK_DEBUG
	    if(debug.DEBUG || debug.PRT_BALANCE) {
	      last_moist = moist[lindex];
	      /** Update debugging storage terms **/
	      debug.outflow[dist][band][lindex+1] += Q12[lindex-1];
	      debug.inflow[dist][band][lindex+2]  += Q12[lindex-1];
	    }
#endif // LINK_DEBUG
	    
	    /** Compute relative moisture **/
	    rel_moist = (moist[lindex]-resid_moist[lindex])
	      / (soil_con->max_moist[lindex]-resid_moist[lindex]);
	    
	    /** Compute baseflow as function of relative moisture **/
	    frac = Dsmax * soil_con->Ds / soil_con->Ws;
	    dt_baseflow = frac * rel_moist;
	    if (rel_moist > soil_con->Ws) {
	      frac = (rel_moist - soil_con->Ws) / (1 - soil_con->Ws);
	      dt_baseflow += Dsmax * (1 - soil_con->Ds / soil_con->Ws)
		* pow(frac,soil_con->c);
	    }
	    
	    /** Make sure baseflow isn't negative **/
	    if(dt_baseflow < 0) dt_baseflow = 0;
	    
	    /** Extract baseflow from the bottom soil layer **/ 
	    
	    moist[lindex] += Q12[lindex-1] - (evap[lindex][frost_area] + dt_baseflow);
	    
	    /** Check Lower Sub-Layer Moistures **/
	    tmp_moist = 0;
	    
	    /* If baseflow > 0 and soil moisture has gone below residual,
	     * take water out of baseflow and add back to soil to make up the difference
	     * Note: if baseflow is small, soil moisture may still be < residual after this */
	    if(dt_baseflow > 0 && (moist[lindex]+ice[lindex]) < resid_moist[lindex]) {
	      if ( dt_baseflow > resid_moist[lindex] - (moist[lindex]+ice[lindex]) ) {
		dt_baseflow -= resid_moist[lindex] - (moist[lindex]+ice[lindex]);
		moist[lindex] += resid_moist[lindex] - (moist[lindex]+ice[lindex]);
	      }
	      else {
		moist[lindex] += dt_baseflow;
		dt_baseflow = 0.0;
	      }
	    }
	    
	    if((moist[lindex]+ice[lindex]) > max_moist[lindex]) {
	      /* soil moisture above maximum */
	      tmp_moist = ((moist[lindex]+ice[lindex]) - max_moist[lindex]);
	      moist[lindex] = max_moist[lindex] - ice[lindex];
	      tmplayer = lindex;
	      while(tmp_moist > 0) {
		tmplayer--;
#if LINK_DEBUG
		if(debug.PRT_BALANCE) {
		  /** Update debugging storage terms **/
		  debug.inflow[dist][band][lindex+2]  -= tmp_moist;
		  debug.outflow[dist][band][lindex+1] -= tmp_moist;
		}
#endif // LINK_DEBUG
		if(tmplayer<0) {
		  /** If top layer saturated, add to runoff **/
		  runoff[frost_area] += tmp_moist;
		  tmp_moist = 0;
		}
		else {
		  /** else if sublayer exists, add excess soil moisture **/
		  moist[tmplayer] += tmp_moist ;
		  if ( ( moist[tmplayer] + ice[tmplayer]) 
		       > max_moist[tmplayer] ) {
		    tmp_moist = ((moist[tmplayer] + ice[tmplayer])
				 - max_moist[tmplayer]);
		    moist[tmplayer] = max_moist[tmplayer] - ice[tmplayer];
		  }
		  else tmp_moist=0;
		}
	      }
	    }
	    
	    baseflow[frost_area] += dt_baseflow;
	    
	  } /* end of hourly time step loop */
#if EXCESS_ICE
	}//end if subsidence did not occur or non-simple scenario for subsidence
#endif
	
	/*************************************************
	    Store runoff, baseflow and soil layer moisture
	*************************************************/
	
	if ( baseflow[frost_area] < 0 ) {
	  layer[lindex].evap   += baseflow[frost_area];
	  baseflow[frost_area]  = 0;
	}

	for ( lindex = 0; lindex < options.Nlayer; lindex++ ) 
#if SPATIAL_FROST
	  layer[lindex].moist += ((moist[lindex] + ice[lindex]) 
				  * frost_fract[frost_area]); 
	*runoff_ptr   += runoff[frost_area] * frost_fract[frost_area];
	*baseflow_ptr += baseflow[frost_area] * frost_fract[frost_area];
#else
	layer[lindex].moist = moist[lindex] + ice[lindex];      
	*runoff_ptr   += runoff[frost_area];
	*baseflow_ptr += baseflow[frost_area];
#endif // SPATIAL_FROST

#if LINK_DEBUG
	if(debug.PRT_BALANCE) {
	  debug.outflow[dist][band][options.Nlayer+2] = (runoff[frost_area] + baseflow[frost_area]);
	  debug.outflow[dist][band][options.Nlayer+1] = *baseflow;
	}
#endif // LINK_DEBUG
#if SPATIAL_FROST
      }
#endif // SPATIAL_FROST

    } /* if mu>0 */
  } /** Loop over wet and dry fractions **/

  /** Recompute Thermal Parameters Based on New Moisture Distribution **/
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    
    for(lindex=0;lindex<options.Nlayer;lindex++) {
      tmp_layer = find_average_layer(&(layer_wet[lindex]), 
				     &(layer_dry[lindex]), 
				     soil_con->depth[lindex], tmp_mu);
      moist[lindex] = tmp_layer.moist;
    }
    
#if EXCESS_ICE     
    if(SubsidenceUpdate == 0 ){
#endif
      ErrorFlag = distribute_node_moisture_properties(energy->moist, energy->ice,
						      energy->kappa_node, energy->Cs_node,
						      soil_con->dz_node, soil_con->Zsum_node, energy->T,
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
						      moist, soil_con->depth, 
						      soil_con->soil_density,
						      soil_con->bulk_density,
						      soil_con->quartz, Nnodes, 
						      options.Nlayer, soil_con->FS_ACTIVE);
      if ( ErrorFlag == ERROR ) return (ERROR);
#if EXCESS_ICE
    }
#endif
  }
  return (0);

}
