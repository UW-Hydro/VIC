#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void runoff(layer_data_struct *layer_wet,
	    layer_data_struct *layer_dry,
            energy_bal_struct *energy,
            soil_con_struct    soil_con,
            double            *runoff_wet, 
	    double            *runoff_dry, 
	    double            *baseflow_wet,
	    double            *baseflow_dry,
	    double            *ppt, 
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
    5/22/96	Routine modified to account for spatially varying
		precipitation, and it's effects on runoff.	KAC
    11/96		Code modified to account for extra model layers
  		needed for frozen soils modeling.		KAC
    1/9/97	Infiltration and other rate parameters modified
		for time scales of less than 1 day.		KAC
    4-1-98 Soil moisture transport is now done on an hourly time
           step, irregardless to the model time step, to prevent
           numerical stabilities in the solution	Dag and KAC

**********************************************************************/
{  
  extern option_struct options;
  extern debug_struct  debug;

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
  double             ex, A, i_0, basis, frac;
  double             inflow;
  double             last_moist;
  double            *resid_moist;
  double             submoist[MAX_LAYERS][3];
  double             subice[MAX_LAYERS][3];
  double             sublayer[MAX_LAYERS][3];
  double             submax_moist[MAX_LAYERS][3];
  double             max_infil;
  double             Ksat[MAX_LAYERS];
  double             Q12[MAX_LAYERS][3];
  double            *kappa;
  double            *Cs;
  double            *M;
  double             Dsmax;
  double             top_moist;
  double             top_max_moist;
  double             bot_moist;
  double             tmp_inflow;
  double             tmp_moist;
  double             dt_inflow, dt_outflow;
  double             dt_runoff;
  double            *runoff;
  double            *baseflow;
  double             tmp_mu;
  layer_data_struct *layer;

  /** Set Residual Moisture **/
  resid_moist = (double *)calloc(options.Nlayer,sizeof(double));
  if(options.FULL_ENERGY) 
    for(i=0;i<options.Nlayer;i++) resid_moist[i] = soil_con.resid_moist[i] 
				    * soil_con.depth[i] * 1000.;
  else for(i=0;i<options.Nlayer;i++) resid_moist[i] = 0.;

  /** Initialize Other Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  tmp_mu = mu;
  
  /** Allocate and Set Values for Soil Sublayers **/
  for(dist=0;dist<Ndist;dist++) {

    if(dist>0) {
      layer    = layer_dry;
      runoff   = runoff_dry;
      baseflow = baseflow_dry;
      mu       = (1. - mu);
    }
    else {
      layer    = layer_wet;
      runoff   = runoff_wet;
      baseflow = baseflow_wet;
    }

    if(mu>0.) {
      
      /** ppt = amount of liquid water coming to the surface **/
      inflow = ppt[dist];
      
      /**************************************************
	Initialize Variables
      **************************************************/
      for(lindex=0;lindex<options.Nlayer;lindex++) {
	Ksat[lindex]         = soil_con.Ksat[lindex] / 24.;

	/** Set Sublayer Fraction **/
	sublayer[lindex][0]     
	  = (layer[lindex].tdepth) / soil_con.depth[lindex];
	sublayer[lindex][1]     
	  = (layer[lindex].fdepth - layer[lindex].tdepth) 
	  / soil_con.depth[lindex];
	sublayer[lindex][2]     
	  = (soil_con.depth[lindex] - layer[lindex].fdepth) 
	  / soil_con.depth[lindex];

	/** Set Sublayer Unfrozen Moisture Content **/
	submoist[lindex][0] = layer[lindex].moist_thaw;
	if(submoist[lindex][0]<0) {
	  sprintf(ErrStr,
		  "Layer %i thawed sublayer has negative soil moisture, %lf", 
		  lindex, submoist[lindex][0]);
	  vicerror(ErrStr);
	}
	submoist[lindex][1] = layer[lindex].moist_froz;
	if(submoist[lindex][1]<0) {
	  sprintf(ErrStr,
		  "Layer %i frozen sublayer has negative soil moisture, %lf", 
		  lindex, submoist[lindex][1]);
	  vicerror(ErrStr);
	}
	submoist[lindex][2] = layer[lindex].moist;
	if(submoist[lindex][2]<0) {
	  sprintf(ErrStr,
		  "Layer %i unfrozen sublayer has negative soil moisture, %lf", 
		  lindex, submoist[lindex][2]);
	  vicerror(ErrStr);
	}

	/** Set Sublayer Ice Content **/
	subice[lindex][0]       = 0.;
	subice[lindex][1]       = layer[lindex].ice;
	subice[lindex][2]       = 0.;

	/** Set Sublayer Maximum Moisture Content **/
	submax_moist[lindex][0] = soil_con.max_moist[lindex];
	submax_moist[lindex][1] = soil_con.max_moist[lindex];
	submax_moist[lindex][2] = soil_con.max_moist[lindex];
	
      }
      
      /********************************************************************
        Store original soil moisture for bottom layers for baseflow calcs
      ********************************************************************/
      bot_moist = 0;
      lindex    = options.Nlayer-1;
      for(sub=0;sub<3;sub++) {
	bot_moist += (submoist[lindex][sub]) * sublayer[lindex][sub];
      }
      
      /******************************************
	Check if top soil layer is frozen solid
      ******************************************/
      *runoff = -99;
/*****
      if(options.FROZEN_SOIL &&  
	 (submax_moist[0][1])/(soil_con.depth[0]*1000.) < 0.13
	 && sublayer[0][1]*soil_con.depth[0]>0.05 && sublayer[0][0]==0) {
	*runoff = inflow;
      }
*****/

      /**********************************************************
        Check if thawed layer exists above solidly frozen layer
      **********************************************************/
/*****
      else if(options.FROZEN_SOIL &&  
	      (submax_moist[0][1])/(soil_con.depth[0]*1000.) < 0.13
	      && sublayer[0][1]*soil_con.depth[0]>0.05 && sublayer[0][0]>0) {
	top_moist = submoist[0][0] * sublayer[0][0];
	top_max_moist = submax_moist[0][0] * sublayer[0][0];
      }
*****/

      /******************************************************
        Runoff Based on Soil Moisture Level of Upper Layers
      ******************************************************/

/*****
      else {
*****/
	top_moist = 0.;
	top_max_moist=0.;
	for(lindex=0;lindex<options.Nlayer-1;lindex++) {
	  for(sub=0;sub<3;sub++) {
	    top_moist += (submoist[lindex][sub] + subice[lindex][sub]) 
	      * sublayer[lindex][sub];
	    top_max_moist += submax_moist[lindex][sub] * sublayer[lindex][sub];
	  }
	}
	if(top_moist>top_max_moist) top_moist = top_max_moist;
	if(layer[0].tdepth>0) sub = 0;
	else if(layer[0].fdepth>0) sub = 1;
	else sub = 2;
/*****
      }
*****/

      /**************************************************
        Calculate Runoff from Surface
      **************************************************/
      
      /** Runoff Calculations for Top Layer Only **/
      /** A and i_0 as in Wood et al. in JGR 97, D3, 1992 equation (1) **/
      
      if(*runoff < 0) {
	max_infil = (1.0+soil_con.b_infilt) * top_max_moist;
	
	ex        = soil_con.b_infilt / (1.0 + soil_con.b_infilt);
	A         = 1.0 - pow((1.0 - top_moist / top_max_moist),ex);
	i_0       = max_infil * (1.0 - pow((1.0 - A),(1.0 
						      / soil_con.b_infilt))); 
	/* Maximum Inflow */
	
	/** equation (3a) Wood et al. **/
	
	if (inflow == 0.0) *runoff = 0.0;
	else if (max_infil == 0.0) *runoff = inflow;
	else if ((i_0 + inflow) > max_infil) 
	  *runoff = inflow - top_max_moist + top_moist;
	
	/** equation (3b) Wood et al. (wrong in paper) **/
	else {
	  basis = 1.0 - (i_0 + inflow) / max_infil;
	  *runoff = (inflow - top_max_moist + top_moist + top_max_moist
		     * pow(basis,1.0*(1.0+soil_con.b_infilt)));
	}
	if(*runoff<0.) *runoff=0.;
      }
      
      /**************************************************
        Compute Flow Between Soil Layers (using an hourly time step)
      **************************************************/
      
      dt_inflow  =  inflow / (double) dt;
      dt_outflow =  0.0;
      for(sub=0;sub<3;sub++) 
	if(sublayer[options.Nlayer-1][sub]>0) storesub = sub;
      
      for (time_step = 0; time_step < dt; time_step++) {
	inflow   = dt_inflow;
	last_cnt = 0;

	/*************************************
          Compute Drainage between Sublayers 
        *************************************/
	
	for(lindex=0;lindex<options.Nlayer;lindex++) {
	  for(sub=0;sub<3;sub++) {
	    
	    if(sublayer[lindex][sub] > 0.) {
	      
	      /** Brooks & Corey relation for hydraulic conductivity **/
	      /** If Saturated Moisture Content - Ice Content < 0.13, and
		  frozen layer is thicker than 5cm, layer 
		  is assumed impermiable **/
	      froz_solid[last_cnt] = FALSE;

	      if((tmp_moist=submoist[lindex][sub]-layer[lindex].evap) 
		 < resid_moist[lindex])
		tmp_moist = resid_moist[lindex];

	      if(submoist[lindex][sub] > resid_moist[lindex]) {
		Q12[lindex][sub] 
		  = Ksat[lindex] * pow(((tmp_moist
					 - resid_moist[lindex]) 
					/ (soil_con.max_moist[lindex] 
					   - resid_moist[lindex])),
				       soil_con.expt[lindex]); 
	      }
	      else Q12[lindex][sub] = 0.;
	      last_layer[last_cnt] = lindex;
	      last_sub[last_cnt] = sub;
	      last_cnt++;
	    }
	  }
	}

	/*************************************************
	  Stop drainage from bottom sublayer
        *************************************************/
	Q12[options.Nlayer-1][storesub] = 0.;
	
	/**************************************************
          Solve for Current Soil Layer Moisture, and
          Check Versus Maximum and Minimum Moisture
          Contents.  
	**************************************************/
	
	firstlayer = TRUE;
	last_index = 0;
	for(lindex=0;lindex<options.Nlayer;lindex++) {

	  if(lindex==0) dt_runoff = *runoff / (double) dt;

	  /* Store moistures for water balance debugging */
	  if(debug.PRT_BALANCE) {
	    if(time_step==0) {
	      if(firstlayer)
		debug.inflow[dist][band][lindex+2] = inflow - dt_runoff;
	      else {
		debug.inflow[dist][band][lindex+2] = inflow;
		debug.outflow[dist][band][lindex+1] = inflow;
	      }
	    }
	    else {
	      if(firstlayer)
		debug.inflow[dist][band][lindex+2]  += inflow - dt_runoff;
	      else {
		debug.inflow[dist][band][lindex+2]  += inflow;
		debug.outflow[dist][band][lindex+1] += inflow;
	      }
	    }
	  }

	  /* transport moisture for all sublayers **/
	  for(sub=0;sub<3;sub++) {
	    
	    if(sublayer[lindex][sub] > 0.) {
	      
	      if(debug.DEBUG || debug.PRT_BALANCE) 
		last_moist = submoist[lindex][sub];
	      
	      tmp_inflow = 0.;

	      /** Update soil layer moisture content **/
	      submoist[lindex][sub] = submoist[lindex][sub] 
		+ (inflow - dt_runoff)
		/ sublayer[lindex][sub] - (Q12[lindex][sub]
					   + layer[lindex].evap/(double)dt);

	      dt_runoff = 0;

	      /** Verify that soil sublayer moisture content falls within 
		  acceptable range **/
	      if((submoist[lindex][sub]+subice[lindex][sub]) 
		 > submax_moist[lindex][sub]
		 && (!froz_solid[last_index+1])) {
		tmp_inflow = (submoist[lindex][sub]+subice[lindex][sub])
		  - submax_moist[lindex][sub];
		submoist[lindex][sub] = submax_moist[lindex][sub] 
		  - subice[lindex][sub];

		/* Layer is thick, so excess soil moisture is trapped above 
                   (except for top thin layer, where excess is included in infiltration) **/
		if(lindex>0 && sublayer[lindex][sub]*soil_con.depth[lindex] > 0.10) {
		  tmp_inflow *= sublayer[lindex][sub];
		  tmpsub   = sub;
		  tmplayer = lindex;
		  while(tmp_inflow > 0) {
		    tmpsub--;
		    if(tmpsub<0) {
		      tmplayer--;
		      tmpsub=2;
		      if(debug.PRT_BALANCE) {
			/** Update debugging storage terms **/
			debug.inflow[dist][band][lindex+2]  -= tmp_inflow;
			debug.outflow[dist][band][lindex+1] -= tmp_inflow;
		      }
		    }
		    if(tmplayer<0) {
		      /** If top layer saturated, add to top layer **/
		      *runoff += tmp_inflow;
		      tmp_inflow = 0;
		    }
		    else if(sublayer[tmplayer][tmpsub] > 0) {
		      /** else if sublayer exists, add excess soil moisture **/
		      submoist[tmplayer][tmpsub] += tmp_inflow 
			/ sublayer[tmplayer][tmpsub];
		      if((submoist[tmplayer][tmpsub]+subice[tmplayer][tmpsub]) 
			 > submax_moist[tmplayer][tmpsub]) {
			tmp_inflow = ((submoist[tmplayer][tmpsub]
				      +subice[tmplayer][tmpsub])
				     - submax_moist[tmplayer][tmpsub]) 
			  * sublayer[tmplayer][tmpsub];
			submoist[tmplayer][tmpsub] 
			  = submax_moist[tmplayer][tmpsub] 
			  - subice[tmplayer][tmpsub];
		      }
		      else tmp_inflow=0;
		    }
		  }
		} /** end trapped excess moisture **/
	      } 
	      
	      firstlayer=FALSE;
	      
	      /** check threshold values for current soil layer **/
	      if ((submoist[lindex][sub]+subice[lindex][sub])
		  < resid_moist[lindex]) {
		/** moisture cannot fall below residual moisture content **/
		Q12[lindex][sub]
		  += submoist[lindex][sub] - resid_moist[lindex];
		submoist[lindex][sub] = resid_moist[lindex];
	      }
	      
	      inflow = (Q12[lindex][sub]+tmp_inflow) * sublayer[lindex][sub];
	      Q12[lindex][sub] += tmp_inflow;
	      
	      last_index++;
	      
	    }
	  }
	  
	  if(sublayer[lindex][0] > 0)
	    layer[lindex].moist_thaw = submoist[lindex][0];
	  else layer[lindex].moist_thaw = 0.;
	  
	  if(sublayer[lindex][1] > 0)
	    layer[lindex].moist_froz = submoist[lindex][1];
	  else layer[lindex].moist_froz = 0.;
	  
	  if(sublayer[lindex][2] > 0)
	    layer[lindex].moist = submoist[lindex][2];
	  else layer[lindex].moist = 0.;

	}
	dt_outflow += inflow;
	
      } /** End hourly time step loop **/
      inflow = dt_outflow;
      
      /**************************************************
        Compute Baseflow
      **************************************************/
      
      /** ARNO model for the bottom soil layer **/
      
      lindex = options.Nlayer-1;
      submoist[lindex][storesub] += inflow / sublayer[lindex][storesub];
      Dsmax = soil_con.Dsmax * (double)dt / 24.0;
      if(debug.DEBUG || debug.PRT_BALANCE) last_moist = submoist[lindex][2];
      
      frac = soil_con.Ds * Dsmax / (soil_con.Ws * soil_con.max_moist[lindex]);
      *baseflow = frac * bot_moist;
      if (bot_moist > soil_con.Ws * soil_con.max_moist[lindex]) {
	frac = (bot_moist - soil_con.Ws 
		* soil_con.max_moist[lindex]) / (soil_con.max_moist[lindex] 
						 - soil_con.Ws 
						 * soil_con.max_moist[lindex]);
	*baseflow += (Dsmax - soil_con.Ds * Dsmax / soil_con.Ws) 
	  * pow(frac,soil_con.c);
      }
      if(*baseflow < 0) *baseflow = 0;
      
      /** Extract basflow from all sublayers in the bottom soil layer **/ 
      
      bot_moist = 0;
      for(sub=0;sub<3;sub++) 
	bot_moist += (submoist[lindex][sub]) * sublayer[lindex][sub];
      for(sub=0;sub<3;sub++) 
	if(sublayer[lindex][sub] > 0)
	  submoist[lindex][sub] -= (*baseflow 
				    * (submoist[lindex][sub] 
				       * sublayer[lindex][sub]) 
				    / bot_moist) / sublayer[lindex][sub];
      
      /** Check Lower Sub-Layer Moistures **/
      tmp_moist = 0;

      for(sub=1;sub<3;sub++) {

	if (sublayer[lindex][sub] > 0) {
	  if((submoist[lindex][sub]+subice[lindex][sub]) 
	     < resid_moist[lindex]) {
	    /* Not enough moisture left */
	    *baseflow  += ((submoist[lindex][sub]
			    - resid_moist[lindex]) * sublayer[lindex][sub]);
	    submoist[lindex][sub]  = resid_moist[lindex];
	  }

	  if((submoist[lindex][sub]+subice[lindex][sub]) 
	     > submax_moist[lindex][sub]) {
	    /* Too much soil moisture */
	    tmp_moist = ((submoist[lindex][sub]+subice[lindex][sub]) 
			  - submax_moist[lindex][sub]) * sublayer[lindex][sub];
	    submoist[lindex][sub] = submax_moist[lindex][sub] 
	      - subice[lindex][sub];
	    tmpsub   = sub;
	    tmplayer = lindex;
	    while(tmp_moist > 0) {
	      tmpsub--;
	      if(tmpsub<0) {
		tmplayer--;
		tmpsub=2;
		if(debug.PRT_BALANCE) {
		  /** Update debugging storage terms **/
		  debug.inflow[dist][band][lindex+2]  -= tmp_moist;
		  debug.outflow[dist][band][lindex+1] -= tmp_moist;
		}
	      }
	      if(tmplayer<0) {
		/** If top layer saturated, add to top layer **/
		*runoff += tmp_moist;
		tmp_moist = 0;
	      }
	      else if(sublayer[tmplayer][tmpsub] > 0) {
		/** else if sublayer exists, add excess soil moisture **/
		submoist[tmplayer][tmpsub] += tmp_moist 
		  / sublayer[tmplayer][tmpsub];
		if((submoist[tmplayer][tmpsub]+subice[tmplayer][tmpsub]) 
		   > submax_moist[tmplayer][tmpsub]) {
		  tmp_moist = ((submoist[tmplayer][tmpsub]
				+subice[tmplayer][tmpsub])
		    - submax_moist[tmplayer][tmpsub]) 
		    * sublayer[tmplayer][tmpsub];
		  submoist[tmplayer][tmpsub] = submax_moist[tmplayer][tmpsub] 
		    - subice[tmplayer][tmpsub];
		}
		else tmp_moist=0;
	      }
	    }
	  }
	}

      }
      
      for(lindex=0;lindex<options.Nlayer;lindex++) {

	layer[lindex].moist_thaw = submoist[lindex][0];
	layer[lindex].moist_froz = submoist[lindex][1];
	layer[lindex].moist      = submoist[lindex][2];

      }
      
      if(debug.PRT_BALANCE) {
	debug.outflow[dist][band][options.Nlayer+2] = (*runoff + *baseflow);
	debug.outflow[dist][band][options.Nlayer+1] = *baseflow;
      }
    }
  } /** Loop over wet and dry fractions **/

  /** Recompute Thermal Parameters Based on New Moisture Distribution **/
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    layer = (layer_data_struct *)calloc(options.Nlayer,
					sizeof(layer_data_struct));
    for(lindex=0;lindex<options.Nlayer;lindex++)
      layer[lindex] = find_average_layer(layer_wet[lindex], layer_dry[lindex],
					 soil_con.depth[lindex], tmp_mu);
    kappa = NULL;
    Cs = NULL;
    M = NULL;
    soil_thermal_calc(soil_con,layer,*energy,kappa,Cs,M,options.Nlayer,
		      Nnodes);
    
    /***** WARNING This will not work if dz or layers are changed *****/
    for(lindex=0;lindex<options.Nlayer;lindex++) {
      layer_wet[lindex].kappa = layer[lindex].kappa;
      layer_wet[lindex].Cs    = layer[lindex].Cs;
      if(Ndist>1) {
	layer_dry[lindex].kappa = layer[lindex].kappa;
	layer_dry[lindex].Cs    = layer[lindex].Cs;
      }
    }
    energy->kappa[0] = layer[0].kappa;
    energy->Cs[0] = layer[0].Cs;
    energy->kappa[1] = layer[1].kappa;
    energy->Cs[1] = layer[1].Cs;
    free((char*)layer);
  }

  free((char *)resid_moist);

}
