#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

void runoff(layer_data_struct *layer,
            energy_bal_struct *energy,
            soil_con_struct soil_con,
            double *runoff, double *baseflow,
	    double ppt, int dt,
            int Tlayer)
/**********************************************************************
	runoff.c	Keith Cherkauer		May 18, 1996

  This subroutine calculates infiltration and runoff from the surface,
  gravity driven drainage between all soil layers, and generates 
  baseflow from the bottom layer..
  
  Modifications:
  5/22/96	Routine modified to account for spatially varying
		precipitation, and it's effects on runoff.	KAC
  11/96		Code modified to account for extra model layers
		needed for frozen soils modeling.		KAC
  1/9/97	Infiltration and other rate parameters modified
		for time scales of less than 1 day.		KAC

  UNITS:	Ksat (mm/day)
		Q12  (mm/time step)
		moist (mm)
		inflow (mm)
                runoff (mm)

  Variables:
	ppt	incoming precipitation and snow melt
	mu	fraction of area that receives precipitation
	inflow	incoming water corrected for fractional area of precip (mu)

  NOTE: If not using distributed precipitation, x, and mu should be 
  set to 1.0.

  NOTE: Infiltration parameter should probably be adjusted for 
    temperature, since it affects Ksat, and interlayer flow.

  MODIFICATIONS:
    4-1-98 Soil moisture transport is now done on an hourly time
           step, irregardless to the model time step, to prevent
           numerical stabilities in the solution	Dag and KAC

**********************************************************************/
{  
  extern option_struct options;
  extern debug_struct debug;

  int firstlayer, lindex, sub;
  int num_it;
  int i;
  int *last_layer;
  int *last_sub;
  int *froz_solid;
  int last_index;
  int last_cnt;
  int tmp_index;
  int time_step;
  double ex, A, i_0, basis, frac;
  double inflow;
  double inte_f, ratio;
  double integral_a;
  double last_moist;
  double *resid_moist;
  double **submoist, **subice, **sublayer, **subT;
  double **submax_moist;
  double max_infil;
  double *Ksat;
  double **Q12;
  double *kappa;
  double *Cs;
  double *M;
  double Dsmax;
  double top_moist;
  double top_max_moist;
  double top_depth;
  double tmp_inflow;
  double dt_inflow, dt_outflow;
  double dt_runoff;

  resid_moist = (double *)calloc(options.Nlayer,sizeof(double));
  if(options.FULL_ENERGY) 
    for(i=0;i<options.Nlayer;i++) resid_moist[i] = soil_con.resid_moist[i];
  else for(i=0;i<options.Nlayer;i++) resid_moist[i] = 0.;

  last_layer   = (int *)calloc(options.Nlayer*3,sizeof(int));
  last_sub     = (int *)calloc(options.Nlayer*3,sizeof(int));
  froz_solid   = (int *)calloc(options.Nlayer*3,sizeof(int));
  Ksat         = (double *)calloc(options.Nlayer,sizeof(double));
  Q12          = (double **)calloc(options.Nlayer,sizeof(double *));
  submoist     = (double **)calloc(options.Nlayer,sizeof(double *));
  subice       = (double **)calloc(options.Nlayer,sizeof(double *));
  submax_moist = (double **)calloc(options.Nlayer,sizeof(double *));
  sublayer     = (double **)calloc(options.Nlayer,sizeof(double *));
  subT         = (double **)calloc(options.Nlayer,sizeof(double *));

  /** ppt = amount of liquid water coming to the surface **/
  inflow = ppt;

  /**************************************************
    Initialize Variables
  **************************************************/
  for(lindex=0;lindex<options.Nlayer;lindex++) {
    Ksat[lindex]         = soil_con.Ksat[lindex] / 24.;
    Q12[lindex]          = (double *)calloc(3,sizeof(double));
    submoist[lindex]     = (double *)calloc(3,sizeof(double));
    submax_moist[lindex] = (double *)calloc(3,sizeof(double));
    sublayer[lindex]     = (double *)calloc(3,sizeof(double));
    subT[lindex]         = (double *)calloc(3,sizeof(double));

    if(options.FROZEN_SOIL) {
      subT[lindex][0] = layer[lindex].T_thaw;
      subT[lindex][1] = layer[lindex].T_froz;
      subT[lindex][2] = layer[lindex].T;
    }
    sublayer[lindex][0]     = (layer[lindex].tdepth) / soil_con.depth[lindex];
    sublayer[lindex][1]     = (layer[lindex].fdepth - layer[lindex].tdepth) 
                            / soil_con.depth[lindex];
    sublayer[lindex][2]     = (soil_con.depth[lindex] - layer[lindex].fdepth) 
                            / soil_con.depth[lindex];
    submoist[lindex][0]     = layer[lindex].moist_thaw;
    submoist[lindex][1]     = layer[lindex].moist_froz;
    submoist[lindex][2]     = layer[lindex].moist;
    subice[lindex][0]       = 0.;
    subice[lindex][1]       = layer[lindex].ice;
    subice[lindex][2]       = 0.;
    submax_moist[lindex][0] = soil_con.max_moist[lindex];
    submax_moist[lindex][1] = soil_con.max_moist[lindex];
    submax_moist[lindex][2] = soil_con.max_moist[lindex];
    if(submoist[lindex][1]>submax_moist[lindex][1])
      submoist[lindex][1] = submax_moist[lindex][1];

  }

  /**************************************************
    Runoff Based on Soil Moisture Level of Top Two Layers
  **************************************************/
  top_moist = 0.;
  top_max_moist=0.;
  if(options.Nlayer>2) {
    if(options.FROZEN_SOIL) top_depth = soil_con.depth[0] + soil_con.depth[1];
    for(sub=0;sub<3;sub++) {
      top_moist     += (submoist[0][sub]+subice[0][sub])*sublayer[0][sub];
      top_max_moist += submax_moist[0][sub]*sublayer[0][sub];
      top_moist     += (submoist[1][sub]+subice[1][sub])*sublayer[0][sub];
      top_max_moist += submax_moist[1][sub]*sublayer[0][sub];
    }
  }
  else {
    for(sub=0;sub<3;sub++) {
      top_moist     += (submoist[0][sub]+subice[0][sub])*sublayer[0][sub];
      top_max_moist += submax_moist[0][sub]*sublayer[0][sub];
    }
  }
  if(top_moist>top_max_moist) top_moist=top_max_moist;
  if(layer[0].tdepth>0) sub=0;
  else if(layer[0].fdepth>0) sub=1;
  else sub=2;

  /**************************************************
    Calculate Runoff from Surface
  **************************************************/

  /** Runoff Calculations for Top Layer Only **/
  /** A and i_0 as in Wood et al. in JGR 97, D3, 1992 equation (1) **/

  max_infil = (1.0+soil_con.b_infilt) * top_max_moist;

  ex        = soil_con.b_infilt / (1.0 + soil_con.b_infilt);
  A         = 1.0 - pow((1.0 - top_moist / top_max_moist),ex);
  i_0       = max_infil * (1.0 - pow((1.0 - A),(1.0 
            / soil_con.b_infilt))); /* Maximum Inflow */

  /** equation (3a) Wood et al. **/
  
  if ( inflow == 0.0 ) *runoff = 0.0;
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

  /**************************************************
    Compute Flow Between Soil Layers (using an hourly time step)
  **************************************************/

  dt_inflow  = inflow / (double) dt;
  dt_runoff  = *runoff / (double) dt;
  dt_outflow = 0.0;

  for (time_step = 0; time_step < dt; time_step++) {
    inflow = dt_inflow;
    last_cnt = 0;

    for(lindex=0;lindex<options.Nlayer;lindex++) {
      for(sub=0;sub<3;sub++) {

        if(sublayer[lindex][sub] > 0. && (lindex!=options.Nlayer-1
            || (lindex==options.Nlayer-1 && sub!=2))) {

          /** Brooks & Corey relation for hydraulic conductivity **/
          /** If Saturated Moisture Content - Ice Content < 0.13, and
	      frozen layer is thicker than 5cm, layer 
              is assumed impermiable **/
          froz_solid[last_cnt] = FALSE;
          if(options.FROZEN_SOIL && sub==1 && 
	     (submax_moist[lindex][sub])/(soil_con.depth[lindex]*1000.) < 0.13 
	     && sublayer[lindex][1]*soil_con.depth[lindex]>0.05) {
            Q12[lindex][sub] = 0.0;
            froz_solid[last_cnt] = TRUE;
            if(last_cnt>0)
              Q12[last_layer[last_cnt-1]][last_sub[last_cnt-1]] = 0.;
          }
          else if(options.FROZEN_SOIL) {
	    if(sub==1)
	      Q12[lindex][sub] = modify_Ksat(subT[lindex][sub]) 
		             *  Ksat[lindex] 
                 	     * pow((submoist[lindex][sub])
                             / (soil_con.max_moist[lindex]),
                               soil_con.expt[lindex]); 
	    else {
	      Q12[lindex][sub] = modify_Ksat(subT[lindex][sub]) 
		             *  Ksat[lindex] 
                             * pow((submoist[lindex][sub] 
                             - resid_moist[lindex] 
                             * soil_con.depth[lindex] * 1000.)
                             / (soil_con.max_moist[lindex] 
                             - resid_moist[lindex] 
                             * soil_con.depth[lindex] * 1000.),
                               soil_con.expt[lindex]); 
              if(submoist[lindex][sub] <= 
                 resid_moist[lindex] * soil_con.depth[lindex] * 1000.)
                Q12[lindex][sub] = 0.;
            }
          }
          else {
            if(submoist[lindex][sub]
                > (resid_moist[lindex] * soil_con.depth[lindex] * 1000.)) {
              Q12[lindex][sub] = Ksat[lindex] * pow(((submoist[lindex][sub]
                               - resid_moist[lindex] * soil_con.depth[lindex] 
                               * 1000.) / (soil_con.max_moist[lindex] 
                               - resid_moist[lindex]
                               * soil_con.depth[lindex] * 1000.)),
                                 soil_con.expt[lindex]); 
            }
            else Q12[lindex][sub] = 0.;
          }
          last_layer[last_cnt] = lindex;
          last_sub[last_cnt] = sub;
          last_cnt++;
        }
      }
    }
    
    /**************************************************
      Solve for Current Soil Layer Moisture, and
      Check Versus Maximum and Minimum Moisture
      Contents.  
    **************************************************/
  
    firstlayer = TRUE;
    last_index = 0;
    for(lindex=0;lindex<options.Nlayer;lindex++) {
      for(sub=0;sub<3;sub++) {
  
        if(sublayer[lindex][sub] > 0. 
	   && !(lindex==options.Nlayer-1 && sub==2)) {
  
          if(debug.DEBUG || debug.PRT_BALANCE) 
	    last_moist = submoist[lindex][sub];
  
          tmp_inflow = 0.;
          if(firstlayer) {
            submoist[lindex][sub] = submoist[lindex][sub] + inflow
                / sublayer[lindex][sub] - dt_runoff
                / sublayer[lindex][sub] 
                - (Q12[lindex][sub] + layer[lindex].evap/(double)dt);
            if((submoist[lindex][sub]+subice[lindex][sub])
	       > submax_moist[lindex][sub] 
	       && !froz_solid[last_index+1]) {
              tmp_inflow = (submoist[lindex][sub]+subice[lindex][sub]) 
			 - submax_moist[lindex][sub];
              submoist[lindex][sub] = submax_moist[lindex][sub] 
		                    - subice[lindex][sub];
            }
	    else if((submoist[lindex][sub]+subice[lindex][sub])
		    > submax_moist[lindex][sub] 
		    && froz_solid[last_index+1]) {
              tmp_inflow=0.;
              *runoff += ((submoist[lindex][sub]+subice[lindex][sub]) 
                       - submax_moist[lindex][sub]) 
            	       * sublayer[lindex][sub];
              submoist[lindex][sub] = submax_moist[lindex][sub] 
                                    - subice[lindex][sub];
	    }
          }
          else {
            submoist[lindex][sub] = submoist[lindex][sub] + inflow
                                  / sublayer[lindex][sub] - (Q12[lindex][sub]
                                  + layer[lindex].evap/(double)dt);
            if((submoist[lindex][sub]+subice[lindex][sub]) 
	       > submax_moist[lindex][sub]
	       && !froz_solid[last_index+1]) {
              tmp_inflow = (submoist[lindex][sub]+subice[lindex][sub])
                         - submax_moist[lindex][sub];
              submoist[lindex][sub] = submax_moist[lindex][sub] 
		                    - subice[lindex][sub];
            } 
	    else if((submoist[lindex][sub]+subice[lindex][sub]) 
		    > submax_moist[lindex][sub] 
		    && froz_solid[last_index+1]) {
              tmp_inflow=0.;
              tmp_index = last_index;
              while(tmp_index>0 
		    && (submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
			+ subice[last_layer[tmp_index]][last_sub[tmp_index]])
		    > submax_moist[last_layer[tmp_index]][last_sub[tmp_index]]) {
		Q12[last_layer[tmp_index-1]][last_sub[tmp_index]] 
		  -= ((submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
                  + subice[last_layer[tmp_index]][last_sub[tmp_index]])
		  - submax_moist[last_layer[tmp_index]][last_sub[tmp_index]]) 
		  * sublayer[last_layer[tmp_index]][last_sub[tmp_index]]
		  / sublayer[last_layer[tmp_index-1]][last_sub[tmp_index]];
		submoist[last_layer[tmp_index-1]][last_sub[tmp_index]] 
		  += ((submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
                  + subice[last_layer[tmp_index]][last_sub[tmp_index]])
		  - submax_moist[last_layer[tmp_index]][last_sub[tmp_index]]) 
		  * sublayer[last_layer[tmp_index]][last_sub[tmp_index]]
		  / sublayer[last_layer[tmp_index-1]][last_sub[tmp_index]];
		submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
		  = submax_moist[last_layer[tmp_index]][last_sub[tmp_index]] 
                  - subice[last_layer[tmp_index]][last_sub[tmp_index]];
                tmp_index--;
              }
	      if(tmp_index==0 
		 && (submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
                 + subice[last_layer[tmp_index]][last_sub[tmp_index]])
		 > submax_moist[last_layer[tmp_index]][last_sub[tmp_index]]) {
		*runoff 
		  += ((submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
                  + subice[last_layer[tmp_index]][last_sub[tmp_index]])
		  - submax_moist[last_layer[tmp_index]][last_sub[tmp_index]])
		  * sublayer[last_layer[tmp_index]][last_sub[tmp_index]];
		submoist[last_layer[tmp_index]][last_sub[tmp_index]] 
		  = submax_moist[last_layer[tmp_index]][last_sub[tmp_index]]
                  - subice[last_layer[tmp_index]][last_sub[tmp_index]];
	      }
	    }
          }

          firstlayer=FALSE;

          /** check threshold values for current soil layer **/
	  if ((submoist[lindex][sub]+subice[lindex][sub])
	      < resid_moist[lindex] * soil_con.depth[lindex] * 1000.) {
	    /** moisture cannot fall below residual moisture content **/
	    Q12[lindex][sub] += submoist[lindex][sub] - resid_moist[lindex] 
	      * soil_con.depth[lindex] * 1000.;
	    submoist[lindex][sub] = resid_moist[lindex] 
	      * soil_con.depth[lindex] * 1000.;
	  }

          inflow = (Q12[lindex][sub]+tmp_inflow) * sublayer[lindex][sub];
  
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

    if(debug.PRT_BALANCE) {
      if(time_step==0) {
        for(lindex=0;lindex<options.Nlayer-1;lindex++) {
          sub=2;
          while(sublayer[lindex][sub]<=0. && sub>=0) sub--;
          debug.inflow[lindex+3] = Q12[lindex][sub]/sublayer[lindex][sub];
          if(lindex!=2)
            debug.outflow[lindex+2] = Q12[lindex][sub]/sublayer[lindex][sub];
        }
      }
      else {
        for(lindex=0;lindex<options.Nlayer-1;lindex++) {
          sub=2;
          while(sublayer[lindex][sub]<=0. && sub>=0) sub--;
          debug.inflow[lindex+3] += Q12[lindex][sub]/sublayer[lindex][sub];
          if(lindex!=2)
            debug.outflow[lindex+2] += Q12[lindex][sub]/sublayer[lindex][sub];
        }
      }
    }
  }
  inflow = dt_outflow;

  /**************************************************
    Compute Baseflow
  **************************************************/

  /** ARNO model for the bottom soil layer **/

  lindex = options.Nlayer-1;
  Dsmax = soil_con.Dsmax * (double)dt / 24.0;
  if(debug.DEBUG || debug.PRT_BALANCE) last_moist = submoist[lindex][2];

  frac = soil_con.Ds * Dsmax / (soil_con.Ws * soil_con.max_moist[lindex]);
  *baseflow = frac * submoist[lindex][2];
  if (submoist[lindex][2] > soil_con.Ws * soil_con.max_moist[lindex]) {
    frac = (submoist[lindex][2] - soil_con.Ws 
         * soil_con.max_moist[lindex]) / (soil_con.max_moist[lindex] 
         - soil_con.Ws * soil_con.max_moist[lindex]);
    *baseflow += (Dsmax - soil_con.Ds * Dsmax / soil_con.Ws) 
         * pow(frac,soil_con.c);
  }

  /** set threshold values for second soil layer **/ 

  submoist[lindex][2] += inflow / sublayer[lindex][2] - layer[lindex].evap
                      - *baseflow;

  /** Check Lower Layer Moisture **/
  if ((submoist[lindex][2]+subice[lindex][2]) < resid_moist[lindex] 
      * soil_con.depth[lindex] * 1000. ) {
    *baseflow  += submoist[lindex][2]
                - resid_moist[lindex] * soil_con.depth[lindex] * 1000.;
    submoist[lindex][2]  = resid_moist[lindex] * soil_con.depth[lindex] 
                         * 1000.;
    if(*baseflow<0.) {
      submoist[lindex][2] += *baseflow;
      *baseflow=0.;
    }
  }

  if((submoist[lindex][2]+subice[lindex][2]) > submax_moist[lindex][2]) {
    *baseflow += (submoist[lindex][2]+subice[lindex][2]) 
               - submax_moist[lindex][2];
    submoist[lindex][2] = submax_moist[lindex][2] - subice[lindex][2];
  }

  layer[lindex].moist = submoist[lindex][2];

  *baseflow *= sublayer[lindex][2];

  /** Recompute Thermal Parameters Based on New Moisture Distribution **/
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    kappa = NULL;
    Cs = NULL;
    M = NULL;
    soil_thermal_calc(soil_con,layer,*energy,kappa,Cs,M,M,M,options.Nlayer,
        Tlayer);

    /***** WARNING This will not work if dz or layers are changed *****/
    energy->kappa[0] = layer[0].kappa;
    energy->Cs[0] = layer[0].Cs;
    energy->kappa[1] = layer[1].kappa;
    energy->Cs[1] = layer[1].Cs;
  }

  for(lindex=0;lindex<options.Nlayer;lindex++) {
    free((char *)Q12[lindex]);
    free((char *)submoist[lindex]);
    free((char *)subice[lindex]);
    free((char *)submax_moist[lindex]);
    free((char *)sublayer[lindex]);
    free((char *)subT[lindex]);
  }
  free((char *)Ksat);
  free((char *)Q12);
  free((char *)submoist);
  free((char *)subice);
  free((char *)submax_moist);
  free((char *)sublayer);
  free((char *)subT);
  free((char *)last_layer);
  free((char *)last_sub);
  free((char *)froz_solid);
  free((char *)resid_moist);

  if(debug.PRT_BALANCE) {
    debug.inflow[2] = ppt - *runoff;
    debug.outflow[options.Nlayer+2] = *runoff + *baseflow;
    debug.outflow[options.Nlayer+1] = *baseflow;
  }
}
