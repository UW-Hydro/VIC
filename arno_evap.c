#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/****************************************************************************
                                                                           
  ARNO/ARNO Model of Evaporation                                           
                                                                           
  Routine to compute evaporation based on the assumption that              
  evaporation is at the potential for the area which is saturated,         
  and at some percentage of the potential for the area which is partial    
  saturated.                                                               
                                                                           
  Evaporation from bare soil calculated only from uppermost layer.         
                                                                           
  Evaporation is in mm/(time step)  --> usually 1 day or 1 hour            

  modifications:
  6-8-2000  modified to make use of spatially distributed frost   KAC
       
****************************************************************************/

double arno_evap(layer_data_struct *layer_wet,
		 layer_data_struct *layer_dry,
		 double             rad,
		 double             air_temp,
		 double             vpd,
		 double             net_short,
		 double             D1,
		 double             max_moist,
		 double             elevation,
		 double             b_infilt,
		 double             displacement,
		 double             roughness,
		 double             ref_height,
		 double             ra,
		 double             delta_t,
		 double             mu,
#if SPATIAL_FROST
		 double             moist_resid,
		 double            *frost_fract)
#else
		 double             moist_resid)
#endif
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int    num_term;
  int    i;
  int    Ndist;
  int    dist;
#if SPATIAL_FROST
  int    frost_area;
#endif
  double tmp,beta_asp,dummy;
  double ratio,as;
  double Epot;		/* potential bare soil evaporation */
  double moist;
  double ice;
  double evap;
  double max_infil;
  double Evap;
  double tmpsum;
  layer_data_struct *layer;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;

  Evap = 0;

  for(dist=0;dist<Ndist;dist++) {

    if(dist>0) {
      mu = (1. - mu);
      layer = layer_dry;
    }
    else {
      layer = layer_wet;
    }

#if SPATIAL_FROST
    moist = 0;
    ice   = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      moist += (layer[0].moist - layer[0].ice[frost_area])
	* frost_fract[frost_area];
      ice   += layer[0].ice[frost_area] * frost_fract[frost_area];
    }
#else
    moist = layer[0].moist - layer[0].ice;
    ice   = layer[0].moist - layer[0].ice;
#endif
    if ( moist > max_moist ) moist = max_moist;

    if ( moist > 0 ) {
    
      /* Calculate the potential bare soil evaporation (mm/time step) */
    
      Epot = penman(rad, vpd, ra, 0.0, 0.0, 1.0, 1.0, 
	  	    air_temp, net_short, elevation, 0) * delta_t / SEC_PER_DAY;
    
      /**********************************************************************/
      /*  Compute temporary infiltration rate based on given soil_moist.    */
      /**********************************************************************/
      max_infil = (1.0 + b_infilt) * max_moist;
      if(b_infilt == -1.0)
        tmp = max_infil;
      else {
        ratio = 1.0 - (moist) / (max_moist);
        /*****if(ratio < SMALL && ratio > -SMALL) ratio = 0.;*****/
        if(ratio > 1.0) {
	  printf("\n  ERROR: SOIL RATIO GREATER THAN 1.0\n");
  	  printf("moisture %f   max_moisture %f -> ratio = %f\n",
	         moist,max_moist,ratio);
	  exit(0);
        }
        else {
	  if(ratio < 0.0) {
	    printf("\n  ERROR: SOIL RATIO LESS THAN 0.0\n");
	    printf("moisture %f   max_moisture %f -> ratio = %e\n",
		   moist,max_moist,ratio);
	    exit(0);
	  }
	  else
	    ratio = pow(ratio,(1.0 / (b_infilt + 1.0)));
        }
        tmp = max_infil*(1.0 - ratio);
      }
    
      /**********************************************************************/
      /* Evaporate at potential rate, i.e., Eq.(10) in Liang's derivation.  */
      /**********************************************************************/
    
      if(tmp >= max_infil)
        evap = Epot;     
      else {                 
      
        /********************************************************************/
        /*  Compute As. 'As' is % area saturated, '1-As' is % area          */
        /*  that is unsaturated.                                            */
        /********************************************************************/
      
        ratio = tmp/max_infil; 
        ratio = 1.0 - ratio;
      
        if(ratio > 1.0) {
	  printf("\n ARNO ERROR: EVAP RATIO GREATER THAN 1.0");
	  exit(0);
        }
        else {
	  if(ratio < 0.0) {
	    printf("\n ARNO ERROR: EVAP RATIO LESS THAN 0.0");
	    exit(0);
	  }
	  else {
	    if(ratio != 0.0)
	      ratio = pow(ratio,b_infilt);
	  }
        }
      
        as = 1 - ratio;
      
        /********************************************************************/
        /*  Compute the beta function in the ARNO evaporation model using   */
        /*  the first 30 terms in the power expansion expression.           */
        /********************************************************************/ 
      
        ratio = pow(ratio,(1.0/b_infilt));
      
        dummy = 1.0;
        for(num_term=1;num_term<=30;num_term++) {
	  tmpsum = ratio;
	  for ( i = 1; i < num_term; i++ ) tmpsum *= ratio;
	  dummy += b_infilt * tmpsum / (b_infilt + num_term);
        }
      
        beta_asp = as+(1.0-as)*(1.0-ratio)*dummy;
        evap = Epot*beta_asp;
      }
	
      /***********************************************************************/
      /*  Evaporation cannot exceed available soil moisture.                 */
      /*  Evaporation second soil layer = 0.0                                */
      /***********************************************************************/

      moist_resid *= D1 * 1000.;
      if ( ice > 0 ) {
        if ( ice > moist_resid ) {
 	  // ice content greater than wilting point can use all unfrozen moist
          if ( evap > moist ) evap = moist;
        }
        else {
	  // ice content less than wilting point restrict loss of unfrozen moist
          if ( evap > ( moist + ice ) - moist_resid )
            evap = ( moist + ice ) -  moist_resid;
        }
      }
      else {
        if ( evap > ( moist + ice ) - moist_resid)
          evap = ( moist + ice ) -  moist_resid;
      }

      layer[0].evap = evap;
      Evap += evap /1000. / delta_t * mu;

    }
    else {
      Evap = 0;
      layer[0].evap = 0;
    }

  }
  
  return(Evap);

}
