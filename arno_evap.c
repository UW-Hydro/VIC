#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************/
/*                                                                           */
/*  ARNO/ARNO Model of Evaporation                                           */
/*                                                                           */
/*  Routine to compute evaporation based on the assumption that              */
/*  evaporation is at the potential for the area which is saturated,         */
/*  and at some percentage of the potential for the area which is partial    */
/*  saturated.                                                               */
/*                                                                           */
/*  Evaporation from bare soil calculated only from uppermost layer.         */
/*                                                                           */
/*  Evaporation is in mm/(time step)  --> usually 1 day or 1 hour            */
/*****************************************************************************/

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
		 double             Tsurf,
		 double             displacement,
		 double             roughness,
		 double             ref_height,
		 double             ra,
		 double             dt,
		 double             mu)
{
  extern option_struct options;
  extern debug_struct debug;

  int    num_term;
  int    i;
  int    Ndist;
  int    dist;
  double tmp,beta_asp,dummy;
  double ratio,as;
  double Epot;		/* potential bare soil evaporation */
  double evap_temp;
  double moist;
  double evap;
  double moist_resid;
  double max_infil;
  double Evap;
  layer_data_struct *layer;

  evap_temp = Tsurf;

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
    
    moist = layer[0].moist_thaw*(layer[0].tdepth)/D1
          + layer[0].moist_froz 
          * (layer[0].fdepth-layer[0].tdepth)/D1
          + layer[0].moist*(D1-layer[0].fdepth)/D1;
    if(moist>max_moist) moist=max_moist;

    /* Calculate the potential bare soil evaporation (mm/time step) */
  
    Epot = penman(rad, vpd * 1000., ra, 0.0, 0.0, 1.0, 1.0, 
		  air_temp, evap_temp, net_short,
		  elevation, 0) * dt / 24.0;

    /**********************************************************************/
    /*  Compute temporary infiltration rate based on given soil_moist.    */
    /**********************************************************************/
    max_infil = (1.0+b_infilt)*max_moist;
    if(b_infilt == -1.0)
      tmp = max_infil;
    else {
      ratio = 1.0 - (moist) / (max_moist);
      /*****if(ratio < SMALL && ratio > -SMALL) ratio = 0.;*****/
      if(ratio > 1.0) {
	printf("\n  ERROR: SOIL RATIO GREATER THAN 1.0\n");
	printf("moisture %lf   max_moisture %lf -> ratio = %lf\n",
	       moist,max_moist,ratio);
	exit(0);
      }
      else {
	if(ratio < 0.0) {
	  printf("\n  ERROR: SOIL RATIO LESS THAN 0.0\n");
	  printf("moisture %lf   max_moisture %lf -> ratio = %le\n",
		 moist,max_moist,ratio);
	  exit(0);
	}
	else
	  ratio = pow(ratio,(1.0 / (b_infilt + 1.0)));
      }
      tmp = max_infil*(1.0 - ratio);
    }

    /************************************************************************/
    /* Evaporate at potential rate, i.e., Eq.(10) in Liang's derivation.    */
    /************************************************************************/

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
      for(num_term=1;num_term<=30;num_term++)
	dummy += b_infilt * pow(ratio,(double)num_term)/
	  (b_infilt + num_term);
    
      beta_asp = as+(1.0-as)*(1.0-ratio)*dummy;
      evap = Epot*beta_asp;
    }
	
    /***********************************************************************/
    /*  Evaporation cannot exceed available soil moisture.                 */
    /*  Evaporation second soil layer = 0.0                                */
    /***********************************************************************/

    if(options.FULL_ENERGY) moist_resid = MOIST_RESID;
    else moist_resid = 0.;

    moist = layer[0].moist_thaw*(layer[0].tdepth)/D1
          + layer[0].moist_froz 
          * (layer[0].fdepth-layer[0].tdepth)/D1
          + layer[0].moist*(D1-layer[0].fdepth)/D1;
    if((evap) > moist)
      evap = moist -  moist_resid * D1 * 1000.;

    layer[0].evap = evap;
    Evap += evap / 1000. / dt / 3600. * mu;

  }
  
  return(Evap);

}
