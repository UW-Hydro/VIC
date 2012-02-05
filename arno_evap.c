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
  2-24-03   moved unit conversion of moist_resid outside of distributed
            precipitation loop.  Moist_resid was being multiplied by
            D1 * 1000 twice for the dry fraction.                 KAC
  04-Jun-04 Moved unit conversion of moist_resid back inside distributed
	    precipitation loop in such a way that it does not get multiplied
	    by D1 * 1000 twice.					TJB
  04-Jun-04 Changed logic of evap limit check to avoid creating spurious
	    condensation.  Previously, when liquid moisture < residual
	    moisture, (liquid moisture - residual moisture) would be
	    negative.  Any non-negative evap would be greater than this,
	    resulting in evap getting set to (liquid moisture - residual
	    moisture), which would be negative (i.e. condensation).
	    This artificially created condensation in whatever amount
	    necessary to bring liquid moisture up to residual, causing
	    1) large latent heat flux, 2) incorrect surface temperatures,
	    3) occasional inability for calc_surf_energy_bal to converge
	    in root_brent, and 4) spuriously high runoff and baseflow.
	    Now there is an added condition that liquid moisture > residual
	    moisture for evap to be capped at (liquid moisture - residual
	    moisture).	In addition, the previous logic for capping evap
	    involved an incorrect calculation of the soil's ice content.
	    Since the new logic doesn't require calculation of ice content,
	    this calculation has been removed altogether.		TJB
  2007-Apr-06 Modified to handle grid cell errors by returning ERROR
	      status that can be trapped by calling routines.			GCT/KAC
  2009-Mar-16 Modified to use min_liq (minimum allowable liquid water
	      content) instead of resid_moist.  For unfrozen soil,
	      min_liq = resid_moist.						TJB
  2009-Jun-09 Moved computation of canopy resistance rc out of penman()
	      and into separate function calc_rc().				TJB
  2009-Nov-15 Changed D1 to depth1 to avoid confusion with the D1 used in
	      func_surf_energy_bal() (the parent function), which refers
	      to the depth of the 1st soil thermal node in some cases.		TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.				TJB
  2010-Apr-28 Removed net_short, displacement, roughness, and ref_height from
	      arg list as they are no longer used.				TJB
  2012-Jan-16 Removed LINK_DEBUG code						BN
****************************************************************************/

double arno_evap(layer_data_struct *layer_wet,
		 layer_data_struct *layer_dry,
		 double             rad,
		 double             air_temp,
		 double             vpd,
		 double             depth1,
		 double             max_moist,
		 double             elevation,
		 double             b_infilt,
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

    /* moist = liquid soil moisture */
#if SPATIAL_FROST
    moist = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      moist += (layer[0].moist - layer[0].ice[frost_area])
	* frost_fract[frost_area];
    }
#else
    moist = layer[0].moist - layer[0].ice;
#endif
    if ( moist > max_moist ) moist = max_moist;

    /* Calculate the potential bare soil evaporation (mm/time step) */
  
    Epot = penman(air_temp, elevation, rad, vpd, ra, 0.0, 0.0) * delta_t / SEC_PER_DAY;
  
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
	return( ERROR );
      }
      else {
	if(ratio < 0.0) {
	  printf("\n  ERROR: SOIL RATIO LESS THAN 0.0\n");
	  printf("moisture %f   max_moisture %f -> ratio = %e\n",
	         moist,max_moist,ratio);
	  return( ERROR );
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
	return ( ERROR );
      }
      else {
	if(ratio < 0.0) {
	  printf("\n ARNO ERROR: EVAP RATIO LESS THAN 0.0");
	  return( ERROR );
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

    /* only consider positive evaporation; we won't put limits on condensation */
    if (evap > 0.0) {
      if (moist > moist_resid * depth1 * 1000.) {
        /* there is liquid moisture available; cap evap at available liquid moisture */
        if (evap > moist -  moist_resid * depth1 * 1000.) {
          evap = moist -  moist_resid * depth1 * 1000.;
        }
      }
      else {
        /* no moisture available; cap evap at 0 */
        evap = 0.0;
      }
    }

    layer[0].evap = evap;
    Evap += evap / 1000. / delta_t * mu;

  }
  
  return(Evap);

}
