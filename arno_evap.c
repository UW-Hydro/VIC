#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

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

void arno_evap(atmos_data_struct *atmos, 
               layer_data_struct *layer,
               soil_con_struct    soil_con,
               double  Tsurf,
               double displacement,
               double roughness,
               double ref_height,
               double ra,
               global_param_struct global)
{
  extern option_struct options;
  extern debug_struct debug;

  int num_term;
  int i;
  double tmp,beta_asp,dummy;
  double ratio,as;
  double Epot;		/* potential bare soil evaporation */
  double evap_temp;
  double moist;
  double max_moist;
  double evap;
  double moist_resid;

  evap_temp = Tsurf;

  moist = layer[0].moist_thaw*(layer[0].tdepth)/soil_con.depth[0]
      + layer[0].moist_froz*(layer[0].fdepth-layer[0].tdepth)/soil_con.depth[0]
      + layer[0].moist*(soil_con.depth[0]-layer[0].fdepth)/soil_con.depth[0];
  max_moist = soil_con.max_moist[0];
  if(moist>max_moist) moist=max_moist;	/** correct for machine errors **/

  /* Calculate the potential bare soil evaporation (mm/time step) */
  
  Epot = penman(atmos->rad, atmos->vpd * 1000., ra, 0.0, 0.0, 1.0, 1.0, 
		atmos->air_temp, evap_temp, atmos->net_short,
		soil_con.elevation) * (double)global.dt / 24.0;

/**************************************************************************/
/*  Compute temporary infiltration rate based on given soil_moist.        */
/**************************************************************************/

  if(soil_con.b_infilt == -1.0)
    tmp = soil_con.max_infil;
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
	ratio = pow(ratio,(1.0 / (soil_con.b_infilt + 1.0)));
    }
    tmp = soil_con.max_infil*(1.0 - ratio);
  }

/****************************************************************************/
/* Evaporate at potential rate, i.e., Eq.(10) in Liang's derivation.        */
/****************************************************************************/

  if(tmp >= soil_con.max_infil)
    evap = Epot;     
  else {                 
    
/****************************************************************************/
/*  Compute As. 'As' is % area saturated, '1-As' is % area                  */
/*  that is unsaturated.                                                    */
/****************************************************************************/
    
    ratio = tmp/soil_con.max_infil; 
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
	  ratio = pow(ratio,soil_con.b_infilt);
      }
    }
    
    as = 1 - ratio;

/****************************************************************************/
/*  Compute the beta function in the ARNO evaporation model using the first */
/*  30 terms in the power expansion expression.                             */
/****************************************************************************/ 
	
    ratio = pow(ratio,(1.0/soil_con.b_infilt));
    
    dummy = 1.0;
    for(num_term=1;num_term<=30;num_term++)
      dummy += soil_con.b_infilt * pow(ratio,(double)num_term)/
	      (soil_con.b_infilt + num_term);
    
    beta_asp = as+(1.0-as)*(1.0-ratio)*dummy;
    evap = Epot*beta_asp;
  }
	
/***************************************************************************/
/*  Evaporation cannot exceed available soil moisture.                     */
/*  Evaporation second soil layer = 0.0                                    */
/***************************************************************************/

  layer[0].evap = evap;
  if(options.FULL_ENERGY) moist_resid = MOIST_RESID;
  else moist_resid = 0.;

  moist = layer[0].moist_thaw*(layer[0].tdepth)/soil_con.depth[0]
      + layer[0].moist_froz*(layer[0].fdepth-layer[0].tdepth)/soil_con.depth[0]
      + layer[0].moist*(soil_con.depth[0]-layer[0].fdepth)/soil_con.depth[0];
  if((layer[0].evap) > moist /*- soil_con.Wpwp[0] moist_resid * soil_con.depth[0] * 1000.*/ )
    layer[0].evap = moist - /*soil_con.Wpwp[0]*/ moist_resid * soil_con.depth[0] * 1000.;

}
