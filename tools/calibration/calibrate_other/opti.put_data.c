#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: put_data.c,v 2.9 1998/10/02 01:34:50 vicadmin Exp $";

void put_data(dist_prcp_struct  *prcp,
	      atmos_data_struct *atmos,
	      veg_con_struct    *veg_con,
              outfiles_struct    outfiles,
              double            *depth,
	      double            *AreaFract,
	      dmy_struct        *dmy,
              int                rec,
	      int                dt)
/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed presipitation data structures KAC

**********************************************************************/
{
  extern veg_lib_struct  *veg_lib;
  extern option_struct    options;
  extern debug_struct     debug;

  static int day_cnt;
  static double           runoff_sum, baseflow_sum;

  out_data_struct        *out_data;

  int                     veg, vegnum;
  int                     index;
  int                     Ndist;
  int                     dist;
  int                     band;
  int                     Nbands;
  double                  mu;
  double                  tmp_evap;
  double                  tmp_moist;
  double                  tmp_ice;
  double                  rad_temp = 0.0;
  double                  inflow, outflow, storage;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(rec==0) { 
    day_cnt = 0; 
    runoff_sum = baseflow_sum = 0;
  }
  
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = options.SNOW_BAND;

  out_data = (out_data_struct *) calloc(1,sizeof(out_data_struct)); 
  /*************************************************
    Store Output for Precipitation Distribution Type
    *************************************************/

  cell    = prcp[0].cell;
  veg_var = prcp[0].veg_var;
  snow    = prcp[0].snow;
  energy  = prcp[0].energy;
  
  /**************************************
    Store Output for all Vegetation Types
    **************************************/
  for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++) {
    
    /** record total evaporation **/
    for ( dist = 0; dist < Ndist; dist++ ) {
      if(dist==0) 
	mu = prcp[0].mu[veg];
      else 
	mu = 1. - prcp[0].mu[veg];

      for(band=0;band<Nbands;band++) {
	if(AreaFract[band] > 0.) {

	  
	  /** record runoff **/
	  out_data->runoff   += cell[dist][veg][band].runoff 
	    * veg_con[veg].Cv * mu * AreaFract[band];
	  
	  /** record baseflow **/
	  out_data->baseflow += cell[dist][veg][band].baseflow 
	    * veg_con[veg].Cv * mu * AreaFract[band]; 
	  
	}
      }
    }
    
  }
  
  /***********************
    Store Bare Soil Output
    ***********************/

  /** record evaporation for bare soil **/
  vegnum = veg_con[0].vegetat_type_num;
  for ( dist = 0; dist < Ndist; dist++ ) {
    if(dist==0) 
      mu = prcp[0].mu[veg];
    else 
      mu = 1. - prcp[0].mu[veg];

    for(band=0;band<Nbands;band++) {
      if(AreaFract[band]>0.) {
	
	/** record runoff for bare soil **/
	out_data->runoff   += cell[dist][vegnum][band].runoff 
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
	/** record baseflow for bare soil **/
	out_data->baseflow += cell[dist][vegnum][band].baseflow 
	  * (1.0 - veg_con[0].Cv_sum) * mu * AreaFract[band];
	
      }
    }
  }


  day_cnt++;
  if(day_cnt==24/dt) {
    out_data->runoff += runoff_sum;
    out_data->baseflow += baseflow_sum;
    write_data(out_data, outfiles, dmy, dt);
    day_cnt=0;
    runoff_sum = baseflow_sum = 0;
  }
  else {
    runoff_sum += out_data->runoff;
    baseflow_sum += out_data->baseflow;
  }

   free((char *)out_data);  
}
