#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void dist_prec(atmos_data_struct   *atmos,
               dist_prcp_struct    *prcp,
               soil_con_struct     *soil_con,
               veg_con_struct      *veg_con,
               dmy_struct          *dmy,
               global_param_struct *global_param,
               outfiles_struct     *outfiles,
               int                  rec,
               int                  cellnum,
               char                 NEWCELL,
               char                 LASTREC) {
/**********************************************************************
  dist_prec		Keith Cherkauer		October 9, 1997

  This subroutine calls the solution routines for a single grid cell
  for one time step.  It also controls the distribution of precipitation
  and reassembles grid cell data for output.

  The fractional coverage of precipitation over an area or grid cell, 
  mu, is estimated using the equation from Fan et. al. (1996).  The 
  coefficient, 0.6, was selected for the Arkansas - Red River Basin and
  was found using precipitation records on a 100km x 100km area.  It
  may not be applicable to all regions, please check the reference

  References:

  Modifications:
  11-30-98 Added counter to assure that a storm has been stopped
           for at least one day, before allowing the model to 
	   average soil moisture when a new precipitation event
	   arrives.                                             KAC

**********************************************************************/

  extern option_struct   options;
  extern veg_lib_struct *veg_lib; 
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  static char STILL_STORM;
  static int  DRY_TIME;

  char    ANY_SNOW;
  int     veg, i;
  int     month;
  double  Wdmax;
  double  NEW_MU;

#if SAVE_STATE

  /************************************
    Save model state at assigned date
  ************************************/

  if ( dmy[rec].hour == 0 && dmy[rec].year == global_param->stateyear
       && dmy[rec].month == global_param->statemonth 
       && dmy[rec].day == global_param->stateday )
    write_model_state(prcp, global_param, veg_con[0].vegetat_type_num, 
		      soil_con->gridcel, outfiles, soil_con);

#endif

  if(options.DIST_PRCP) {

    /*******************************************
      Controls Distributed Precipitation Model
    *******************************************/
     
    NEW_MU = 1.0 - exp(-options.PREC_EXPT*atmos->prec[NR]);
    for(veg=0; veg<=veg_con[0].vegetat_type_num; veg++) {
      ANY_SNOW = FALSE;
      for(i=0; i<options.SNOW_BAND; i++)
        /* Check for snow on ground or falling */
	if(prcp->snow[veg][i].swq > 0 
	   || prcp->snow[veg][i].snow_canopy > 0.) 
	  ANY_SNOW = TRUE;
      if(ANY_SNOW || atmos->snowflag[NR]) {
        /* If snow present, mu must be set to 1. */
	NEW_MU = 1.;
	if(rec == 0) {
          /* Set model variables if first time step */
	  prcp->mu[veg]=NEW_MU;
	  if(atmos->prec[NR] > 0) 
	    STILL_STORM=TRUE;
	  else 
	    STILL_STORM=FALSE;
          DRY_TIME = 0;
	} 
	ANY_SNOW = TRUE;
      }
      else {
	if(rec==0) {
	  if(atmos->prec[NR] == 0) {
	    /* If first time step has no rain, than set mu to 1. */
	    prcp->mu[veg] = 1.;
	    NEW_MU=1.;
	    STILL_STORM = TRUE;
	    DRY_TIME = 24;
	  }
	  else {
	    /* If first time step has rain, then set mu based on intensity */
	    prcp->mu[veg]=NEW_MU;
	    STILL_STORM=TRUE;
	    DRY_TIME = 0;
	  }
	}
	else if(atmos->prec[NR] == 0 
		&& DRY_TIME >= 24./(float)global_param->dt) {
          /* Check if storm has ended */
	  NEW_MU=prcp->mu[veg];
	  STILL_STORM=FALSE;
          DRY_TIME = 0;
	}
        else if(atmos->prec[NR] == 0) {
	  /* May be pause in storm, keep track of pause length */
	  NEW_MU=prcp->mu[veg];
	  DRY_TIME += global_param->dt;
	}
      }

      if(!STILL_STORM && (atmos->prec[NR] > STORM_THRES || ANY_SNOW)) {
	/** Average soil moisture before a new storm **/
	initialize_new_storm(prcp->cell,prcp->veg_var,
			     veg,veg_con[0].vegetat_type_num,rec,
			     prcp->mu[veg],NEW_MU);
	STILL_STORM=TRUE;
	prcp->mu[veg] = NEW_MU;
      }
      else if(NEW_MU != prcp->mu[veg] && STILL_STORM) {
	/** Redistribute soil moisture during the storm if mu changes **/
	if ( dmy[rec].day == 1 && dmy[rec].hour == 0 ) {
	  month = dmy[rec].month - 2;
	  if ( month < 0 ) month = 11;
	}
	else month = dmy[rec].month - 1;
	if (veg < veg_con[0].vegetat_type_num) 
	  Wdmax = veg_lib[veg_con[veg].veg_class].Wdmax[month];
	else 
	  Wdmax = 0;
	redistribute_during_storm(prcp->cell, prcp->veg_var, veg, 
				  veg_con[0].vegetat_type_num, rec, Wdmax, 
				  prcp->mu[veg], NEW_MU, soil_con->max_moist);
	prcp->mu[veg] = NEW_MU;
      }
    }

    /** Solve model time step **/
    full_energy(rec, atmos, soil_con, veg_con, prcp, dmy, global_param, 
		cellnum, NEWCELL);

  }

  else {

    /**************************************************
      Controls Grid Cell Averaged Precipitation Model
    **************************************************/

    full_energy(rec, atmos, soil_con, veg_con, prcp, dmy, global_param, 
		cellnum, NEWCELL);

  }

  /**************************************************
    Write cell average values for current time step
  **************************************************/

  put_data(prcp, atmos, veg_con, outfiles, soil_con->depth, 
	   soil_con->dz_node, soil_con->dp, soil_con->AreaFract, 
	   &dmy[rec], rec, global_param->dt, options.Nnode, 
	   global_param->skipyear);

}
