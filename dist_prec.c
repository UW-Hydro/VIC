#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void dist_prec(atmos_data_struct   *atmos,
               dist_prcp_struct    *prcp,
               soil_con_struct      soil_con,
               veg_con_struct      *veg_con,
               dmy_struct          *dmy,
               global_param_struct  global_param,
               outfiles_struct      outfiles,
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

  extern option_struct options;
  extern debug_struct debug;

  static char STILL_STORM;
  static int  DRY_TIME;

  char    ANY_SNOW;
  int     veg, i;
  double  prec[(MAX_VEG+1)*2];
  double  rainonly[(MAX_VEG+1)*2];
  double  melt[(MAX_VEG+1)*2];
  double  NEW_MU;

  if(options.DIST_PRCP) {

    /*****************************************
      Controls Distributed Precipitation Model
      *****************************************/
 
    
    NEW_MU = 1.0 - exp(-options.PREC_EXPT*atmos->prec);
    for(veg=0;veg<=veg_con[0].vegetat_type_num;veg++) {
      ANY_SNOW = FALSE;
      for(i=0;i<options.SNOW_BAND;i++)
        /* Check for snow on ground or falling */
	if(prcp[0].snow[veg][i].snow 
	   || prcp[0].snow[veg][i].snow_canopy>0.) ANY_SNOW = TRUE;
      if(options.SNOW_MODEL && (ANY_SNOW || atmos->prec-atmos->rainonly>0)) {
        /* If snow present, mu must be set to 1. */
	NEW_MU = 1.;
	if(rec==0) {
          /* Set model variables if first time step */
	  prcp[0].mu[veg]=NEW_MU;
	  if(atmos->prec>0) STILL_STORM=TRUE;
	  else STILL_STORM=FALSE;
          DRY_TIME = 0;
	} 
	ANY_SNOW = TRUE;
      }
      else {
	if(rec==0) {
	  if(atmos->prec==0) {
	    /* If first time step has no rain, than set mu to 1. */
	    prcp[0].mu[veg] = 1.;
	    NEW_MU=1.;
	    STILL_STORM = TRUE;
	    DRY_TIME = 24;
	  }
	  else {
	    /* If first time step has rain, then set mu based on intensity */
	    prcp[0].mu[veg]=NEW_MU;
	    STILL_STORM=TRUE;
	    DRY_TIME = 0;
	  }
	}
	else if(atmos->prec==0 && DRY_TIME >= 24./(float)global_param.dt) {
          /* Check if storm has ended */
	  NEW_MU=prcp[0].mu[veg];
	  STILL_STORM=FALSE;
          DRY_TIME = 0;
	}
        else if(atmos->prec==0) {
	  /* May be pause in storm, keep track of pause length */
	  NEW_MU=prcp[0].mu[veg];
	  DRY_TIME += global_param.dt;
	}
      }

      if(!STILL_STORM && (atmos->prec>STORM_THRES || ANY_SNOW)) {
	/** Average soil moisture before a new storm **/
	initialize_new_storm(prcp[0].cell,prcp[0].veg_var,
			     veg,veg_con[0].vegetat_type_num,rec,
			     prcp[0].mu[veg],NEW_MU);
	STILL_STORM=TRUE;
	prcp[0].mu[veg] = NEW_MU;
      }
      else if(NEW_MU != prcp[0].mu[veg] && STILL_STORM) {
	/** Redistribute soil moisture during the storm if mu changes **/
	redistribute_during_storm(prcp[0].cell,prcp[0].veg_var,veg,
				  veg_con[0].vegetat_type_num,rec,
				  prcp[0].mu[veg],NEW_MU,
				  soil_con.max_moist);
	prcp[0].mu[veg] = NEW_MU;
      }
    }

    /** Prepare precipitation and melt forcing data for wet 
	and dry fractions **/
    for(i=0; i<(veg_con[0].vegetat_type_num+1); i++) {
      prec[i*2] = atmos->prec;
      prec[(i*2)+1] = 0;
      rainonly[i*2] = atmos->rainonly;
      rainonly[(i*2)+1] = 0;
      if(options.FULL_ENERGY || options.SNOW_MODEL) {
	prec[i*2] /= prcp[0].mu[i];
	rainonly[i*2] /= prcp[0].mu[i];
	melt[i*2] = 0.;
	melt[i*2+1] = 0.;
      }
      else {
	prec[i*2] = (atmos->prec - atmos->melt) / prcp[0].mu[i] + atmos->melt;
	melt[i*2] = atmos->melt;
	melt[i*2+1] = atmos->melt;
      }
    }

    /** Solve model time step **/
    full_energy(rec,atmos,soil_con,veg_con,prcp,dmy,
		global_param,cellnum,NEWCELL,prec,
		rainonly,melt);

  }

  else {

    /************************************************
      Controls Grid Cell Averaged Precipitation Model
      ************************************************/

    for(i=0; i<veg_con[0].vegetat_type_num+1; i++) {
      prec[i*2]     = atmos->prec;
      prec[i*2+1]     = 0.;
      rainonly[i*2] = atmos->rainonly;
      rainonly[i*2+1] = 0.;
      melt[i*2]     = atmos->melt;
      melt[i*2+1]     = 0;
    }
    full_energy(rec,atmos,soil_con,veg_con,prcp,dmy,
		global_param,cellnum,NEWCELL,prec,rainonly,melt);

  }

  put_data(prcp, atmos, veg_con, outfiles, soil_con.depth, 
	   soil_con.AreaFract, &dmy[rec], rec, global_param.dt,
	   global_param.skipyear);

}
