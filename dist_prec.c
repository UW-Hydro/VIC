#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

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

**********************************************************************/

  extern option_struct options;
  extern debug_struct debug;

  static char STILL_STORM;

  char    ANY_SNOW;
  int     veg, i;
  double *prec;
  double *rainonly;
  double *melt;
  double  NEW_MU;

  prec     = (double *)calloc((veg_con[0].vegetat_type_num+1)*2,
			      sizeof(double));
  rainonly = (double *)calloc((veg_con[0].vegetat_type_num+1)*2,
			      sizeof(double));
  melt     = (double *)calloc((veg_con[0].vegetat_type_num+1)*2,
			      sizeof(double));

  if(options.DIST_PRCP) {

    /*****************************************
      Controls Distributed Precipitation Model
      *****************************************/
 
    NEW_MU = 1.0 - exp(-0.6*atmos->prec);
    for(veg=0;veg<=veg_con[0].vegetat_type_num;veg++) {
      ANY_SNOW = FALSE;
      for(i=0;i<options.SNOW_BAND;i++)
	if(prcp[0].snow[veg][i].snow 
	   || prcp[0].snow[veg][i].snow_canopy>0.) ANY_SNOW = TRUE;
      if(options.SNOW_MODEL && (ANY_SNOW || atmos->prec-atmos->rainonly>0)) {
	NEW_MU = 1.;
	if(rec==0) {
	  prcp[0].mu[veg]=NEW_MU;
	  if(atmos->prec>0) STILL_STORM=TRUE;
	  else STILL_STORM=FALSE;
	} 
      }
      else {
	if(atmos->prec==0 && rec==0) {
	  prcp[0].mu[veg] = 1.;
	  NEW_MU=1.;
	  STILL_STORM = TRUE;
	}
	else if(rec==0) {
	  prcp[0].mu[veg]=NEW_MU;
	  STILL_STORM=TRUE;
	}
	else if(atmos->prec==0) {
	  NEW_MU=prcp[0].mu[veg];
	  STILL_STORM=FALSE;
	}
      }

      if(!STILL_STORM && atmos->prec>STORM_THRES) {
	initialize_new_storm(prcp[0].cell,prcp[0].veg_var,
			     veg,veg_con[0].vegetat_type_num,rec,
			     prcp[0].mu[veg],NEW_MU);
	STILL_STORM=TRUE;
	prcp[0].mu[veg] = NEW_MU;
      }
      else if(NEW_MU != prcp[0].mu[veg] && STILL_STORM) {
	redistribute_during_storm(prcp[0].cell,prcp[0].veg_var,veg,
				  veg_con[0].vegetat_type_num,rec,
				  prcp[0].mu[veg],NEW_MU);
	prcp[0].mu[veg] = NEW_MU;
      }
    }

    /** Solve Wet Fraction (mu) of Grid Cell **/
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
    full_energy(rec,atmos,soil_con,veg_con,prcp,dmy,
		global_param,cellnum,NEWCELL,prec,
		rainonly,melt);

/***
    atmos->prec = prec;
    atmos->rainonly = rainonly;
    atmos->melt = melt;
***/
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

  free((char *)prec);
  free((char *)rainonly);
  free((char *)melt);

  put_data(prcp, atmos, veg_con, outfiles, soil_con.depth, soil_con.AreaFract,
	   &dmy[rec], rec);

}
