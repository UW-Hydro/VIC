#include <stdio.h>
#include <vicNl.h>
#include <math.h>

void dist_prec(atmos_data_struct *atmos,
               dist_prcp_struct *prcp,
               soil_con_struct soil_con,
               veg_con_struct *veg_con,
               dmy_struct *dmy,
               global_param_struct global_param,
               outfiles_struct outfiles,
               int rec,
               int cellnum,
               int *STILL_STORM,
               char NEWCELL,
               char LASTREC) {
/**********************************************************************
  dist_prec		Keith Cherkauer		October 9, 1997

  This subroutine calls the solution routines for a single grid cell
  for one time step.  It also controls the distribution of precipitation
  and reassembles grid cell data for output.

**********************************************************************/

  extern option_struct options;
  extern debug_struct debug;

  double prec;
  double rainonly;
  double melt;
  double NEW_STORM;

  if(options.DIST_PRCP) {

    /*****************************************
      Controls Distributed Precipitation Model
      *****************************************/
 
    NEW_STORM = 1.0 - exp(-0.6*atmos->prec);
    if(atmos->prec==0 && rec==0) {
      prcp[0].mu = 1.;
      NEW_STORM=1.;
    }
    else if(rec==0) prcp[0].mu=NEW_STORM;
    else if(atmos->prec==0) NEW_STORM=prcp[0].mu;

    if(NEW_STORM != prcp[0].mu) {
      initialize_new_storm(prcp,veg_con[0].vegetat_type_num,rec,
			   prcp[0].mu,NEW_STORM,soil_con,
			   global_param.Ulayer+global_param.Llayer+2);
      prcp[0].mu = NEW_STORM;
    }

    /** Solve Wet Fraction (mu) of Grid Cell **/
    prec = atmos->prec;
    rainonly = atmos->rainonly;
    melt = atmos->melt;
    if(options.FULL_ENERGY || options.SNOW_MODEL) {
      atmos->prec /= prcp[0].mu;
      atmos->rainonly /= prcp[0].mu;
    }
    else {
      atmos->prec = (atmos->prec - atmos->melt) / prcp[0].mu + melt;
      atmos->melt = melt;
    }
    full_energy(rec,atmos,soil_con,veg_con,&prcp[0].dist[0],dmy,
		global_param,prcp[0].mu,0,cellnum,NEWCELL);

    /** Solve Dry Fraction (1-mu) of Grid Cell **/
    atmos->rainonly = 0.;
    if(options.FULL_ENERGY || options.SNOW_MODEL) {
      atmos->prec = 0.;
      atmos->melt = 0.;
    }
    else {
      atmos->prec = melt;
      atmos->melt = melt;
    }
    full_energy(rec,atmos,soil_con,veg_con,&prcp[0].dist[1],dmy,
		global_param,1.0-prcp[0].mu,1,cellnum,NEWCELL);
    atmos->prec = prec;
    atmos->rainonly = rainonly;
    atmos->melt = melt;
  }

  else {

    /************************************************
      Controls Grid Cell Averaged Precipitation Model
      ************************************************/

    full_energy(rec,atmos,soil_con,veg_con,&prcp[0].dist[0],dmy,
		global_param,1.,0,cellnum,NEWCELL);

  }

  put_data(prcp, atmos, veg_con, outfiles, soil_con.depth, &dmy[rec], rec);

}
