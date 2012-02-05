#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  dist_prec(atmos_data_struct   *atmos,
               dist_prcp_struct    *prcp,
               soil_con_struct     *soil_con,
               veg_con_struct      *veg_con,
	       lake_con_struct     *lake_con,
               dmy_struct          *dmy,
               global_param_struct *global_param,
               filep_struct        *filep,
               out_data_file_struct *out_data_files,
               out_data_struct     *out_data,
               save_data_struct    *save_data,
               int                  rec,
               int                  cellnum,
               char                 NEWCELL,
               char                 LASTREC,
	       char                *init_STILL_STORM,
	       int                 *init_DRY_TIME) {
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
  03-05-01 Fixed error in which distributed precipitation accounting
           variables (DRY_TIME, STILL_STORM, ANY_SNOW) were used 
           within the vegetation loop, but did not store separate
           values for each vegetation type.                     KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  03-27-03 Modified calculation of DRY_TIME.  Originally the check
           to see if a new storm was warranted checked if DRY_TIME
           was greater than 24/dt.  However, DRY_TIME is incremented
           by dt, so it was checking hours against time steps.  The
           division by dt has been removed, so a new storm starts if
           the cell has been drying for a full 24 hours.     RS & KAC
  04-10-03 Modified to store STILL_STORM and DRY_TIME in the model
           statefile, so that full conditions will be preserved.  KAC
  01-Nov-04 Added support for state files containing SPATIAL_FROST and
	    LAKE_MODEL state variables.					TJB
  02-Feb-05 Modified to save state file at the end of the final timestep
	    of the date indicated by STATEYEAR, STATEMONTH, and STATEDAY
	    in the global parameter file.				TJB
  2005-Mar-24 Modified parameter list of put_data() to accomodate support
	      for ALMA variables.					TJB
  2006-Sep-23 Implemented flexible output configuration; uses new out_data,
	      out_data_files, and save_data structures.			TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
              main subroutine, rather than ending the simulation.	GCT/KAC
  2008-Oct-23 Modified call to put_data() to store ErrorFlag.		TJB
  2009-Mar-03 Modified routine to store put_data() error in ErrorFlag2 and 
	      return a single ERROR value if an error occurs.		KAC via TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Sep-28 Added logic for initial (pre-simulation) call to put_data.TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
**********************************************************************/

  extern option_struct   options;
  extern veg_lib_struct *veg_lib; 

  static char STILL_STORM[MAX_VEG];
  static int  DRY_TIME[MAX_VEG];

  char    ANY_SNOW[MAX_VEG];
  int     veg, i;
  int     month;
  int     ErrorFlag, ErrorFlag2;
  double  Wdmax;
  double  NEW_MU;

  /**************************************************
    If rec < 0, initialize the storage terms for water and energy balances
  **************************************************/
  if (rec < 0) {
    ErrorFlag2 = put_data(prcp, atmos, soil_con, veg_con,
			  lake_con, out_data_files, out_data, save_data,
			  &dmy[0], rec);
    if ( ErrorFlag2 == ERROR ) ErrorFlag = ERROR;
    return (0);
  }

  /**************************************************
    If rec >= 0, proceed with simulation
  **************************************************/
  // check if state file has been used to initialize storm tracking
  if ( init_DRY_TIME >= 0 ) {
    // initialize storm tracking
    for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
      DRY_TIME[veg] = init_DRY_TIME[veg];
      STILL_STORM[veg] = init_STILL_STORM[veg];
    }
  }

  if(options.DIST_PRCP) {

    /*******************************************
      Controls Distributed Precipitation Model
    *******************************************/
     
    NEW_MU = 1.0 - exp(-options.PREC_EXPT*atmos->prec[NR]);
    for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
      ANY_SNOW[veg] = FALSE;
      for ( i = 0; i < options.SNOW_BAND; i++ )
        /* Check for snow on ground or falling */
	if ( prcp->snow[veg][i].swq > 0 
	     || prcp->snow[veg][i].snow_canopy > 0. ) 
	  ANY_SNOW[veg] = TRUE;
      if ( ANY_SNOW[veg] || atmos->snowflag[NR] ) {
        /* If snow present, mu must be set to 1. */
	NEW_MU = 1.;
	if ( rec == 0 ) {
          /* Set model variables if first time step */
	  prcp->mu[veg] = NEW_MU;
	  if ( atmos->prec[NR] > 0 ) 
	    STILL_STORM[veg] = TRUE;
	  else 
	    STILL_STORM[veg] = FALSE;
          DRY_TIME[veg] = 0;
	} 
	ANY_SNOW[veg] = TRUE;
      }
      else {
	if ( rec == 0 ) {
	  if ( atmos->prec[NR] == 0 ) {
	    /* If first time step has no rain, than set mu to 1. */
	    prcp->mu[veg]    = 1.;
	    NEW_MU           = 1.;
	    STILL_STORM[veg] = TRUE;
	    DRY_TIME[veg]    = 24;
	  }
	  else {
	    /* If first time step has rain, then set mu based on intensity */
	    prcp->mu[veg]    = NEW_MU;
	    STILL_STORM[veg] = TRUE;
	    DRY_TIME[veg]    = 0;
	  }
	}
	else if(atmos->prec[NR] == 0 && DRY_TIME[veg] >= 24.) {
          /* Check if storm has ended */
	  NEW_MU           = prcp->mu[veg];
	  STILL_STORM[veg] = FALSE;
          DRY_TIME[veg]    = 0;
	}
        else if ( atmos->prec[NR] == 0 ) {
	  /* May be pause in storm, keep track of pause length */
	  NEW_MU         = prcp->mu[veg];
	  DRY_TIME[veg] += global_param->dt;
	}
      }

      if ( !STILL_STORM[veg] && (atmos->prec[NR] > STORM_THRES 
				 || ANY_SNOW[veg] ) ) {
	/** Average soil moisture before a new storm **/
	ErrorFlag = initialize_new_storm(prcp->cell,prcp->veg_var,
			     veg,veg_con[0].vegetat_type_num,rec,
			     prcp->mu[veg],NEW_MU);
        if ( ErrorFlag == ERROR ) return ( ERROR );

	STILL_STORM[veg] = TRUE;
	prcp->mu[veg]    = NEW_MU;
      }
      else if ( NEW_MU != prcp->mu[veg] && STILL_STORM[veg] ) {
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
    ErrorFlag = full_energy(NEWCELL, cellnum, rec, atmos, prcp, dmy, global_param, 
		lake_con, soil_con, veg_con);

  }

  else {

    /**************************************************
      Controls Grid Cell Averaged Precipitation Model
    **************************************************/

    ErrorFlag = full_energy(NEWCELL, cellnum, rec, atmos, prcp, dmy, global_param, 
		lake_con, soil_con, veg_con);

  }

  /**************************************************
    Write cell average values for current time step
  **************************************************/

  ErrorFlag2 = put_data(prcp, atmos, soil_con, veg_con,
			lake_con, out_data_files, out_data, save_data,
			&dmy[rec], rec);
  if ( ErrorFlag2 == ERROR ) ErrorFlag = ERROR;

  /************************************
    Save model state at assigned date
    (after the final time step of the assigned date)
  ************************************/

  if ( filep->statefile != NULL
       &&  ( dmy[rec].year == global_param->stateyear
	     && dmy[rec].month == global_param->statemonth 
	     && dmy[rec].day == global_param->stateday
	     && ( rec+1 == global_param->nrecs
		  || dmy[rec+1].day != global_param->stateday ) ) )
    write_model_state(prcp, global_param, veg_con[0].vegetat_type_num, 
		      soil_con->gridcel, filep, soil_con,
		      STILL_STORM, DRY_TIME, *lake_con);

  return ( ErrorFlag );

}
