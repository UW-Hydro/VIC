#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

int initialize_model_state(all_vars_struct     *all_vars,
			   all_vars_struct     *all_vars_crop,
			   dmy_struct           dmy,
			   global_param_struct *global_param,
			   filep_struct         filep,
			   int                  cellnum,
			   int                  Nveg,
			   int                  Nnodes,
			   double               surf_temp, 
			   soil_con_struct     *soil_con,
			   veg_con_struct      *veg_con,
			   lake_con_struct      lake_con)
/**********************************************************************
  initialize_model_state      Keith Cherkauer	    April 17, 2000

  This routine initializes the model state (energy balance, water balance,
  and snow components).  If a state file is provided to the model than its
  contents are checked to see if it agrees with the current simulation
  set-up, if so it is used to initialize the model state.  If no state
  file is provided the model initializes all variables with defaults and
  the user should expect to throw out the beginning of the simulation 
  period as model start-up.

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  Modifications:
  4-17-00 Modified from initialize_energy_bal.c and initialize_snow.c
          to provide a single controlling routine for initializing the
          model state.
  9-00    Fixed bug where initialization of soil node temperatures 
          and moitures was within two vegetation loops, thus only
          the first vegetation type was properly initialized.     KAC
  2-19-03 Modified to initialize soil and vegetation parameters for
          the dry grid cell fraction, if distributed precipitation
          is activated.                                           KAC
  11-18-02 Modified to initialize lake and wetland algorithms 
          variables.                                              LCB
  2-10-03 Fixed looping problem with initialization of soil moisture. KAC
  3-12-03 Modified so that soil layer ice content is only calculated 
          when frozen soil is implemented and active in the current 
          grid cell.                                                KAC
  04-10-03 Modified to read storm parameters from model state file.  KAC
  04-25-03 Modified to work with vegetation type specific storm 
           parameters.                                              KAC
  07-May-04 Initialize soil_con->dz_node[Nnodes] to 0.0, since it is
	    accessed in set_node_parameters().				TJB
  01-Nov-04 Added support for state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.				TJB
  2006-Apr-21 Replaced Cv (uninitialized) with lake_con.Cl[0] in
	      surfstor calculation.					TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              save_data structure to track changes in moisture storage
              over each time step; this needs initialization here.	TJB
  2006-Oct-10 Added snow[veg][band].snow_canopy to save_data.swe.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included removing the unused init_snow file.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-24 Added EXP_TRANS option.					JCA
  2007-Apr-24 Zsum_node loaded into soil_con structure for later use
	      without having to recalculate.				JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2007-Aug-21 Return value of ErrorFlag if error in
	      distribute_node_moisture_properties.			JCA
  2007-Sep-18 Check for soil moist exceeding max moist moved from
	      read_initial_model_state to here.				JCA
  2007-Oct-24 Modified initialize_lake() to return ErrorFlag.		TJB
  2008-Mar-01 Reinserted missing logic for QUICK_FS in calls to
	      distribute_node_moisture_properties() and
	      estimate_layer_ice_content().				TJB
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
  2009-Mar-15 Modified to not call estimate_layer_ice_content() if
	      not modeling frozen soil.					KAC via TJB
  2009-Mar-16 Added resid_moist to argument list of
	      estimate_layer_ice_content().  This allows computation
	      of min_liq, the minimum allowable liquid water content
	      in each layer as a function of temperature.		TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-26 Added initial estimate of incoming longwave at surface
	      (LongUnderOut) for use in canopy snow T iteration.	TJB
  2009-Jul-31 Removed extra lake/wetland veg tile.			TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-19 Made initialization of Tfoliage more accurate for snow bands.	TJB
  2009-Sep-28 Added initialization of energy structure.			TJB
  2009-Nov-15 Added check to ensure that depth of first thermal node
	      is <= depth of first soil layer.				TJB
  2009-Dec-11 Removed initialization of save_data structure, since this
	      is now performed by the initial call to put_data().	TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.			TJB
  2010-Nov-11 Updated call to initialize_lake() to accommodate new
	      skip_hydro flag.						TJB
  2011-Mar-01 Updated calls to initialize_soil() and initialize_lake()
	      to accommodate new arguments.  Added more detailed validation
	      of soil moisture.						TJB
  2011-Mar-05 Added validation of initial soil moisture, ice, and snow
	      variables to make sure they are self-consistent.		TJB
  2011-May-31 Removed options.GRND_FLUX.  Now soil temperatures and 
	      ground flux are always computed.				TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Jul-05 Changed logic initializing soil temperatures so that
	      type of initialization depends solely on options.QUICK_FLUX;
	      options.Nnodes is no longer automatically reset here.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Jan-28 Added stability check for case of (FROZEN_SOIL=TRUE &&
	      IMPLICIT=FALSE).						TJB
  2013-Jul-25 Fixed incorrect condition on lake initialization.		TJB
  2013-Jul-25 Moved computation of tmp_moist argument of
	      compute_runoff_and_asat() so that it would always be
	      initialized.						TJB
  2013-Dec-26 Removed EXCESS_ICE option.						TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-27 Removed QUICK_FS option.							TJB
  2014-Jan-13 Added validation of Nnodes and dp for EXP_TRANS=TRUE.	TJB
  2014-Feb-09 Made non-spinup initial temperatures more consistent with
	      annual average air temperature and bottom boundary
	      temperature.											TJB
  2014-Mar-28 Removed DIST_PRCP option.							TJB
**********************************************************************/
{
  extern option_struct options;
  extern veg_lib_struct *veg_lib;

  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      i, j, ii, veg, index;
  int      idx;
  int      cidx;
  int      lidx;
  int      nidx;
  double   tmp_moist[MAX_LAYERS];
  double   tmp_runoff;
  int      band;
  int      frost_area;
  int      ErrorFlag;
  double   Cv;
  double   Zsum, dp;
  double   tmpdp, tmpadj, Bexp;
  double   Tair;
  double   tmp;
  double  *M;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS][MAX_FROST_AREAS];
  double   Clake;
  double   surf_swq;
  double   pack_swq;
  double   TreeAdjustFactor[MAX_BANDS];
  double dt_thresh;
  int tmp_lake_idx;

  cell_data_struct      **cell;
  energy_bal_struct     **energy;
  lake_var_struct        *lake_var;
  snow_data_struct      **snow;
  veg_var_struct        **veg_var;

  cell    = all_vars->cell;
  energy  = all_vars->energy;
  lake_var = &all_vars->lake_var;
  snow    = all_vars->snow;
  veg_var = all_vars->veg_var;

  // Initialize soil depths
  dp = soil_con->dp;

  FIRST_VEG = TRUE;

  // increase initial soil surface temperature if air is very cold
  Tair = surf_temp;
  if ( surf_temp < -1. ) surf_temp = -1.;
  
  /********************************************
    Initialize all snow pack variables 
    - some may be reset if state file present
  ********************************************/

  initialize_snow(snow, Nveg, cellnum);

  /********************************************
    Initialize all soil layer variables 
    - some may be reset if state file present
  ********************************************/

  initialize_soil(cell, soil_con, veg_con, Nveg);

  /********************************************
    Initialize all vegetation variables 
    - some may be reset if state file present
  ********************************************/

  initialize_veg(veg_var, veg_con, global_param, Nveg);

  /********************************************
    Initialize all lake variables 
  ********************************************/

  if ( options.LAKES ) {
    tmp_lake_idx = lake_con.lake_idx;
    if (tmp_lake_idx < 0) tmp_lake_idx = 0;
    ErrorFlag = initialize_lake(lake_var, lake_con, soil_con, &(cell[tmp_lake_idx][0]), surf_temp, 0);
    if (ErrorFlag == ERROR) return(ErrorFlag);
  }

  /********************************************
    Initialize all spatial frost variables 
  ********************************************/

  for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
    if ( options.Nfrost == 1 ) soil_con->frost_fract[frost_area] = 1.;
    else if (options.Nfrost == 2 ) soil_con->frost_fract[frost_area] = 0.5;
    else {
      soil_con->frost_fract[frost_area] = 1. / (options.Nfrost - 1);
      if ( frost_area == 0 || frost_area == options.Nfrost-1 ) 
	soil_con->frost_fract[frost_area] /= 2.;
    }
  }

  /********************************************************
    Compute grid cell fractions for all subareas used in 
    spatial distribution of soil frost routines.
  ********************************************************/

  /************************************************************************
    CASE 1: Not using quick ground heat flux, and initial conditions files 
    provided
  ************************************************************************/

  if(options.INIT_STATE) {

    read_initial_model_state(filep.init_state, all_vars, global_param,  
			     Nveg, options.SNOW_BAND, cellnum, soil_con,
			     lake_con);

    /******Check that soil moisture does not exceed maximum allowed************/
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      for( band = 0; band < options.SNOW_BAND; band++ ) {
        for( lidx = 0; lidx < options.Nlayer; lidx++ ) {	  

	  if ( cell[veg][band].layer[lidx].moist > soil_con->max_moist[lidx] ) {
            fprintf( stderr, "WARNING: Initial soil moisture (%f mm) exceeds maximum (%f mm) in layer %d for veg tile %d and snow band%d.  Resetting to maximum.\n", cell[veg][band].layer[lidx].moist, soil_con->max_moist[lidx], lidx, veg, band );
            for ( frost_area = 0; frost_area < options.Nfrost; frost_area++)
              cell[veg][band].layer[lidx].ice[frost_area] *= soil_con->max_moist[lidx]/cell[veg][band].layer[lidx].moist;
            cell[veg][band].layer[lidx].moist = soil_con->max_moist[lidx];
	  }

          for ( frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            if (cell[veg][band].layer[lidx].ice[frost_area] > cell[veg][band].layer[lidx].moist)
              cell[veg][band].layer[lidx].ice[frost_area] = cell[veg][band].layer[lidx].moist;
          }
          tmp_moist[lidx] = cell[veg][band].layer[lidx].moist;

	}
        compute_runoff_and_asat(soil_con, tmp_moist, 0, &(cell[veg][band].asat), &tmp_runoff);
      }

      // Override possible bad values of soil moisture under lake coming from state file
      // (ideally we wouldn't store these in the state file in the first place)
      if (options.LAKES && veg == lake_con.lake_idx) {
        for( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          lake_var->soil.layer[lidx].moist = soil_con->max_moist[lidx];
          for ( frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            if (lake_var->soil.layer[lidx].ice[frost_area] > lake_var->soil.layer[lidx].moist)
              lake_var->soil.layer[lidx].ice[frost_area] = lake_var->soil.layer[lidx].moist;
          }
        }
      }
    }


    /****** initialize moist and ice ************/
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  for( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[veg][band].layer[lidx].moist;

	    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
	      ice[veg][band][lidx][frost_area] = cell[veg][band].layer[lidx].ice[frost_area];
	  }
	}
      }
    }

    /******Check that snow pack terms are self-consistent************/
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      for ( band = 0 ; band < options.SNOW_BAND ; band++ ) {
        if (snow[veg][band].swq > MAX_SURFACE_SWE) {
          pack_swq = snow[veg][band].swq-MAX_SURFACE_SWE;
          surf_swq = MAX_SURFACE_SWE;
        }
        else {
          pack_swq = 0;
          surf_swq = snow[veg][band].swq;
          snow[veg][band].pack_temp = 0;
        }
        if (snow[veg][band].surf_water > LIQUID_WATER_CAPACITY*surf_swq) {
          snow[veg][band].pack_water += snow[veg][band].surf_water - (LIQUID_WATER_CAPACITY*surf_swq);
          snow[veg][band].surf_water = LIQUID_WATER_CAPACITY*surf_swq;
        }
        if (snow[veg][band].pack_water > LIQUID_WATER_CAPACITY*pack_swq) {
          snow[veg][band].pack_water = LIQUID_WATER_CAPACITY*pack_swq;
        }
      }
    }

  }
  
  /************************************************************************
    CASE 2: Initialize soil if using quick heat flux, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(options.QUICK_FLUX) {
    Nnodes = options.Nnode;

    /* Initialize soil node thicknesses */
    soil_con->dz_node[0]   = soil_con->depth[0];
    soil_con->dz_node[1]   = soil_con->depth[0];
    soil_con->dz_node[2]   = 2. * (dp - 1.5 * soil_con->depth[0]);    
    soil_con->Zsum_node[0] = 0;
    soil_con->Zsum_node[1] = soil_con->depth[0];
    soil_con->Zsum_node[2] = dp;

    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {

	  /* Initialize soil node temperatures */
	  energy[veg][band].T[0] = surf_temp;
	  energy[veg][band].T[1] = surf_temp;
	  energy[veg][band].T[2] = soil_con->avg_temp;

	  /* Initialize soil layer moisture and ice contents */
	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[veg][band].layer[lidx].moist;
	    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
	  }
	}
      }
    }
  }

  /*****************************************************************
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
  *****************************************************************/
  else if(!options.QUICK_FLUX) {
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;
      
      if ( Cv > 0 ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  
	  if(!options.EXP_TRANS){  
	    /* Initialize soil node temperatures and thicknesses 
	       Nodes set at surface, the depth of the first layer,
	       twice the depth of the first layer, and at the
	       damping depth.  Extra nodes are placed equal distance
	       between the damping depth and twice the depth of the
	       first layer. */
	    
	    soil_con->dz_node[0] = soil_con->depth[0];
	    soil_con->dz_node[1] = soil_con->depth[0];
	    soil_con->dz_node[2] = soil_con->depth[0];

	    soil_con->Zsum_node[0] = 0;
	    soil_con->Zsum_node[1] = soil_con[0].depth[0];
	    Zsum   = 2. * soil_con[0].depth[0];
	    soil_con->Zsum_node[2] = Zsum;
	    tmpdp  = dp - soil_con[0].depth[0] * 2.5;
	    tmpadj = 3.5;
	    for ( index = 3; index < Nnodes-1; index++ ) {
	      if ( FIRST_VEG ) {
		soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
	      }
	      Zsum += (soil_con->dz_node[index]
		       +soil_con->dz_node[index-1])/2.;
	      soil_con->Zsum_node[index] = Zsum;
	    }
	    energy[veg][band].T[0] = surf_temp;
	    for ( index = 1; index < Nnodes; index++ ) {
	      energy[veg][band].T[index] = soil_con->avg_temp;
	    }
	    if ( FIRST_VEG ) {
	      FIRST_VEG = FALSE;
	      soil_con->dz_node[Nnodes-1] = (dp - Zsum 
					     - soil_con->dz_node[Nnodes-2] 
					     / 2. ) * 2.;
	      Zsum += (soil_con->dz_node[Nnodes-2]
		       +soil_con->dz_node[Nnodes-1])/2.;
	      soil_con->Zsum_node[Nnodes-1] = Zsum;
	      if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
		sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
		nrerror(ErrStr);
	      }
	    }
	  }
	  else{ /* exponential grid transformation, EXP_TRANS = TRUE*/
	    
	    /*calculate exponential function parameter */
	    if ( FIRST_VEG ) {
	      Bexp = logf(dp+1.)/(double)(Nnodes-1); //to force Zsum=dp at bottom node
              /* validate Nnodes by requiring that there be at least 3 nodes in the top 50cm */
              if (Nnodes < 5*logf(dp+1.)+1) {
		sprintf(ErrStr,"The number of soil thermal nodes (%d) is too small for the supplied damping depth (%f) with EXP_TRANS set to TRUE, leading to fewer than 3 nodes in the top 50 cm of the soil column.  For EXP_TRANS=TRUE, Nnodes and dp must follow the relationship:\n5*ln(dp+1)<Nnodes-1\nEither set Nnodes to at least %d in the global param file or reduce damping depth to %f in the soil parameter file.  Or set EXP_TRANS to FALSE in the global parameter file.",Nnodes,dp,(int)(5*logf(dp+1.))+2,exp(0.2*(Nnodes-1))+1);
                nrerror(ErrStr);
              }
	      for ( index = 0; index <= Nnodes-1; index++ )
		soil_con->Zsum_node[index] = expf(Bexp*index)-1.;
	      if(soil_con->Zsum_node[0] > soil_con->depth[0]) {
		sprintf(ErrStr,"Depth of first thermal node (%f) in initialize_model_state is greater than depth of first soil layer (%f); increase the number of nodes or decrease the thermal damping depth dp (%f)",soil_con->Zsum_node[0],soil_con->depth[0],dp);
		nrerror(ErrStr);
	      }
	    }	    
	    
	    //top node	  
	    index=0;
	    if ( FIRST_VEG )
	      soil_con->dz_node[index] = soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index];
	    energy[veg][band].T[index] = surf_temp;
	    //middle nodes
	    for ( index = 1; index < Nnodes-1; index++ ) {
	      if ( FIRST_VEG ) {
		soil_con->dz_node[index] = (soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index])/2.+(soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1])/2.;
	      }
//	      energy[veg][band].T[index] = exp_interp(soil_con->Zsum_node[index],0.,soil_con[0].dp,
//						      surf_temp,soil_con[0].avg_temp);
              energy[veg][band].T[index] = soil_con->avg_temp;
	    }
	    //bottom node
	    index=Nnodes-1;
	    if ( FIRST_VEG )
	      soil_con->dz_node[index] = soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1];
	    energy[veg][band].T[index] = soil_con->avg_temp;

	  } // end if !EXP_TRANS
	  
	  //initialize moisture and ice for each soil layer
	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[veg][band].layer[lidx].moist;
	    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
	  }
	}
      }
    }
  }

  /*********************************
    CASE 4: Unknown option
  *********************************/
  else {
    for ( veg = 0 ; veg <= Nveg ; veg++ ) {
      // Initialize soil for existing vegetation types
      Cv = veg_con[veg].Cv;

      if ( Cv > 0 ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  // Initialize soil for existing snow elevation bands
	  if ( soil_con->AreaFract[band] > 0. ) {	  
	    for ( index = 0; index < options.Nlayer; index++ ) {
	      soil_con->dz_node[index] = 1.;
	    }
	  }
	}
      }
    }
  }

  /******************************************
    Initialize soil thermal node properties 
  ******************************************/

  FIRST_VEG = TRUE;
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    // Initialize soil for existing vegetation types
    Cv = veg_con[veg].Cv;

    if ( Cv > 0 ) {
      for( band = 0; band < options.SNOW_BAND; band++ ) {
	// Initialize soil for existing snow elevation bands
	if ( soil_con->AreaFract[band] > 0. ) {
	    
	  /** Set soil properties for all soil nodes **/
	  if(FIRST_VEG) {
	    FIRST_VEG = FALSE;
	    set_node_parameters(soil_con->dz_node, soil_con->Zsum_node, soil_con->max_moist_node,
				soil_con->expt_node, soil_con->bubble_node,
				soil_con->alpha, soil_con->beta,
				soil_con->gamma, soil_con->depth,
				soil_con->max_moist, soil_con->expt, 
				soil_con->bubble, soil_con->quartz, 
				Nnodes, options.Nlayer, soil_con->FS_ACTIVE);	  
	  }
	
	  /* set soil moisture properties for all soil thermal nodes */
	  ErrorFlag = distribute_node_moisture_properties(energy[veg][band].moist,
						energy[veg][band].ice,
						energy[veg][band].kappa_node,
						energy[veg][band].Cs_node,
						soil_con->Zsum_node,
						energy[veg][band].T,
						soil_con->max_moist_node,
						soil_con->expt_node,
						soil_con->bubble_node,
						moist[veg][band], 
						soil_con->depth,
						soil_con->soil_dens_min,
						soil_con->bulk_dens_min,
						soil_con->quartz,
						soil_con->soil_density,
						soil_con->bulk_density,
						soil_con->organic,
						Nnodes, options.Nlayer,
						soil_con->FS_ACTIVE);
	  if ( ErrorFlag == ERROR ) return ( ErrorFlag );

          /* Check node spacing v time step */
          /* (note this is only approximate since heat capacity and conductivity can change considerably during the simulation depending on soil moisture and ice content) */
          if ((options.FROZEN_SOIL && !options.QUICK_FLUX) && !options.IMPLICIT) {
            dt_thresh = 0.5*energy[veg][band].Cs_node[1]/energy[veg][band].kappa_node[1]*pow((soil_con->dz_node[1]),2)/3600; // in hours
            if (global_param->dt > dt_thresh) {
              sprintf(ErrStr,"ERROR: You are currently running FROZEN SOIL with an explicit method (IMPLICIT is set to FALSE).  For the explicit method to be stable, time step %d hours is too large for the given thermal node spacing %f m, soil heat capacity %f J/m3/K, and soil thermal conductivity %f J/m/s/K.  Either set IMPLICIT to TRUE in your global parameter file (this is the recommended action), or decrease time step length to <= %f hours, or decrease the number of soil thermal nodes.",global_param->dt,soil_con->dz_node[1],energy[veg][band].Cs_node[1],energy[veg][band].kappa_node[1],dt_thresh);
              nrerror(ErrStr);
            }
          }

	  /* initialize layer moistures and ice contents */
	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    cell[veg][band].layer[lidx].moist = moist[veg][band][lidx];
	    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
              cell[veg][band].layer[lidx].ice[frost_area] = ice[veg][band][lidx][frost_area];
	  }
          if (options.QUICK_FLUX) {
            ErrorFlag = estimate_layer_ice_content_quick_flux(cell[veg][band].layer,
					 soil_con->depth, soil_con->dp,
					 energy[veg][band].T[0], energy[veg][band].T[1],
					 soil_con->avg_temp, soil_con->max_moist, 
					 soil_con->expt, soil_con->bubble, 
					 soil_con->frost_fract, soil_con->frost_slope, 
					 soil_con->FS_ACTIVE);
          }
          else {
	    ErrorFlag = estimate_layer_ice_content(cell[veg][band].layer,
						     soil_con->Zsum_node,
						     energy[veg][band].T,
						     soil_con->max_moist_node,
						     soil_con->expt_node,
						     soil_con->bubble_node,
						     soil_con->depth,
						     soil_con->max_moist,
						     soil_con->expt,
						     soil_con->bubble,
						     soil_con->frost_fract, 
						     soil_con->frost_slope, 
						     Nnodes, options.Nlayer, 
						     soil_con->FS_ACTIVE);
          }

	  /* Find freezing and thawing front depths */
	  if(!options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	    find_0_degree_fronts(&energy[veg][band], soil_con->Zsum_node, energy[veg][band].T, Nnodes);
	}
      }
    }
  }	

  // initialize miscellaneous energy balance terms
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    for ( band = 0; band < options.SNOW_BAND; band++ ) {
      /* Set fluxes to 0 */
      energy[veg][band].advected_sensible = 0.0;
      energy[veg][band].advection         = 0.0;
      energy[veg][band].AtmosError        = 0.0;
      energy[veg][band].AtmosLatent       = 0.0;
      energy[veg][band].AtmosLatentSub    = 0.0;
      energy[veg][band].AtmosSensible     = 0.0;
      energy[veg][band].canopy_advection  = 0.0;
      energy[veg][band].canopy_latent     = 0.0;
      energy[veg][band].canopy_latent_sub = 0.0;
      energy[veg][band].canopy_refreeze   = 0.0;
      energy[veg][band].canopy_sensible   = 0.0;
      energy[veg][band].deltaCC           = 0.0;
      energy[veg][band].deltaH            = 0.0;
      energy[veg][band].error             = 0.0;
      energy[veg][band].fusion            = 0.0;
      energy[veg][band].grnd_flux         = 0.0;
      energy[veg][band].latent            = 0.0;
      energy[veg][band].latent_sub        = 0.0;
      energy[veg][band].longwave          = 0.0;
      energy[veg][band].LongOverIn        = 0.0;
      energy[veg][band].LongUnderIn       = 0.0;
      energy[veg][band].LongUnderOut      = 0.0;
      energy[veg][band].melt_energy       = 0.0;
      energy[veg][band].NetLongAtmos      = 0.0;
      energy[veg][band].NetLongOver       = 0.0;
      energy[veg][band].NetLongUnder      = 0.0;
      energy[veg][band].NetShortAtmos     = 0.0;
      energy[veg][band].NetShortGrnd      = 0.0;
      energy[veg][band].NetShortOver      = 0.0;
      energy[veg][band].NetShortUnder     = 0.0;
      energy[veg][band].out_long_canopy   = 0.0;
      energy[veg][band].out_long_surface  = 0.0;
      energy[veg][band].refreeze_energy   = 0.0;
      energy[veg][band].sensible          = 0.0;
      energy[veg][band].shortwave         = 0.0;
      energy[veg][band].ShortOverIn       = 0.0;
      energy[veg][band].ShortUnderIn      = 0.0;
      energy[veg][band].snow_flux         = 0.0;
      /* Initial estimate of LongUnderOut for use by snow_intercept() */
      tmp = energy[veg][band].T[0] + KELVIN;
      energy[veg][band].LongUnderOut = STEFAN_B * tmp * tmp * tmp * tmp;
      energy[veg][band].Tfoliage     = Tair + soil_con->Tfactor[band];
    }
  }

  // initialize Tfallback counters
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    for ( band = 0; band < options.SNOW_BAND; band++ ) {
      energy[veg][band].Tfoliage_fbcount = 0;
      energy[veg][band].Tcanopy_fbcount = 0;
      energy[veg][band].Tsurf_fbcount = 0;
      for ( index = 0; index < Nnodes-1; index++ ) {
	energy[veg][band].T_fbcount[index] = 0;
      }
    }
  }

  // Compute treeline adjustment factors
  for ( band = 0; band < options.SNOW_BAND; band++ ) {
    if ( soil_con->AboveTreeLine[band] ) {
      Cv = 0;
      for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++ ) {
        if ( veg_lib[veg_con[veg].veg_class].overstory )
          Cv += veg_con[veg].Cv;
      }
      TreeAdjustFactor[band] = 1. / ( 1. - Cv );
    }
    else TreeAdjustFactor[band] = 1.;
  }

  // Initialize crop structures
  for ( veg = 0 ; veg < veg_con[0].vegetat_type_num ; veg++ ) {
    if (veg_con[veg].crop_frac_active) {
      for ( band = 0; band < options.SNOW_BAND; band++ ) {
      for (idx = veg_con[veg].crop_frac_idx; idx<veg_con[veg].crop_frac_idx+2; idx++) {

        // Copy veg_var state data
        if (idx % 2 == 0)
          all_vars_crop->veg_var[idx][band].Wdew = 0;
        else
          all_vars_crop->veg_var[idx][band].Wdew = veg_var[veg][band].Wdew;
        if (options.CARBON) {
          for (cidx=0; cidx<options.Ncanopy; cidx++) {
            all_vars_crop->veg_var[idx][band].NscaleFactor[cidx] = veg_var[veg][band].NscaleFactor[cidx];
            all_vars_crop->veg_var[idx][band].aPARLayer[cidx] = veg_var[veg][band].aPARLayer[cidx];
            all_vars_crop->veg_var[idx][band].CiLayer[cidx] = veg_var[veg][band].CiLayer[cidx];
            all_vars_crop->veg_var[idx][band].rsLayer[cidx] = veg_var[veg][band].rsLayer[cidx];
          }
        }
        all_vars_crop->veg_var[idx][band].Ci = veg_var[veg][band].Ci;
        all_vars_crop->veg_var[idx][band].rc = veg_var[veg][band].rc;
        all_vars_crop->veg_var[idx][band].NPPfactor = veg_var[veg][band].NPPfactor;
        all_vars_crop->veg_var[idx][band].AnnualNPP = veg_var[veg][band].AnnualNPP;
        all_vars_crop->veg_var[idx][band].AnnualNPPPrev = veg_var[veg][band].AnnualNPPPrev;

        // Copy veg_var flux data
        all_vars_crop->veg_var[idx][band].canopyevap = veg_var[veg][band].canopyevap;
        all_vars_crop->veg_var[idx][band].throughfall = veg_var[veg][band].throughfall;
        all_vars_crop->veg_var[idx][band].aPAR = veg_var[veg][band].aPAR;
        all_vars_crop->veg_var[idx][band].GPP = veg_var[veg][band].GPP;
        all_vars_crop->veg_var[idx][band].Rphoto = veg_var[veg][band].Rphoto;
        all_vars_crop->veg_var[idx][band].Rdark = veg_var[veg][band].Rdark;
        all_vars_crop->veg_var[idx][band].Rmaint = veg_var[veg][band].Rmaint;
        all_vars_crop->veg_var[idx][band].Rgrowth = veg_var[veg][band].Rgrowth;
        all_vars_crop->veg_var[idx][band].Raut = veg_var[veg][band].Raut;
        all_vars_crop->veg_var[idx][band].NPP = veg_var[veg][band].NPP;
        all_vars_crop->veg_var[idx][band].Litterfall = veg_var[veg][band].Litterfall;

        // Copy cell state data
	for ( lidx = 0; lidx < 2; lidx++ ) {
          all_vars_crop->cell[idx][band].aero_resist[lidx] = cell[veg][band].aero_resist[lidx];
        }
        all_vars_crop->cell[idx][band].asat = cell[veg][band].asat;
        all_vars_crop->cell[idx][band].CLitter = cell[veg][band].CLitter;
        all_vars_crop->cell[idx][band].CInter = cell[veg][band].CInter;
        all_vars_crop->cell[idx][band].CSlow = cell[veg][band].CSlow;
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          all_vars_crop->cell[idx][band].layer[lidx].bare_evap_frac = cell[veg][band].layer[lidx].bare_evap_frac;
          all_vars_crop->cell[idx][band].layer[lidx].Cs = cell[veg][band].layer[lidx].Cs;
          all_vars_crop->cell[idx][band].layer[lidx].kappa = cell[veg][band].layer[lidx].kappa;
          all_vars_crop->cell[idx][band].layer[lidx].moist = cell[veg][band].layer[lidx].moist;
	  for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
            all_vars_crop->cell[idx][band].layer[lidx].ice[frost_area] = cell[veg][band].layer[lidx].ice[frost_area];
          all_vars_crop->cell[idx][band].layer[lidx].phi = cell[veg][band].layer[lidx].phi;
          all_vars_crop->cell[idx][band].layer[lidx].T = cell[veg][band].layer[lidx].T;
          all_vars_crop->cell[idx][band].layer[lidx].zwt = cell[veg][band].layer[lidx].zwt;
        }

        // Copy cell flux data
        all_vars_crop->cell[idx][band].baseflow = cell[veg][band].baseflow;
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
          all_vars_crop->cell[idx][band].layer[lidx].evap = cell[veg][band].layer[lidx].evap;
        }
        all_vars_crop->cell[idx][band].inflow = cell[veg][band].inflow;
	for ( lidx = 0; lidx < N_PET_TYPES; lidx++ ) {
          all_vars_crop->cell[idx][band].pot_evap[lidx] = cell[veg][band].pot_evap[lidx];
        }
        all_vars_crop->cell[idx][band].runoff = cell[veg][band].runoff;
        all_vars_crop->cell[idx][band].RhLitter = cell[veg][band].RhLitter;
        all_vars_crop->cell[idx][band].RhLitter2Atm = cell[veg][band].RhLitter2Atm;
        all_vars_crop->cell[idx][band].RhInter = cell[veg][band].RhInter;
        all_vars_crop->cell[idx][band].RhSlow = cell[veg][band].RhSlow;
        all_vars_crop->cell[idx][band].RhTot = cell[veg][band].RhTot;
        all_vars_crop->cell[idx][band].rootmoist = cell[veg][band].rootmoist;
        all_vars_crop->cell[idx][band].wetness = cell[veg][band].wetness;
        all_vars_crop->cell[idx][band].zwt = cell[veg][band].zwt;
        all_vars_crop->cell[idx][band].zwt_lumped = cell[veg][band].zwt_lumped;

        // Copy snow state data
        all_vars_crop->snow[idx][band].albedo = snow[veg][band].albedo;
        all_vars_crop->snow[idx][band].canopy_albedo = snow[veg][band].canopy_albedo;
        all_vars_crop->snow[idx][band].coldcontent = snow[veg][band].coldcontent;
        all_vars_crop->snow[idx][band].coverage = snow[veg][band].coverage;
        all_vars_crop->snow[idx][band].density = snow[veg][band].density;
        all_vars_crop->snow[idx][band].depth = snow[veg][band].depth;
        all_vars_crop->snow[idx][band].last_snow = snow[veg][band].last_snow;
        all_vars_crop->snow[idx][band].max_snow_depth = snow[veg][band].max_snow_depth;
        all_vars_crop->snow[idx][band].MELTING = snow[veg][band].MELTING;
        all_vars_crop->snow[idx][band].pack_temp = snow[veg][band].pack_temp;
        all_vars_crop->snow[idx][band].pack_water = snow[veg][band].pack_water;
        all_vars_crop->snow[idx][band].snow = snow[veg][band].snow;
        if (idx % 2 == 0)
          all_vars_crop->snow[idx][band].snow_canopy = 0;
        else
          all_vars_crop->snow[idx][band].snow_canopy = snow[veg][band].snow_canopy;
        all_vars_crop->snow[idx][band].store_coverage = snow[veg][band].store_coverage;
        all_vars_crop->snow[idx][band].store_snow = snow[veg][band].store_snow;
        all_vars_crop->snow[idx][band].store_swq = snow[veg][band].store_swq;
        all_vars_crop->snow[idx][band].surf_temp = snow[veg][band].surf_temp;
        all_vars_crop->snow[idx][band].surf_temp_fbcount = snow[veg][band].surf_temp_fbcount;
        all_vars_crop->snow[idx][band].surf_temp_fbflag = snow[veg][band].surf_temp_fbflag;
        all_vars_crop->snow[idx][band].surf_water = snow[veg][band].surf_water;
        all_vars_crop->snow[idx][band].swq = snow[veg][band].swq;
        all_vars_crop->snow[idx][band].snow_distrib_slope = snow[veg][band].snow_distrib_slope;
        all_vars_crop->snow[idx][band].tmp_int_storage = snow[veg][band].tmp_int_storage;

        // Copy snow flux data
        all_vars_crop->snow[idx][band].blowing_flux = snow[veg][band].blowing_flux;
        all_vars_crop->snow[idx][band].canopy_vapor_flux = snow[veg][band].canopy_vapor_flux;
        all_vars_crop->snow[idx][band].mass_error = snow[veg][band].mass_error;
        all_vars_crop->snow[idx][band].melt = snow[veg][band].melt;
        all_vars_crop->snow[idx][band].Qnet = snow[veg][band].Qnet;
        all_vars_crop->snow[idx][band].surface_flux = snow[veg][band].surface_flux;
        all_vars_crop->snow[idx][band].transport = snow[veg][band].transport;
        all_vars_crop->snow[idx][band].vapor_flux = snow[veg][band].albedo;

        // Copy energy state data
        all_vars_crop->energy[idx][band].AlbedoLake = energy[veg][band].AlbedoLake;
        all_vars_crop->energy[idx][band].AlbedoOver = energy[veg][band].AlbedoOver;
        all_vars_crop->energy[idx][band].AlbedoUnder = energy[veg][band].AlbedoUnder;
        all_vars_crop->energy[idx][band].frozen = energy[veg][band].frozen;
        all_vars_crop->energy[idx][band].Nfrost = energy[veg][band].Nfrost;
        all_vars_crop->energy[idx][band].Nthaw = energy[veg][band].Nthaw;
        all_vars_crop->energy[idx][band].T1_index = energy[veg][band].T1_index;
        all_vars_crop->energy[idx][band].Tcanopy = energy[veg][band].Tcanopy;
        all_vars_crop->energy[idx][band].Tcanopy_fbcount = energy[veg][band].Tcanopy_fbcount;
        all_vars_crop->energy[idx][band].Tcanopy_fbflag = energy[veg][band].Tcanopy_fbflag;
        all_vars_crop->energy[idx][band].Tfoliage = energy[veg][band].Tfoliage;
        all_vars_crop->energy[idx][band].Tfoliage_fbcount = energy[veg][band].Tfoliage_fbcount;
        all_vars_crop->energy[idx][band].Tfoliage_fbflag = energy[veg][band].Tfoliage_fbflag;
        all_vars_crop->energy[idx][band].Tsurf = energy[veg][band].Tsurf;
        all_vars_crop->energy[idx][band].Tsurf_fbcount = energy[veg][band].Tsurf_fbcount;
        all_vars_crop->energy[idx][band].Tsurf_fbflag = energy[veg][band].Tsurf_fbflag;
        all_vars_crop->energy[idx][band].unfrozen = energy[veg][band].unfrozen;
        for (lidx=0; lidx<2; lidx++) {
          all_vars_crop->energy[idx][band].Cs[lidx] = energy[veg][band].Cs[lidx];
          all_vars_crop->energy[idx][band].kappa[lidx] = energy[veg][band].kappa[lidx];
        }
        for (index=0; index<Nnodes; index++) {
          all_vars_crop->energy[idx][band].Cs_node[index] = energy[veg][band].Cs_node[index];
          all_vars_crop->energy[idx][band].fdepth[index] = energy[veg][band].fdepth[index];
          all_vars_crop->energy[idx][band].ice[index] = energy[veg][band].ice[index];
          all_vars_crop->energy[idx][band].kappa_node[index] = energy[veg][band].kappa_node[index];
          all_vars_crop->energy[idx][band].moist[index] = energy[veg][band].moist[index];
          all_vars_crop->energy[idx][band].T[index] = energy[veg][band].T[index];
          all_vars_crop->energy[idx][band].T_fbcount[index] = energy[veg][band].T_fbcount[index];
          all_vars_crop->energy[idx][band].T_fbflag[index] = energy[veg][band].T_fbflag[index];
          all_vars_crop->energy[idx][band].tdepth[index] = energy[veg][band].tdepth[index];
        }

        // Copy energy flux data
        all_vars_crop->energy[idx][band].advected_sensible = energy[veg][band].advected_sensible;
        all_vars_crop->energy[idx][band].advection = energy[veg][band].advection;
        all_vars_crop->energy[idx][band].AtmosError = energy[veg][band].AtmosError;
        all_vars_crop->energy[idx][band].AtmosLatent = energy[veg][band].AtmosLatent;
        all_vars_crop->energy[idx][band].AtmosLatentSub = energy[veg][band].AtmosLatentSub;
        all_vars_crop->energy[idx][band].AtmosSensible = energy[veg][band].AtmosSensible;
        all_vars_crop->energy[idx][band].canopy_advection = energy[veg][band].canopy_advection;
        all_vars_crop->energy[idx][band].canopy_latent = energy[veg][band].canopy_latent;
        all_vars_crop->energy[idx][band].canopy_latent_sub = energy[veg][band].canopy_latent_sub;
        all_vars_crop->energy[idx][band].canopy_refreeze = energy[veg][band].canopy_refreeze;
        all_vars_crop->energy[idx][band].canopy_sensible = energy[veg][band].canopy_sensible;
        all_vars_crop->energy[idx][band].deltaCC = energy[veg][band].deltaCC;
        all_vars_crop->energy[idx][band].deltaH = energy[veg][band].deltaH;
        all_vars_crop->energy[idx][band].error = energy[veg][band].error;
        all_vars_crop->energy[idx][band].fusion = energy[veg][band].fusion;
        all_vars_crop->energy[idx][band].grnd_flux = energy[veg][band].grnd_flux;
        all_vars_crop->energy[idx][band].latent = energy[veg][band].latent;
        all_vars_crop->energy[idx][band].latent_sub = energy[veg][band].latent_sub;
        all_vars_crop->energy[idx][band].longwave = energy[veg][band].longwave;
        all_vars_crop->energy[idx][band].LongOverIn = energy[veg][band].LongOverIn;
        all_vars_crop->energy[idx][band].LongUnderIn = energy[veg][band].LongUnderIn;
        all_vars_crop->energy[idx][band].LongUnderOut = energy[veg][band].LongUnderOut;
        all_vars_crop->energy[idx][band].melt_energy = energy[veg][band].melt_energy;
        all_vars_crop->energy[idx][band].NetLongAtmos = energy[veg][band].NetLongAtmos;
        all_vars_crop->energy[idx][band].NetLongOver = energy[veg][band].NetLongOver;
        all_vars_crop->energy[idx][band].NetLongUnder = energy[veg][band].NetLongUnder;
        all_vars_crop->energy[idx][band].NetShortAtmos = energy[veg][band].NetShortAtmos;
        all_vars_crop->energy[idx][band].NetShortGrnd = energy[veg][band].NetShortGrnd;
        all_vars_crop->energy[idx][band].NetShortOver = energy[veg][band].NetShortOver;
        all_vars_crop->energy[idx][band].NetShortUnder = energy[veg][band].NetShortUnder;
        all_vars_crop->energy[idx][band].out_long_canopy = energy[veg][band].out_long_canopy;
        all_vars_crop->energy[idx][band].out_long_surface = energy[veg][band].out_long_surface;
        all_vars_crop->energy[idx][band].refreeze_energy = energy[veg][band].refreeze_energy;
        all_vars_crop->energy[idx][band].sensible = energy[veg][band].sensible;
        all_vars_crop->energy[idx][band].shortwave = energy[veg][band].shortwave;
        all_vars_crop->energy[idx][band].ShortOverIn = energy[veg][band].ShortOverIn;
        all_vars_crop->energy[idx][band].ShortUnderIn = energy[veg][band].ShortUnderIn;
        all_vars_crop->energy[idx][band].snow_flux = energy[veg][band].snow_flux;

      }
      }
    }
  }

  return(0);
}


int update_thermal_nodes(all_vars_struct     *all_vars,
			 int                  Nveg,
			 int                  Nnodes,
			 soil_con_struct     *soil_con,
			 veg_con_struct      *veg_con)
/**********************************************************************
  update_thermal_nodes           Jennifer Adam        August 16, 2007

  This routine is run after subsidence occurs (used only for EXCESS_ICE option).
  This routine updates the node depths and interpolates the current
  node temperatures to the new depths, then recalculates the nodal
  thermal properties.  Much of this routine is taken directly from
  initialize_model_state.

  Modifications:
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
**********************************************************************/
{
  extern option_struct options;
  extern veg_lib_struct *veg_lib;
  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      veg, index;
  int      lidx;
  int      band;
  int      ErrorFlag;
  double   Cv;
  double   Zsum, dp;
  double   tmpdp, tmpadj, Bexp;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];

  cell_data_struct      **cell;
  energy_bal_struct     **energy;

  double Tnode_prior[MAX_NODES];
  double Zsum_prior[MAX_NODES];

  cell    = all_vars->cell;
  energy  = all_vars->energy;
  
  dp = soil_con->dp;

  FIRST_VEG = TRUE;

  /*****************************************************************
    Update soil thermal node depths, thicknesses, and temperatures.
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
  *****************************************************************/

  /*****************************************************************
    Update soil thermal node depths and thicknesses.
  *****************************************************************/
  //set previous Zsum
  for ( index = 0; index < Nnodes; index++ ) 
    Zsum_prior[index] = soil_con->Zsum_node[index];

  if(!options.EXP_TRANS){  
    /* Nodes set at surface, the depth of the first layer,
       twice the depth of the first layer, and at the
       damping depth.  Extra nodes are placed equal distance
       between the damping depth and twice the depth of the
       first layer. */
    
    soil_con->dz_node[0] = soil_con->depth[0];
    soil_con->dz_node[1] = soil_con->depth[0];
    soil_con->dz_node[2] = soil_con->depth[0];	  
    soil_con->Zsum_node[0] = 0;
    soil_con->Zsum_node[1] = soil_con[0].depth[0];
    Zsum   = 2. * soil_con[0].depth[0];
    soil_con->Zsum_node[2] = Zsum;
    tmpdp  = dp - soil_con[0].depth[0] * 2.5;
    tmpadj = 3.5;
    for ( index = 3; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
      Zsum += (soil_con->dz_node[index]
	       +soil_con->dz_node[index-1])/2.;
      soil_con->Zsum_node[index] = Zsum;
    }
    soil_con->dz_node[Nnodes-1] = (dp - Zsum 
				   - soil_con->dz_node[Nnodes-2] 
				   / 2. ) * 2.;
    Zsum += (soil_con->dz_node[Nnodes-2]
	     +soil_con->dz_node[Nnodes-1])/2.;
    soil_con->Zsum_node[Nnodes-1] = Zsum;
    if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
      sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
      nrerror(ErrStr);
    }
  }
  else{ /* exponential grid transformation, EXP_TRANS = TRUE*/
    
    /*calculate exponential function parameter */
    Bexp = logf(dp+1.)/(double)(Nnodes-1); //to force Zsum=dp at bottom node
    for ( index = 0; index <= Nnodes-1; index++ )
      soil_con->Zsum_node[index] = expf(Bexp*index)-1.;
    
    //top node	  
    index=0;
    soil_con->dz_node[index] = soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index];
    //middle nodes
    for ( index = 1; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = (soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index])/2.+(soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1])/2.;
    }
    //bottom node
    index=Nnodes-1;
    soil_con->dz_node[index] = soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1];
  }
#if VERBOSE
  fprintf(stderr,"More updated parameters in soil_con: dz_node and Zsum_node.\n");
#endif

  /******************************************
    Update soil thermal node temperatures via linear interpolation.
  ******************************************/
  for ( veg = 0 ; veg <= Nveg ; veg++ ) {
    // Initialize soil for existing vegetation types
    Cv = veg_con[veg].Cv;
    
    if ( Cv > 0 ) {
      for( band = 0; band < options.SNOW_BAND; band++ ) {
	if ( soil_con->AreaFract[band] > 0. ) {
	  //set previous temperatures
	  for ( index = 0; index < Nnodes; index++ ) 
	    Tnode_prior[index] = energy[veg][band].T[index];
	  //top node: no need to update surface temperature
	  //remaining nodes
	  for ( index = 1; index < Nnodes; index++ ) {
	    energy[veg][band].T[index] = linear_interp(soil_con->Zsum_node[index],Zsum_prior[index-1],Zsum_prior[index],Tnode_prior[index-1],Tnode_prior[index]);	
	  }//node
	}	
      }//band
    }
  }//veg

  /******************************************
    Update soil thermal node properties 
  ******************************************/  
  FIRST_VEG = TRUE;
  for ( veg = 0 ; veg <= Nveg ; veg++) {
    // Initialize soil for existing vegetation types
    Cv = veg_con[veg].Cv;

    if ( Cv > 0 ) {
      for( band = 0; band < options.SNOW_BAND; band++ ) {
	// Initialize soil for existing snow elevation bands
	if ( soil_con->AreaFract[band] > 0. ) {
	  /** Set soil properties for all soil nodes **/
	  if(FIRST_VEG) {
	    FIRST_VEG = FALSE;
	    set_node_parameters(soil_con->dz_node, soil_con->Zsum_node, soil_con->max_moist_node,
				  soil_con->expt_node, soil_con->bubble_node,
				  soil_con->alpha, soil_con->beta,
				  soil_con->gamma, soil_con->depth,
				  soil_con->max_moist, soil_con->expt, 
				  soil_con->bubble, soil_con->quartz, 
				  Nnodes, options.Nlayer, soil_con->FS_ACTIVE);	  
	  }

	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) 
	    moist[veg][band][lidx] = cell[veg][band].layer[lidx].moist;

	  /* set soil moisture properties for all soil thermal nodes */
	  if ( !( options.LAKES && veg_con->LAKE != 0 ) ) {
	    ErrorFlag = distribute_node_moisture_properties(energy[veg][band].moist,
						  energy[veg][band].ice,
						  energy[veg][band].kappa_node,
						  energy[veg][band].Cs_node,
						  soil_con->Zsum_node,
						  energy[veg][band].T,
						  soil_con->max_moist_node,
						  soil_con->expt_node,
						  soil_con->bubble_node,
						  moist[veg][band],
						  soil_con->depth,
						  soil_con->soil_dens_min,
						  soil_con->bulk_dens_min,
						  soil_con->quartz,
						  soil_con->soil_density,
						  soil_con->bulk_density,
						  soil_con->organic,
						  Nnodes, options.Nlayer,
						  soil_con->FS_ACTIVE);
	    if ( ErrorFlag == ERROR ) return ( ErrorFlag );
	  }

	  /* initialize layer moistures and ice contents */
	  if ( !( options.LAKES && veg_con->LAKE != 0 ) ) {
            if (options.QUICK_FLUX) {
              ErrorFlag = estimate_layer_ice_content_quick_flux(cell[veg][band].layer,
					 soil_con->depth, soil_con->dp,
					 energy[veg][band].T[0], energy[veg][band].T[1],
					 soil_con->avg_temp, soil_con->max_moist, 
					 soil_con->expt, soil_con->bubble, 
					 soil_con->frost_fract, soil_con->frost_slope, 
					 soil_con->FS_ACTIVE);
            }
            else {
	      ErrorFlag = estimate_layer_ice_content(cell[veg][band].layer,
						     soil_con->Zsum_node,
						     energy[veg][band].T,
						     soil_con->max_moist_node,
						     soil_con->expt_node,
						     soil_con->bubble_node,
						     soil_con->depth,
						     soil_con->max_moist,
						     soil_con->expt,
						     soil_con->bubble,
						     soil_con->frost_fract, 
						     soil_con->frost_slope, 
						     Nnodes, options.Nlayer, 
						     soil_con->FS_ACTIVE);	      
	    }
	  }
	    
	  /* Find freezing and thawing front depths */
	  if(!options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	    if ( !( options.LAKES && veg_con->LAKE != 0 ) ) 
	      find_0_degree_fronts(&energy[veg][band], soil_con->Zsum_node, energy[veg][band].T, Nnodes);
	}
      }//band
    }
  }//veg	

  return(0);  
}
