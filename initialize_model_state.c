#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

void initialize_model_state(dist_prcp_struct    *prcp,
			    dmy_struct           dmy,
                            global_param_struct *global_param,
                            filep_struct         filep,
			    int                  cellnum,
		            int                  Nveg,
                            int                  Nnodes,
			    int                  Ndist,
			    double               surf_temp, 
                            soil_con_struct     *soil_con,
			    veg_con_struct      *veg_con,
			    lake_con_struct      lake_con,
			    char               **init_STILL_STORM,
			    int                **init_DRY_TIME,
			    save_data_struct    *save_data)
/**********************************************************************
  initialize_model_state       Keith Cherkauer	    April 17, 2000

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
	    accessed in set_node_parameters().			TJB
  01-Nov-04 Added support for state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.			TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              save_data structure to track changes in moisture storage
              over each time step; this needs initialization here.  TJB
  2006-Oct-10 Added snow[veg][band].snow_canopy to save_data.swe. TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included removing the unused init_snow file. TJB
  2006-Nov-07 Removed LAKE_MODEL option. TJB

**********************************************************************/
{
  extern option_struct options;
  extern veg_lib_struct *veg_lib;
#if LINK_DEBUG
  extern debug_struct debug;
#endif
#if QUICK_FS
  extern double temps[];
#endif

  char     tmpstr[MAXSTRING];
  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      i, j, ii, veg, index, dist;
  int      nidx, lidx;
  int      tmpint;
  int      dry;
  int      band;
  int      zindex;
  int      MaxVeg;
#if SPATIAL_FROST
  int      frost_area;
#endif
  double   Cv;
  double   sum, Lsum, Zsum, dp, Ltotal;
  double   tmpdp, tmpadj;
  double   Tair;
  double  *kappa, *Cs, *M;
  double   moist[MAX_VEG][MAX_BANDS][MAX_LAYERS];
#if SPATIAL_FROST
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS][FROST_SUBAREAS];
#else
  double   ice[MAX_VEG][MAX_BANDS][MAX_LAYERS];
#endif // SPATIAL_FROST
  double   unfrozen, frozen;
  double **layer_ice;
  double **layer_tmp;
  double  *EMPTY;
#if QUICK_FS
  double   Aufwc, Bufwc;
#endif
  char    *EMPTY_C;
  double   Clake;
  double   mu;
  double   TreeAdjustFactor[MAX_BANDS];

  cell_data_struct     ***cell;
  energy_bal_struct     **energy;
  lake_var_struct   *lake_var;
  snow_data_struct      **snow;
  veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  energy  = prcp->energy;
  lake_var = &prcp->lake_var;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;
  
  // Initialize soil depths
  dp = soil_con->dp;
  Ltotal = 0;
  for ( index = 0; index < options.Nlayer; index++ ) 
    Ltotal += soil_con->depth[index];

  FIRST_VEG = TRUE;

  if(options.LAKES) 
    MaxVeg = Nveg+1;
  else
    MaxVeg = Nveg;

  // increase initial soil surface temperature if air is very cold
  Tair = surf_temp;
  if ( surf_temp < -1. ) surf_temp = -1.;
  
  // initialize storm parameters to start a new simulation
  (*init_STILL_STORM) = (char *)malloc((MaxVeg+1)*sizeof(char));
  (*init_DRY_TIME)    = (int *)malloc((MaxVeg+1)*sizeof(int));
  for ( veg = 0 ; veg <= MaxVeg ; veg++ )
    (*init_DRY_TIME)[veg] = -999;
  
  /********************************************
    Initialize all snow pack variables 
    - some may be reset if state file present
  ********************************************/

  initialize_snow(snow, MaxVeg, cellnum);

  /********************************************
    Initialize all soil layer variables 
    - some may be reset if state file present
  ********************************************/

  initialize_soil(cell[WET], soil_con, MaxVeg);
  if ( options.DIST_PRCP )
    initialize_soil(cell[DRY], soil_con, MaxVeg);

  /********************************************
    Initialize all vegetation variables 
    - some may be reset if state file present
  ********************************************/

  initialize_veg(veg_var[WET], veg_con, global_param, MaxVeg);
  if ( options.DIST_PRCP )
    initialize_veg(veg_var[DRY], veg_con, global_param, MaxVeg);

  /********************************************
    Initialize all lake variables 
  ********************************************/

  if ( options.LAKES && lake_con.Cl[0] > 0) {
    initialize_lake(lake_var, lake_con, &snow[Nveg+1][0], surf_temp);
  }

  /********************************************
    Initialize all spatial frost variables 
  ********************************************/

#if SPATIAL_FROST
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
    if ( FROST_SUBAREAS == 1 ) soil_con->frost_fract[frost_area] = 1.;
    else if (FROST_SUBAREAS == 2 ) soil_con->frost_fract[frost_area] = 0.5;
    else {
      soil_con->frost_fract[frost_area] = 1. / (FROST_SUBAREAS - 1);
      if ( frost_area == 0 || frost_area == FROST_SUBAREAS-1 ) 
	soil_con->frost_fract[frost_area] /= 2.;
    }
  }
#endif // SPATIAL_FROST

  /********************************************************
    Compute grid cell fractions for all subareas used in 
    spatial distribution of soil frost routines.
  ********************************************************/

#if QUICK_FS
  if(options.FROZEN_SOIL) {

    /***********************************************************
      Prepare table of maximum unfrozen water content values
      - This linearizes the equation for maximum unfrozen water
        content, reducing computation time for the frozen soil
        model.
    ***********************************************************/

    for(lidx=0;lidx<options.Nlayer;lidx++) { 
      for(ii=0;ii<QUICK_FS_TEMPS;ii++) {
	Aufwc = maximum_unfrozen_water(temps[ii], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	Bufwc = maximum_unfrozen_water(temps[ii+1], 1.0, 
				       soil_con->bubble[lidx], 
				       soil_con->expt[lidx]);
	soil_con->ufwc_table_layer[lidx][ii][0] 
	  = linear_interp(0., temps[ii], temps[ii+1], Aufwc, Bufwc);
	soil_con->ufwc_table_layer[lidx][ii][1] 
	  = (Bufwc - Aufwc) / (temps[ii+1] - temps[ii]);
      }
    }
  }  
#endif // QUICK_FS

  /************************************************************************
    CASE 1: Not using quick ground heat flux, and initial conditions files 
    provided
  ************************************************************************/

  if(options.INIT_STATE) {

    read_initial_model_state(filep.init_state, prcp, global_param,  
			     Nveg, options.SNOW_BAND, cellnum, soil_con,
			     Ndist, *init_STILL_STORM, *init_DRY_TIME, lake_con);

    for ( veg = 0 ; veg <= MaxVeg ; veg++ ) {
      // Initialize soil for existing vegetation types
      if ( veg < Nveg ) Cv = veg_con[veg].Cv;
      else Cv = (1.0 - veg_con[0].Cv_sum);

      if ( Cv > 0 || ( veg == MaxVeg && MaxVeg > Nveg ) ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  for( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] 
		= cell[0][veg][band].layer[lidx].ice[frost_area];
#else
	    ice[veg][band][lidx] = cell[0][veg][band].layer[lidx].ice;
#endif
	  }
	}
      }
    }
  }
  
  /************************************************************************
    CASE 2: Initialize soil if using quick heat flux, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(options.QUICK_FLUX) {

    for ( veg = 0 ; veg <= MaxVeg ; veg++ ) {
      // Initialize soil for existing vegetation types
      if ( veg < Nveg ) Cv = veg_con[veg].Cv;
      else Cv = (1.0 - veg_con[0].Cv_sum);
      
      if ( Cv > 0 || ( veg == MaxVeg && MaxVeg > Nveg ) ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {

	  /* Initialize soil node temperatures and thicknesses */
	  
	  soil_con->dz_node[0]   = soil_con->depth[0];
	  soil_con->dz_node[1]   = soil_con->depth[0];
	  soil_con->dz_node[2]   = 2. * (dp - 1.5 * soil_con->depth[0]);
	  energy[veg][band].T[0] = surf_temp;
	  energy[veg][band].T[1] = surf_temp;
	  energy[veg][band].T[2] = soil_con->avg_temp;
	  
	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
#else
	    ice[veg][band][lidx] = 0.;
#endif
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
    for ( veg = 0 ; veg <= MaxVeg ; veg++ ) {
      // Initialize soil for existing vegetation types
      if ( veg < Nveg ) Cv = veg_con[veg].Cv;
      else Cv = (1.0 - veg_con[0].Cv_sum);
      
      if ( Cv > 0 || ( veg == MaxVeg && MaxVeg > Nveg ) ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  
	  /* Initialize soil node temperatures and thicknesses 
	     Nodes set at surface, the depth of the first layer,
	     twice the depth of the first layer, and at the
	     damping depth.  Extra nodes are placed equal distance
	     between the damping depth and twice the depth of the
	     first layer. */
	  
	  energy[veg][band].T[0] = surf_temp;
	  soil_con->dz_node[0] = soil_con->depth[0];
	  soil_con->dz_node[1] = soil_con->depth[0];
	  soil_con->dz_node[2] = soil_con->depth[0];
	  energy[veg][band].T[Nnodes-1] = soil_con->avg_temp;
	  energy[veg][band].T[1] = exp_interp(soil_con->depth[0], 0., dp, 
					      surf_temp, soil_con->avg_temp);
	  energy[veg][band].T[2] = exp_interp(2. * soil_con->depth[0], 0., dp, 
					      surf_temp, soil_con->avg_temp);
	  
	  Zsum   = 2. * soil_con[0].depth[0];
	  tmpdp  = dp - soil_con[0].depth[0] * 2.5;
	  tmpadj = 3.5;
	  for ( index = 3; index < Nnodes-1; index++ ) {
	    if ( FIRST_VEG ) {
	      soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
	    }
	    Zsum += (soil_con->dz_node[index]
		     +soil_con->dz_node[index-1])/2.;
	    energy[veg][band].T[index] = exp_interp(Zsum,0.,soil_con[0].dp,
						    surf_temp,
						    soil_con[0].avg_temp);
	  }
	  if ( FIRST_VEG ) {
	    FIRST_VEG = FALSE;
	    soil_con->dz_node[Nnodes-1] = (dp - Zsum 
					   - soil_con->dz_node[Nnodes-2] 
					   / 2. ) * 2.;
	    Zsum += (soil_con->dz_node[Nnodes-2]
		     +soil_con->dz_node[Nnodes-1])/2.;
	    if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
	      sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,dp);
	      nrerror(ErrStr);
	    }
	  }

	  for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	    moist[veg][band][lidx] = cell[0][veg][band].layer[lidx].moist;
#if SPATIAL_FROST
	    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
	      ice[veg][band][lidx][frost_area] = 0.;
#else
	    ice[veg][band][lidx] = 0.;
#endif
	  }
	}
      }
    }
  }

  /*********************************
    CASE 4: Unknown option
  *********************************/
  else {
    for ( veg = 0 ; veg <= MaxVeg ; veg++ ) {
      // Initialize soil for existing vegetation types
      if ( veg < Nveg ) Cv = veg_con[veg].Cv;
      else Cv = (1.0 - veg_con[0].Cv_sum);

      if ( Cv > 0 || ( veg == MaxVeg && MaxVeg > Nveg ) ) {
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

/* dz_node[Nnodes] is accessed later despite not being set. This can
  cause run-time errors on some platforms. Therefore, set it to zero */
  soil_con->dz_node[Nnodes]=0.0;

  if ( options.GRND_FLUX ) {

    FIRST_VEG = TRUE;
    for ( veg = 0 ; veg <= MaxVeg ; veg++) {
      // Initialize soil for existing vegetation types
      if ( veg < Nveg ) Cv = veg_con[veg].Cv;
      else Cv = (1.0 - veg_con[0].Cv_sum);

      if ( Cv > 0 || ( veg == MaxVeg && MaxVeg > Nveg ) ) {
	for( band = 0; band < options.SNOW_BAND; band++ ) {
	  // Initialize soil for existing snow elevation bands
	  if ( soil_con->AreaFract[band] > 0. ) {
	    
	    /** Set soil properties for all soil nodes **/
	    if(FIRST_VEG) {
	      FIRST_VEG = FALSE;
	      set_node_parameters(soil_con->dz_node, soil_con->max_moist_node,
				  soil_con->expt_node, soil_con->bubble_node,
				  soil_con->alpha, soil_con->beta,
				  soil_con->gamma, soil_con->depth,
				  soil_con->max_moist, soil_con->expt, 
				  soil_con->bubble, soil_con->quartz, 
				  soil_con->layer_node_fract,
#if QUICK_FS
				  soil_con->ufwc_table_node,
#endif // QUICK_FS
				  Nnodes, options.Nlayer, soil_con->FS_ACTIVE);
	  
	      sum = soil_con->dz_node[0]/2. + soil_con->dz_node[Nnodes-1]/2.;
	      for ( nidx = 1; nidx < Nnodes-1; nidx++ ) 
		sum += soil_con->dz_node[nidx];
	    }
	
	    /* set soil moisture properties for all soil thermal nodes */
	    distribute_node_moisture_properties(energy[veg][band].moist,
						energy[veg][band].ice,
						energy[veg][band].kappa_node,
						energy[veg][band].Cs_node,
						soil_con->dz_node,
						energy[veg][band].T,
						soil_con->max_moist_node,
#if QUICK_FS
						soil_con->ufwc_table_node,
#else
						soil_con->expt_node,
						soil_con->bubble_node,
#endif // QUICK_FS
						moist[veg][band], 
						soil_con->depth,
						soil_con->soil_density,
						soil_con->bulk_density,
						soil_con->quartz,
						Nnodes, options.Nlayer,
						soil_con->FS_ACTIVE);
	
	    /* initialize layer moistures and ice contents */
	    for ( dry = 0; dry < Ndist; dry++ ) {
	      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
		cell[dry][veg][band].layer[lidx].moist 
		  = moist[veg][band][lidx];
#if SPATIAL_FROST
		for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )

		  cell[dry][veg][band].layer[lidx].ice[frost_area] = ice[veg][band][lidx][frost_area];
#else
		cell[dry][veg][band].layer[lidx].ice = ice[veg][band][lidx];
#endif
	      }
	      if ( !( options.LAKES && veg == MaxVeg ) ) {
		estimate_layer_ice_content(cell[dry][veg][band].layer,
					   soil_con->dz_node,
					   energy[veg][band].T,
					   soil_con->max_moist_node,
#if QUICK_FS
					 soil_con->ufwc_table_node,
#else
					 soil_con->expt_node,
					 soil_con->bubble_node,
#endif // QUICK_FS
					 soil_con->depth,
					 soil_con->max_moist,
#if QUICK_FS
					 soil_con->ufwc_table_layer,
#else
					 soil_con->expt,
					 soil_con->bubble,
#endif // QUICK_FS
#if SPATIAL_FROST
					 soil_con->frost_fract, 
					 soil_con->frost_slope, 
#endif // SPATIAL_FROST
					 soil_con->bulk_density,
					 soil_con->soil_density,
					 soil_con->quartz, 
					 soil_con->layer_node_fract,
					 Nnodes, options.Nlayer, 
					 soil_con->FS_ACTIVE);
	  
	      }
	    }
	
	    /* Find freezing and thawing front depths */
	    if(!options.QUICK_FLUX && soil_con->FS_ACTIVE) 
	      find_0_degree_fronts(&energy[veg][band], soil_con->dz_node,
				   energy[veg][band].T, Nnodes);
	  }
	}
      }
    }	
  }

  // initialize energy balance terms before iterations begin
  for ( veg = 0 ; veg <= MaxVeg ; veg++) {
    for ( band = 0; band < options.SNOW_BAND; band++ ) {
	energy[veg][band].LongUnderOut = 0.;
	energy[veg][band].Tfoliage     = Tair;
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

  // Save initial moisture storages for differencing at end of time step
  save_data->total_soil_moist = 0;
  save_data->swe = 0;
  save_data->wdew = 0;
  save_data->surfstor = 0;
  for( veg = 0; veg <= Nveg; veg++ ) {

    if ( veg < veg_con[0].vegetat_type_num )
      Cv = veg_con[veg].Cv;
    else
      Cv = (1.0 - veg_con[0].Cv_sum);

    if ( Cv > 0 ) {

      for ( dist = 0; dist < Ndist; dist++ ) {

        if(dist==0)
          mu = prcp[0].mu[veg];
        else
          mu = 1. - prcp[0].mu[veg];

        for( band = 0; band < options.SNOW_BAND; band++ ) {

          if (soil_con->AreaFract[band] > 0.
              && ( veg == veg_con[0].vegetat_type_num
                  || ( !soil_con->AboveTreeLine[band]
                      || (soil_con->AboveTreeLine[band] && !veg_lib[veg_con[veg].veg_class].overstory)))) {

            for(lidx=0;lidx<options.Nlayer;lidx++) {
              save_data->total_soil_moist +=
                (cell[dist][veg][band].layer[lidx].moist) // layer[].moist appears to contain both liquid and ice
                * Cv * mu * soil_con->AreaFract[band] * TreeAdjustFactor[band];
            }
            save_data->wdew += veg_var[dist][veg][band].Wdew * Cv * mu * soil_con->AreaFract[band] * TreeAdjustFactor[band];

          }

        }

      }

      for( band = 0; band < options.SNOW_BAND; band++ ) {
        save_data->swe += (snow[veg][band].swq + snow[veg][band].snow_canopy)
                         * Cv * soil_con->AreaFract[band] * TreeAdjustFactor[band];
      }

    }

  }
  Clake = lake_var->sarea/lake_con.basin[0];
  save_data->surfstor = (lake_var->volume/lake_var->sarea)*1000. * Clake * Cv * soil_con->AreaFract[0] * TreeAdjustFactor[0];

}
