#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

static char vcid[] = "$Id$";

/** Main Program **/

int main(int argc, char *argv[])
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
	  model.                                                  KAC
  02-27-01 added controls for lake model                          KAC
  11-18-02 Updated storage of lake water for water balance 
           calculations.                                          LCB
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  04-10-03 Modified to initialize storm parameters using the state
           file.                                                KAC
  04-10-03 Modified to start the model by skipping records until the
           state file date is found.  This replaces the previous method
           of modifying the global file start date, which can change 
           the interpolation of atmospheric forcing data.        KAC
  04-15-03 Modified to store wet and dry fractions when intializing 
           water balance storage.  This accounts for changes in model
           state initialization, which now stores wet and dry fractions
           rather than just averagedvalues.                      KAC
  29-Oct-03 Modified the version display banner to print the version
	    string defined in global.h.				TJB
  01-Nov-04 Updated arglist for make_dist_prcp(), as part of fix for
	    QUICK_FLUX state file compatibility.		TJB
  02-Nov-04 Updated arglist for read_lakeparam(), as part of fix for
	    lake fraction readjustment.				TJB

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif // LINK_DEBUG
  extern Error_struct Error;
  extern global_param_struct global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                    *init_STILL_STORM;
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                     *init_DRY_TIME;
  int                      RUN_MODEL;
  int                      Ncells;
  int                      cell_cnt;
  int                      startrec;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_con_struct          *veg_con;
  soil_con_struct          soil_con;
  dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */
  filenames_struct         filenames;
  filenames_struct         builtnames;
  infiles_struct           infiles;
#if LAKE_MODEL
  lake_con_struct          lake_con;
#endif // LAKE_MODEL
  outfiles_struct          outfiles;
  
  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

#if VERBOSE
  display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
#endif

  /** Read Global Control File **/
  infiles.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames, infiles.globalparam);

  /** Check and Open Files **/
  check_files(&infiles, &filenames);

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(infiles.veglib,&Nveg_type);

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /** Make Date Data Structure **/
  dmy      = make_dmy(&global_param);

  /** allocate memory for the atmos_data_struct **/
  alloc_atmos(global_param.nrecs, &atmos);

#if SAVE_STATE
  /** open state file if model state is to be saved **/
  if ( strcmp( global_param.statename, "NONE" ) != 0 ) 
    outfiles.statefile = open_state_file(&global_param, options.Nlayer, 
					 options.Nnode);
  else outfiles.statefile = NULL;
#endif

  if ( options.INIT_STATE ) 
    infiles.statefile = check_state_file(filenames.init_state, dmy, 
					 &global_param, options.Nlayer, 
					 options.Nnode, &startrec);
  else
    startrec = 0;

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  cell_cnt=0;
  while(!MODEL_DONE) {
    if(!options.ARC_SOIL) {
      if((fscanf(infiles.soilparam, "%d", &flag))!=EOF) {
	if(flag) RUN_MODEL=TRUE;
	else     RUN_MODEL=FALSE;
      }
      else {
	MODEL_DONE = TRUE;
	RUN_MODEL = FALSE;
      }
      if(!MODEL_DONE) soil_con = read_soilparam(infiles.soilparam, RUN_MODEL);
    }
    else {
      soil_con = read_soilparam_arc(infiles.soilparam, 
				    filenames.soil_dir, &Ncells, 
				    &RUN_MODEL, cell_cnt);
      cell_cnt++;
      if(cell_cnt==Ncells) MODEL_DONE = TRUE;
    }
    if(RUN_MODEL) {
#if LINK_DEBUG
      if(debug.PRT_SOIL) write_soilparam(&soil_con); 
#endif

#if !OUTPUT_FORCE
#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
	}
	for(i=0;i<MAX_NODES;i++) {
	  soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
	}
      }
#endif /* QUICK_FS */

      NEWCELL=TRUE;
      cellnum++;

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = read_vegparam(infiles.vegparam, soil_con.gridcel,
                              Nveg_type);
      calc_root_fractions(veg_con, &soil_con);
#if LINK_DEBUG
      if(debug.PRT_VEGE) write_vegparam(veg_con); 
#endif /* LINK_DEBUG*/
#endif /* !OUTPUT_FORCE */

#if LAKE_MODEL
      if ( options.LAKES ) 
	lake_con = read_lakeparam(infiles.lakeparam, soil_con, 
				  veg_con, global_param.resolution);
#endif // LAKE_MODEL

      /** Build Gridded Filenames, and Open **/
      builtnames = make_in_and_outfiles(&infiles, &filenames, &soil_con,
                   &outfiles);

#if !OUTPUT_FORCE
      /** Read Elevation Band Data if Used **/
      read_snowband(infiles.snowband,soil_con.gridcel,
		    (double)soil_con.elevation, &soil_con.Tfactor, 
		    &soil_con.Pfactor, &soil_con.AreaFract, 
		    &soil_con.AboveTreeLine);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */
#endif /* !OUTPUT_FORCE */

      initialize_atmos(atmos, dmy, infiles.forcing,
		       (double)soil_con.time_zone_lng, (double)soil_con.lng,
		       (double)soil_con.lat, soil_con.elevation,
		       soil_con.annual_prec, global_param.wind_h, 
		       soil_con.rough, soil_con.Tfactor, 
#if OUTPUT_FORCE
		       soil_con.AboveTreeLine, &outfiles); 
#else /* OUTPUT_FORCE */
                       soil_con.AboveTreeLine); 
#endif /* OUTPUT_FORCE */

#if LINK_DEBUG
      if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);
#endif

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

#if !OUTPUT_FORCE
#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
      initialize_model_state(atmos[0].air_temp[NR], Ndist, options.Nnode, 
			     veg_con[0].vegetat_type_num, soil_con.gridcel, 
			     &prcp, dmy[0], &global_param, infiles, 
#if LAKE_MODEL
			     lake_con,
#endif // LAKE_MODEL
			     &soil_con, veg_con, &init_STILL_STORM,
			     &init_DRY_TIME);


#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

      /** Update Error Handling Structure **/
      Error.outfp = outfiles;
      Error.infp = infiles;

      /***************************************************
	Intialize Moisture and Energy Balance Error Checks
        --- As of 4/15/03 this does not properly initialize
            storage from bands above treeline, when the model 
            state is restored from a file.  This can lead to 
            water balance errors in the initial time step but 
            does not impact the actual simulation.  It will
            be addressed in the next release version.  KAC
	***************************************************/
      storage = 0.;
      for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
	if ( veg < veg_con[0].vegetat_type_num ) veg_fract = veg_con[veg].Cv;
	else veg_fract = ( 1.0 - veg_con[0].Cv_sum );
	for ( band = 0; band < options.SNOW_BAND; band++ ) {
	  band_fract = soil_con.AreaFract[band];
	  if ( veg_fract > SMALL && band_fract > SMALL ) {
	    for(index=0;index<options.Nlayer;index++)
	      for ( dist = 0; dist < Ndist; dist ++ )
		storage += prcp.cell[dist][veg][band].layer[index].moist 
		  * veg_fract * band_fract;
	    storage += prcp.snow[veg][band].swq * 1000. * veg_fract 
	      * band_fract;
	    if ( veg != veg_con[0].vegetat_type_num ) {
	      for ( dist = 0; dist < Ndist; dist ++ ) 
		storage += prcp.veg_var[dist][veg][band].Wdew 
		  * veg_fract * band_fract;
	      storage += prcp.snow[veg][band].snow_canopy * 1000. 
		* veg_fract * band_fract;
	    }
	  }
	}
      }

#if LAKE_MODEL
      if ( options.LAKES && lake_con.Cl[0] > 0) {
	/** COMPUTE MOISTURE STORAGE IN LAKE FRACTION **/
	storage += lake_con.Cl[0] * (prcp.lake_var.volume 
				     / lake_con.basin[0]) * 1000.;
	band = 0;
	veg = veg_con[0].vegetat_type_num + 1;
	for(index=0;index<options.Nlayer;index++)
	  storage += lake_con.Cl[0] * prcp.cell[WET][veg][band].layer[index].moist; 
      }
#endif // LAKE_MODEL

      calc_water_balance_error(-global_param.nrecs,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

        if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;

        dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
                  dmy, &global_param, 
#if LAKE_MODEL
		  &lake_con, 
#endif /* LAKE_MODEL */
		  &outfiles, rec, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);
        NEWCELL=FALSE;
	for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
	  init_DRY_TIME[veg] = -999;
	if ( options.LAKES )
	  init_DRY_TIME[veg_con[0].vegetat_type_num+1] = -999;

      }	/* End Rec Loop */

      close_files(&infiles,&outfiles,&builtnames); 

#if QUICK_FS
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_layer[i][j]);
	  free((char *)soil_con.ufwc_table_layer[i]);
	}
	for(i=0;i<MAX_NODES;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_node[i][j]);
	  free((char *)soil_con.ufwc_table_node[i]);
	}
      }
#endif /* QUICK_FS */
      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      free_vegcon(&veg_con);
      free((char *)soil_con.AreaFract);
      free((char *)soil_con.Tfactor);
      free((char *)soil_con.Pfactor);
      free((char *)soil_con.AboveTreeLine);
#endif /* !OUTPUT_FORCE */
      for(index=0;index<=options.Nlayer;index++) 
	free((char*)soil_con.layer_node_fract[index]);
      free((char*)soil_con.layer_node_fract);
    }	/* End Run Model Condition */
  } 	/* End Grid Loop */

  /** cleanup **/
  free_atmos(global_param.nrecs, &atmos);

  return EXIT_SUCCESS;
}	/* End Main Program */
