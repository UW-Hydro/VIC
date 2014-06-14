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
	    string defined in global.h.					TJB
  01-Nov-04 Updated arglist for make_all_vars(), as part of fix for
	    QUICK_FLUX state file compatibility.			TJB
  02-Nov-04 Updated arglist for read_lakeparam(), as part of fix for
	    lake fraction readjustment.					TJB
  2005-Apr-13 OUTPUT_FORCE option now calls close_files().		TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data, out_data_files, and save_data structures.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included merging builtnames into filenames.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-07 Changed statefile to init_state in call to
	      check_state_file().					TJB
  2007-Jan-15 Added PRT_HEADER option; added call to
	      write_header().						TJB
  2007-Apr-04 Added option to continue run after a cell fails. 		GCT/KAC
  2007-Apr-21 Added calls to free_dmy(), free_out_data_files(),
	      free_out_data(), and free_veglib().  Added closing of
	      all parameter files.					TJB
  2007-Aug-21 Return ErrorFlag from initialize_model_state.		JCA
  2007-Sep-14 Excluded calls to free_veglib() and closing of parameter
	      files other than the soil param file for the case
	      when OUTPUT_FORCE=TRUE.					TJB
  2007-Nov-06 Moved computation of cell_area from read_lakeparam() to
	      read_soilparam() and read_soilparam_arc().		TJB
  2008-May-05 Added all_vars fraction (mu) to initial water storage
	      computation.  This solves water balance errors for the
	      case where DIST_PRCP is TRUE.				TJB
  2009-Jan-16 Added soil_con.avgJulyAirTemp to argument list of
	      initialize_atmos().					TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-07 Added soil_con.BandElev[] to read_snowband() arg list.	TJB
  2009-Jul-31 Replaced references to N+1st veg tile with references
	      to index of lake/wetland tile.				TJB
  2009-Sep-28 Replaced initial water/energy storage computations and
	      calls to calc_water_balance_error/calc_energy_balance_error
	      with an initial call to put_data.  Modified the call to
	      read_snowband().						TJB
  2009-Dec-11 Removed save_data structure from argument list of 
	      initialize_model_state().					TJB
  2010-Mar-31 Added cell_area to initialize_atmos().			TJB
  2010-Apr-28 Removed individual soil_con variables from argument list
	      of initialize_atmos() and replaced with *soil_con.	TJB
  2010-Nov-10 Added closing of state files.				TJB
  2011-Jan-04 Made read_soilparam_arc() a sub-function of
	      read_soilparam().						TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct.			TJB
  2014-Mar-24 Removed ARC_SOIL option                           	BN
  2014-Apr-02 Moved "free" statements for soil_con arrays outside the
	      OUTPUT_FORCE condition to avoid memory leak.		TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-Apr-25 Added non-climatological veg parameters.			TJB
**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern Error_struct Error;
  extern global_param_struct global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                     RUN_MODEL;
  char                     ErrStr[MAXSTRING];
  int                      rec, i, j;
  int                      veg;
  int                      band;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                      Ncells;
  int                      startrec;
  int                      ErrorFlag;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  double                   Clake;
  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_hist_struct        **veg_hist;
  veg_con_struct          *veg_con;
  soil_con_struct          soil_con;
  all_vars_struct          all_vars;
  filenames_struct         filenames;
  filep_struct             filep;
  lake_con_struct          lake_con;
  out_data_file_struct     *out_data_files;
  out_data_struct          *out_data;
  save_data_struct         save_data;
  
  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

#if VERBOSE
  display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
#endif

  /** Read Global Control File **/
  filep.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames, filep.globalparam);

  /** Set up output data structures **/
  out_data = create_output_list();
  out_data_files = set_output_defaults(out_data);
  fclose(filep.globalparam);
  filep.globalparam = open_file(filenames.global,"r");
  parse_output_info(&filenames, filep.globalparam, &out_data_files, out_data);

  /** Check and Open Files **/
  check_files(&filep, &filenames);

  if (!options.OUTPUT_FORCE) {
    /** Read Vegetation Library File **/
    veg_lib = read_veglib(filep.veglib,&Nveg_type);
  } /* !OUTPUT_FORCE */

  /** Initialize Parameters **/
  cellnum = -1;

  /** Make Date Data Structure **/
  dmy      = make_dmy(&global_param);

  /** allocate memory for the atmos_data_struct **/
  alloc_atmos(global_param.nrecs, &atmos);

  /** Initial state **/
  startrec = 0;
  if (!options.OUTPUT_FORCE) {

    if ( options.INIT_STATE ) 
      filep.init_state = check_state_file(filenames.init_state, dmy, 
					   &global_param, options.Nlayer, 
					   options.Nnode, &startrec);

    /** open state file if model state is to be saved **/
    if ( options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 )
      filep.statefile = open_state_file(&global_param, filenames, options.Nlayer,
                                           options.Nnode);
    else filep.statefile = NULL;

  } /* !OUTPUT_FORCE */

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  while(!MODEL_DONE) {

    soil_con = read_soilparam(filep.soilparam, &RUN_MODEL, &MODEL_DONE);

    if(RUN_MODEL) {

      NEWCELL=TRUE;
      cellnum++;

      if (!options.OUTPUT_FORCE) {

        /** Read Grid Cell Vegetation Parameters **/
        veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
                                Nveg_type);
        calc_root_fractions(veg_con, &soil_con);

        if ( options.LAKES ) 
	  lake_con = read_lakeparam(filep.lakeparam, soil_con, veg_con);

      } /* !OUTPUT_FORCE */

      /** Build Gridded Filenames, and Open **/
      make_in_and_outfiles(&filep, &filenames, &soil_con, out_data_files);

      if (options.PRT_HEADER) {
        /** Write output file headers **/
        write_header(out_data_files, out_data, dmy, global_param);
      }

      if (!options.OUTPUT_FORCE) {

        /** Read Elevation Band Data if Used **/
        read_snowband(filep.snowband, &soil_con);

        /** Make Top-level Control Structure **/
        all_vars     = make_all_vars(veg_con[0].vegetat_type_num);

        /** allocate memory for the veg_hist_struct **/
        alloc_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);

      } /* !OUTPUT_FORCE */

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

      initialize_atmos(atmos, dmy, filep.forcing, veg_lib, veg_con, veg_hist,
		       &soil_con, out_data_files, out_data); 

      if (!options.OUTPUT_FORCE) {

        /**************************************************
          Initialize Energy Balance and Snow Variables 
        **************************************************/

#if VERBOSE
        fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
        ErrorFlag = initialize_model_state(&all_vars, dmy[0], &global_param, filep, 
			       soil_con.gridcel, veg_con[0].vegetat_type_num,
			       options.Nnode, 
			       atmos[0].air_temp[NR],
			       &soil_con, veg_con, lake_con);
        if ( ErrorFlag == ERROR ) {
	  if ( options.CONTINUEONERROR == TRUE ) {
	    // Handle grid cell solution error
	    fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	    break;
	  } else {
	    // Else exit program on cell solution error as in previous versions
	    sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
	    vicerror(ErrStr);
	  }
        }
      
#if VERBOSE
        fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

        /** Update Error Handling Structure **/
        Error.filep = filep;
        Error.out_data_files = out_data_files;

        /** Initialize the storage terms in the water and energy balances **/
        /** Sending a negative record number (-global_param.nrecs) to put_data() will accomplish this **/
	ErrorFlag = put_data(&all_vars, &atmos[0], &soil_con, veg_con, &lake_con, out_data_files, out_data, &save_data, &dmy[0], -global_param.nrecs);

        /******************************************
	  Run Model in Grid Cell for all Time Steps
	******************************************/

        for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

          if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
          else LASTREC = FALSE;

	  /**************************************************
	    Compute cell physics for 1 timestep
	  **************************************************/
	  ErrorFlag = full_energy(cellnum, rec, &atmos[rec], &all_vars, dmy, &global_param, &lake_con, &soil_con, veg_con, veg_hist);

	  /**************************************************
	    Write cell average values for current time step
	  **************************************************/
	  ErrorFlag = put_data(&all_vars, &atmos[rec], &soil_con, veg_con, &lake_con, out_data_files, out_data, &save_data, &dmy[rec], rec);

	  /************************************
	    Save model state at assigned date
	    (after the final time step of the assigned date)
	  ************************************/
	  if ( filep.statefile != NULL
	       &&  ( dmy[rec].year == global_param.stateyear
		     && dmy[rec].month == global_param.statemonth 
		     && dmy[rec].day == global_param.stateday
		     && ( rec+1 == global_param.nrecs
			  || dmy[rec+1].day != global_param.stateday ) ) )
	    write_model_state(&all_vars, &global_param, veg_con->vegetat_type_num, soil_con.gridcel, &filep, &soil_con, lake_con);


          if ( ErrorFlag == ERROR ) {
            if ( options.CONTINUEONERROR == TRUE ) {
              // Handle grid cell solution error
              fprintf(stderr, "ERROR: Grid cell %i failed in record %i so the simulation has not finished.  An incomplete output file has been generated, check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
              break;
            } else {
	      // Else exit program on cell solution error as in previous versions
              sprintf(ErrStr, "ERROR: Grid cell %i failed in record %i so the simulation has ended. Check your inputs before rerunning the simulation.\n", soil_con.gridcel, rec);
              vicerror(ErrStr);
	    }
          }

          NEWCELL=FALSE;

        } /* End Rec Loop */

      } /* !OUTPUT_FORCE */

      close_files(&filep,out_data_files,&filenames); 

      if (!options.OUTPUT_FORCE) {

        free_veg_hist(global_param.nrecs, veg_con[0].vegetat_type_num, &veg_hist);
        free_all_vars(&all_vars,veg_con[0].vegetat_type_num);
        free_vegcon(&veg_con);
        free((char *)soil_con.AreaFract);
        free((char *)soil_con.BandElev);
        free((char *)soil_con.Tfactor);
        free((char *)soil_con.Pfactor);
        free((char *)soil_con.AboveTreeLine);

      } /* !OUTPUT_FORCE */

    }	/* End Run Model Condition */
  } 	/* End Grid Loop */

  /** cleanup **/
  free_atmos(global_param.nrecs, &atmos);
  free_dmy(&dmy);
  free_out_data_files(&out_data_files);
  free_out_data(&out_data);
  fclose(filep.soilparam);
  if (!options.OUTPUT_FORCE) {
    free_veglib(&veg_lib);
    fclose(filep.vegparam);
    fclose(filep.veglib);
    if (options.SNOW_BAND>1)
      fclose(filep.snowband);
    if (options.LAKES)
      fclose(filep.lakeparam);
    if ( options.INIT_STATE )
      fclose(filep.init_state);
    if ( options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 )
      fclose(filep.statefile);
  } /* !OUTPUT_FORCE */

  return EXIT_SUCCESS;

}	/* End Main Program */
