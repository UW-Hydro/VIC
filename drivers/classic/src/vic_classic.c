#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

static char vcid[] = "$Id$";

// global variables

char *version = "4.2 beta 2014-Feb-25";
char *optstring = "g:vo";
int flag;
int NR;         /* array index for atmos struct that indicates
                   the model step avarage or sum */
int NF;         /* array index loop counter limit for atmos
                   struct that indicates the SNOW_STEP values */
 

global_param_struct global_param;
veg_lib_struct *veg_lib;
option_struct options;
Error_struct Error;
param_set_struct param_set;

/**************************************************************************
  Define some reference landcover types that always exist regardless
  of the contents of the library (mainly for potential evap calculations):
  Non-natural:
      satsoil = saturated bare soil
      h2osurf = open water surface (deep enough to have albedo of 0.08)
      short   = short reference crop (grass)
      tall    = tall reference crop (alfalfa)
  Natural:
      natveg  = current vegetation
      vegnocr = current vegetation with canopy resistance set to 0
  NOTE: bare soil roughness and displacement will be overwritten by the
        values found in the soil parameter file; bare soil wind_h will
        be overwritten by the value specified in the global param file.
**************************************************************************/

/* One element for each non-natural PET type */
char   ref_veg_over[]        = { 0, 0, 0, 0 };
double ref_veg_rarc[]        = { 0.0, 0.0, 25, 25 };
double ref_veg_rmin[]        = { 0.0, 0.0, 100, 100 };
double ref_veg_lai[]         = { 1.0, 1.0, 2.88, 4.45 };
double ref_veg_albedo[]      = { BARE_SOIL_ALBEDO, H2O_SURF_ALBEDO, 0.23, 0.23 };
double ref_veg_rough[]       = { 0.001, 0.001, 0.0148, 0.0615 };
double ref_veg_displ[]       = { 0.0054, 0.0054, 0.08, 0.3333 };
double ref_veg_wind_h[]      = { 10.0, 10.0, 10.0, 10.0 };
double ref_veg_RGL[]         = { 0.0, 0.0, 100, 100 };
double ref_veg_rad_atten[]   = { 0.0, 0.0, 0.0, 0.0 };
double ref_veg_wind_atten[]  = { 0.0, 0.0, 0.0, 0.0 };
double ref_veg_trunk_ratio[] = { 0.0, 0.0, 0.0, 0.0 };
/* One element for each PET type (non-natural or natural) */
char ref_veg_ref_crop[] = { FALSE, FALSE, TRUE, TRUE, FALSE, FALSE };

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
  01-Nov-04 Updated arglist for make_dist_prcp(), as part of fix for
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
  2008-May-05 Added prcp fraction (mu) to initial water storage
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
  2014-Mar-24 Removed ARC_SOIL option         BN
**********************************************************************/
{
 /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                     RUN_MODEL;
  char                    *init_STILL_STORM;
  char                     ErrStr[MAXSTRING];
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                     *init_DRY_TIME;
  int                      Ncells;
  int                      startrec;
  int                      ErrorFlag;
  float                    mu;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  double                   Clake;
  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_con_struct          *veg_con;
  soil_con_struct          soil_con;
  dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */
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
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
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

        /** Make Precipitation Distribution Control Structure **/
        prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);

      } /* !OUTPUT_FORCE */

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

      initialize_atmos(atmos, dmy, filep.forcing, &soil_con, out_data_files, out_data); 

      if (!options.OUTPUT_FORCE) {

        /**************************************************
          Initialize Energy Balance and Snow Variables 
        **************************************************/

#if VERBOSE
        fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
        ErrorFlag = initialize_model_state(&prcp, dmy[0], &global_param, filep, 
			       soil_con.gridcel, veg_con[0].vegetat_type_num,
			       options.Nnode, Ndist, 
			       atmos[0].air_temp[NR],
			       &soil_con, veg_con, lake_con,
			       &init_STILL_STORM, &init_DRY_TIME);
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
        /** Sending a negative record number (-global_param.nrecs) to dist_prec() will accomplish this **/
        ErrorFlag = dist_prec(&atmos[0], &prcp, &soil_con, veg_con,
		    &lake_con, dmy, &global_param, &filep, out_data_files,
		    out_data, &save_data, -global_param.nrecs, cellnum,
                    NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

        /******************************************
	  Run Model in Grid Cell for all Time Steps
	******************************************/

        for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

          if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
          else LASTREC = FALSE;

          ErrorFlag = dist_prec(&atmos[rec], &prcp, &soil_con, veg_con,
		  &lake_con, dmy, &global_param, &filep,
		  out_data_files, out_data, &save_data, rec, cellnum,
                  NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME);

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
	  for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ )
	    init_DRY_TIME[veg] = -999;

        } /* End Rec Loop */

      } /* !OUTPUT_FORCE */

      close_files(&filep,out_data_files,&filenames); 

      if (!options.OUTPUT_FORCE) {

        free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
        free_vegcon(&veg_con);
        free((char *)soil_con.AreaFract);
        free((char *)soil_con.BandElev);
        free((char *)soil_con.Tfactor);
        free((char *)soil_con.Pfactor);
        free((char *)soil_con.AboveTreeLine);
        free((char*)init_STILL_STORM);
        free((char*)init_DRY_TIME);

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
