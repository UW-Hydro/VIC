#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id$";

/********************************************************************/
/*			GLOBAL VARIABLES                            */
/********************************************************************/
int NR;		      /* array index for atmos struct that indicates
			 the model step avarage or sum */
int NF;		      /* array index loop counter limit for atmos
			 struct that indicates the SNOW_STEP values */
 
global_param_struct get_global_param(filenames_struct *names,
                                     FILE             *gp)
/**********************************************************************
  get_global_param	Keith Cherkauer	            March 1998

  This routine reads the VIC model global control file, getting
  values for global parameters, model options, and debugging controls.

  NOTE: any additions or removals of parameters in this file must also
  be made in display_current_settings.c.

  Modifications:
  7-19-96 Modified to read time step		        KAC
  4-5-98  Modified to read model options and debugging
          controls from a single file                   KAC
  01-20-00 modified to work with new radiation estimation routines,
           new simplified frozen soil moisture, and new new open
           format forcing file rad routines.              KAC
  02-27-01 added reads for lake model parameters          KAC
  04-21-03 added parameters for blowing snow algorithm, printing
           lake variables during debugging and reading Bart's 
           new Arno parameters.                           KAC
  11-May-04 Modified to display compile-time and run-time options
	    if VERBOSE is set to TRUE.			TJB
  13-Oct-04 Added validation for GRND_FLUX option.              TJB
  01-Nov-04 Added validation for Nnodes with QUICK_FLUX option, as
	    part of fix for QUICK_FLUX state file compatibility.TJB
  2005-03-08 Added EQUAL_AREA option.				TJB
  2005-03-24 Added ALMA_OUTPUT option.				TJB
  2005-04-07 Fixed state file warning check.			TJB
  2005-Apr-13 Added logic for OUTPUT_FORCE option.		TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW	TJB

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
#if LINK_DEBUG
  extern debug_struct     debug;
#endif
  extern int              NF, NR;

  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  char flgstr[MAXSTRING];
  char ErrStr[MAXSTRING];
  int  file_num;
  int  field;
  int  i;
  global_param_struct global;

  /** Initialize non-global parameters **/
  global.endmonth      = MISSING;
  global.endday        = MISSING;
  global.endyear       = MISSING;
  global.skipyear      = 0;
  for(i = 0; i < 2; i++) {
    global.forcemonth[i] = 1;
    global.forceday[i]   = 1;
    global.forceyear[i]  = MISSING;
    global.forcehour[i]  = 0;
    global.forceskip[i]  = 0;
    strcpy(names->forcing[i],"FALSE");
  }
  file_num             = 0;
  global.nrecs         = MISSING;
#if SAVE_STATE
  global.stateyear = MISSING;
  strcpy(global.statename, "NONE");
#endif

  /** Read through global control file to find parameters **/

  fgets(cmdstr,MAXSTRING,gp);

  while(!feof(gp)) {
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {

      sscanf(cmdstr,"%s",optstr);

      /*******************************
        Get Model Global Parameters
	*****************************/
      if(strcasecmp("NLAYER",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&options.Nlayer);
      }
      else if(strcasecmp("TIME_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.dt);
      }
      else if(strcasecmp("RESOLUTION",optstr)==0) {
        sscanf(cmdstr,"%*s %f",&global.resolution);
      }
      else if(strcasecmp("STARTYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.startyear);
      }
      else if(strcasecmp("STARTMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.startmonth);
      }
      else if(strcasecmp("STARTDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.startday);
      }
      else if(strcasecmp("STARTHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.starthour);
      }
      else if(strcasecmp("ENDYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.endyear);
      }
      else if(strcasecmp("ENDMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.endmonth);
      }
      else if(strcasecmp("ENDDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.endday);
      }
      else if(strcasecmp("SKIPYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.skipyear);
      }
      else if(strcasecmp("FORCEYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.forceyear[file_num]);
      }
      else if(strcasecmp("FORCEMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.forcemonth[file_num]);
      }
      else if(strcasecmp("FORCEDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.forceday[file_num]);
      }
      else if(strcasecmp("FORCEHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.forcehour[file_num]);
      }
      else if(strcasecmp("NRECS",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.nrecs);
      }
      else if(strcasecmp("WIND_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.wind_h);
      }
      else if(strcasecmp("MEASURE_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.measure_h);
      }
      else if(strcasecmp("NODES",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&options.Nnode);
      }
      else if(strcasecmp("MIN_RAIN_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.MIN_RAIN_TEMP);
      }
      else if(strcasecmp("MAX_SNOW_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global.MAX_SNOW_TEMP);
      }

      /********************
        Get Model Options
	******************/
      else if(strcasecmp("FULL_ENERGY",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  options.FULL_ENERGY=TRUE;
 	  options.QUICK_FLUX=TRUE;
	  options.GRND_FLUX=TRUE;
	}
	else options.FULL_ENERGY = FALSE;
      }
      else if(strcasecmp("FROZEN_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  options.FROZEN_SOIL=TRUE;
	  options.QUICK_FLUX=FALSE;
	  options.GRND_FLUX=TRUE;
	}
        else options.FROZEN_SOIL = FALSE;
      }
      else if(strcasecmp("QUICK_SOLVE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.QUICK_SOLVE=TRUE;
        else options.QUICK_SOLVE = FALSE;
      }
      else if(strcasecmp("NOFLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.NOFLUX=TRUE;
        else options.NOFLUX = FALSE;
      }
      else if(strcasecmp("DIST_PRCP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.DIST_PRCP=TRUE;
        else options.DIST_PRCP = FALSE;
      }
      else if(strcasecmp("BLOWING",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.BLOWING=TRUE;
        else options.BLOWING = FALSE;
      }
      else if(strcasecmp("COMPRESS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.COMPRESS=TRUE;
        else options.COMPRESS = FALSE;
      }
      else if(strcasecmp("CORRPREC",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.CORRPREC=TRUE;
        else options.CORRPREC = FALSE;
      }
      else if(strcasecmp("PRT_SNOW_BAND",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.PRT_SNOW_BAND=TRUE;
        else options.PRT_SNOW_BAND = FALSE;
      }
      else if(strcasecmp("GRID_DECIMAL",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&options.GRID_DECIMAL);
      }
      else if(strcasecmp("SNOW_BAND",optstr)==0) {
	sscanf(cmdstr,"%*s %i %s",&options.SNOW_BAND,names->snow_band);
      }
      else if(strcasecmp("BINARY_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.BINARY_OUTPUT=TRUE;
        else options.BINARY_OUTPUT = FALSE;
      }
      else if(strcasecmp("EQUAL_AREA",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.EQUAL_AREA=TRUE;
        else options.EQUAL_AREA = FALSE;
      }
      else if(strcasecmp("ALMA_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ALMA_OUTPUT=TRUE;
        else options.ALMA_OUTPUT = FALSE;
      }
      else if(strcasecmp("ARC_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ARC_SOIL=TRUE;
        else options.ARC_SOIL = FALSE;
      }
      else if(strcasecmp("SNOW_STEP",optstr)==0) {
	sscanf(cmdstr,"%*s %i",&options.SNOW_STEP);
      }
      else if(strcasecmp("ROOT_ZONES",optstr)==0) {
	sscanf(cmdstr,"%*s %i",&options.ROOT_ZONES);
      }
      else if(strcasecmp("PREC_EXPT",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&options.PREC_EXPT);
      }
      else if(strcasecmp("MIN_WIND_SPEED",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&options.MIN_WIND_SPEED);
      }
      else if(strcasecmp("GRND_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.GRND_FLUX=TRUE;
        else options.GRND_FLUX = FALSE;
      }
      else if(strcasecmp("QUICK_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.QUICK_FLUX=TRUE;
        else options.QUICK_FLUX = FALSE;
      }
      else if(strcasecmp("MOISTFRACT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.MOISTFRACT=TRUE;
        else options.MOISTFRACT = FALSE;
      }
      else if(strcasecmp("INIT_STATE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) options.INIT_STATE=FALSE;
        else {
	  options.INIT_STATE = TRUE;
	  strcpy(names->init_state,flgstr);
	}
      }
      else if(strcasecmp("BINARY_STATE_FILE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) options.BINARY_STATE_FILE=FALSE;
	else options.BINARY_STATE_FILE=TRUE;
      }
      else if (strcasecmp("ARNO_PARAMS", optstr)==0) {
	fprintf(stderr,"WARNING: Please change \"ARNO_PARAMS\" to \"NIJSSEN2001_BASEFLOW\" in your global parameter file\n");
	sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.NIJSSEN2001_BASEFLOW=TRUE;
        else options.NIJSSEN2001_BASEFLOW = FALSE;
      }
      else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr)==0) {
	sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.NIJSSEN2001_BASEFLOW=TRUE;
        else options.NIJSSEN2001_BASEFLOW = FALSE;
      }
#if SAVE_STATE
      else if(strcasecmp("STATENAME",optstr)==0) {
        sscanf(cmdstr,"%*s %s",global.statename);
      }
      else if(strcasecmp("STATEYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.stateyear);
      }
      else if(strcasecmp("STATEMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.statemonth);
      }
      else if(strcasecmp("STATEDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&global.stateday);
      }
#endif // SAVE_STATE
#if LAKE_MODEL
      else if(strcasecmp("LAKES",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) options.LAKES = FALSE;
        else {
	  options.LAKES = TRUE;
	  strcpy(names->lakeparam, flgstr);
	}
      }
      else if(strcasecmp("LAKE_PROFILE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) options.LAKE_PROFILE = FALSE;
        else {
	  options.LAKE_PROFILE = TRUE;
	}
      }
#endif /* LAKE_MODEL */


      /************************************
        Get Forcing Data File Information
	**********************************/
      else if(strcasecmp("FORCING1",optstr)==0) {
	if ( strcmp( names->forcing[0], "FALSE" ) != 0 ) 
	  nrerror("Tried to define FORCING1 twice, if you want to use two forcing files, the second must be defined as FORCING2");
        sscanf(cmdstr,"%*s %s", names->forcing[0]);
	file_num = 0;
	field=0;
      }
      else if(strcasecmp("FORCING2",optstr)==0) {
        sscanf(cmdstr,"%*s %s", names->forcing[1]);
	file_num = 1;
	field=0;
      }
      else if(strcasecmp("N_TYPES",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&param_set.N_TYPES[file_num]);
      }
      else if(strcasecmp("FORCE_TYPE",optstr)==0) {
	get_force_type(cmdstr,file_num,&field);
      }
      else if(strcasecmp("FORCE_DT",optstr)==0) {
	sscanf(cmdstr,"%*s %i ", &param_set.FORCE_DT[file_num]);
      }
      else if (strcasecmp("FORCE_ENDIAN",optstr)==0) {
	sscanf(cmdstr, "%*s %s", flgstr);
	if (strcasecmp(flgstr, "LITTLE") == 0)
	  param_set.FORCE_ENDIAN[file_num] = LITTLE;
	else if (strcasecmp(flgstr, "BIG") == 0)
	  param_set.FORCE_ENDIAN[file_num] = BIG;
	else
	  nrerror("FORCE_ENDIAN must be either BIG or LITTLE.");
      }
      else if (strcasecmp("FORCE_FORMAT",optstr)==0) {
	sscanf(cmdstr, "%*s %s", flgstr);
	if (strcasecmp(flgstr, "BINARY") == 0)
	  param_set.FORCE_FORMAT[file_num] = BINARY;
	else if (strcasecmp(flgstr, "ASCII") == 0)
	  param_set.FORCE_FORMAT[file_num] = ASCII;
	else
	  nrerror("FORCE_FORMAT must be either ASCII or BINARY.");
      }

      /************************************
	Get Information for Parameter Files 
	************************************/
      
      else if(strcasecmp("SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil);
      }
      else if(strcasecmp("SOIL_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil_dir);
      }
      else if(strcasecmp("VEGPARAM",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veg);
      }
      else if(strcasecmp("VEGLIB",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veglib);
      }
      else if(strcasecmp("GLOBAL_LAI",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.GLOBAL_LAI=TRUE;
        else options.GLOBAL_LAI = FALSE;
      }
      else if(strcasecmp("RESULT_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->result_dir);
      }

      /******************************
        Get Model Debugging Options
	****************************/

#if LINK_DEBUG
      else if(strcasecmp("PRT_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_FLUX=TRUE;
        else debug.PRT_FLUX = FALSE;
      }
      else if(strcasecmp("PRT_BALANCE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_BALANCE=TRUE;
        else debug.PRT_BALANCE = FALSE;
      }
      else if(strcasecmp("PRT_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SOIL=TRUE;
        else debug.PRT_SOIL = FALSE;
      }
      else if(strcasecmp("PRT_VEGE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_VEGE=TRUE;
        else debug.PRT_VEGE = FALSE;
      }
      else if(strcasecmp("PRT_GLOBAL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_GLOBAL=TRUE;
        else debug.PRT_GLOBAL = FALSE;
      }
      else if(strcasecmp("PRT_ATMOS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_ATMOS=TRUE;
        else debug.PRT_ATMOS = FALSE;
      }
      else if(strcasecmp("PRT_SNOW",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SNOW=TRUE;
        else debug.PRT_SNOW = FALSE;
      }
      else if(strcasecmp("PRT_MOIST",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_MOIST=TRUE;
        else debug.PRT_MOIST = FALSE;
      }
      else if(strcasecmp("PRT_LAKE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_LAKE=TRUE;
        else debug.PRT_LAKE = FALSE;
      }
      else if(strcasecmp("PRT_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_TEMP=TRUE;
        else debug.PRT_TEMP = FALSE;
      }
      else if(strcasecmp("DEBUG_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",debug.debug_dir);
      }
#endif

      /***********************************
        Unrecognized Global Parameter Flag
        ***********************************/
      else {
	fprintf(stderr,"WARNING: Unrecognized option in the global parameter file:\n\t%s is unknown - check your spelling\n", optstr);
      }
    }
    fgets(cmdstr,MAXSTRING,gp);
  }

  /******************************************
    Check for undefined required parameters
  ******************************************/

  if ( strcmp ( names->forcing[0], "FALSE" ) == 0 )
    nrerror("No forcing file has been defined, make sure that the global files defines FORCING1.");

  for ( i = 0; i < 2; i++ ) {
    if ( i == 0 || (i == 1 && param_set.N_TYPES[i] != MISSING) ) {
      if (param_set.N_TYPES[i] == MISSING) {
	sprintf(ErrStr,"Need to specify the number forcing variables types in forcing file %i.", i);
	nrerror(ErrStr);
      }
      if (param_set.FORCE_FORMAT[i] == MISSING) {
	sprintf(ErrStr,"Need to specify the INPUT_FORMAT (ASCII or BINARY) for forcing file %i.",i);
	nrerror(ErrStr);
      }
      if (param_set.FORCE_INDEX[i][param_set.N_TYPES[i]-1] == MISSING) {
	sprintf(ErrStr,"Did not define enough forcing variables in forcing file %i.",i);
	nrerror(ErrStr);
      }
      if(param_set.FORCE_DT[i] == MISSING ) {
	sprintf(ErrStr,"Must define time steps (FORCE_DT <dt>) in control file for focing file %i.",file_num);
	nrerror(ErrStr);
      }
    }
  }

  if(param_set.N_TYPES[1] != MISSING && global.forceyear[1] == MISSING) {
    global.forceyear[1] = global.forceyear[0];
    global.forcemonth[1] = global.forcemonth[0];
    global.forceday[1] = global.forceday[0];
    global.forcehour[1] = global.forcehour[0];
    global.forceskip[1] = 0;
  }

#if !OUTPUT_FORCE

  if(options.ROOT_ZONES<0)
    nrerror("ROOT_ZONES must be defined to a positive integer greater than 0, in the global control file.");
  if(options.Nlayer > MAX_LAYERS) {
    sprintf(ErrStr,"Global file wants more soil moisture layers (%i) than are defined by MAX_LAYERS (%i).  Edit user_def.h and recompile.",options.Nlayer,MAX_LAYERS);
    nrerror(ErrStr);
  }
  if(options.Nnode > MAX_NODES) {
    sprintf(ErrStr,"Global file wants more soil thermal nodes (%i) than are defined by MAX_NODES (%i).  Edit user_def.h and recompile.",options.Nnode,MAX_NODES);
    nrerror(ErrStr);
  }
#if LAKE_MODEL
  if ( options.Nlakenode > MAX_LAKE_NODES) {
    sprintf(ErrStr,"Global file wants more lake thermal nodes (%i) than are defined by MAX_LAKE_NODES (%i).  Edit user_def.h and recompile.", options.Nlakenode, MAX_LAKE_NODES);
    nrerror(ErrStr);
  }
#endif // LAKE_MODEL 
  if(options.QUICK_FLUX) {
    if((options.FULL_ENERGY || options.FROZEN_SOIL) && options.Nnode != 3) {
      sprintf(ErrStr,"To run the model in FULL_ENERGY or FROZEN_SOIL modes with QUICK_FLUX=TRUE, you must define exactly 3 soil thermal nodes.  Currently Nnodes is set to  %i.",options.Nnode);
      nrerror(ErrStr);
    }
    else if (!options.FULL_ENERGY && !options.FROZEN_SOIL && options.Nnode != 1) {
      sprintf(ErrStr,"To run the model with FULL_ENERGY=FALSE, FROZEN_SOIL=FALSE, and QUICK_FLUX=TRUE, you must define exactly 1 soil thermal node.  Currently Nnodes is set to  %i.",options.Nnode);
      nrerror(ErrStr);
    }
  }
  else {
    if(options.Nnode < 4) {
      sprintf(ErrStr,"To run the model with QUICK_FLUX=FALSE, you must define at least 4 soil thermal nodes.  Currently Nnodes is set to  %i.",options.Nnode);
      nrerror(ErrStr);
    }
  }
  if((options.FULL_ENERGY || options.FROZEN_SOIL) && options.Nlayer<3) {
    sprintf(ErrStr,"You must define at least 3 soil moisture layers to run the model in FULL_ENERGY or FROZEN_SOIL modes.  Currently Nlayers is set to  %i.",options.Nlayer);
    nrerror(ErrStr);
  }
  if(!options.FULL_ENERGY && !options.FROZEN_SOIL && options.GRND_FLUX) {
    sprintf(ErrStr,"Both FULL_ENERGY and FROZEN_SOIL are FALSE, but GRND_FLUX is TRUE.\nThis combination of options is not recommended.  Unless you intend\nto use this combination of options, we recommend commenting out\nthe GRND_FLUX entry in your global file.  To do this,place a \"#\"\nat the beginning of the line containing \"GRND_FLUX\".\n");
    nrerror(ErrStr);
  }
  if(options.FULL_ENERGY && !options.GRND_FLUX) {
    sprintf(ErrStr,"FULL_ENERGY is TRUE, but GRND_FLUX is explicitly set to FALSE.\nThis combination of options is not recommended.  Unless you intend\nto use this combination of options, we recommend commenting out\nthe GRND_FLUX entry in your global file.  To do this,place a \"#\"\nat the beginning of the line containing \"GRND_FLUX\".\n");
    nrerror(ErrStr);
  }
  if(options.FROZEN_SOIL && !options.GRND_FLUX) {
    sprintf(ErrStr,"FROZEN_SOIL is TRUE, but GRND_FLUX is explicitly set to FALSE.\nThis combination of options is not recommended.  Unless you intend\nto use this combination of options, we recommend commenting out\nthe GRND_FLUX entry in your global file.  To do this,place a \"#\"\nat the beginning of the line containing \"GRND_FLUX\".\n");
    nrerror(ErrStr);
  }
  if(options.SNOW_BAND > MAX_BANDS) {
    sprintf(ErrStr,"Global file wants more snow bands (%i) than are defined by MAX_BANDS (%i).  Edit user_def.h and recompile.",options.SNOW_BAND,MAX_BANDS);
    nrerror(ErrStr);
  }
#if SAVE_STATE
  if ( strcmp( names->init_state, global.statename ) == 0 && strcmp( names->init_state, "" ) != 0 ) {
    sprintf(ErrStr,"The save state file (%s) has the same name as the initialize state file (%s).  The initialize state file will be destroyed when the save state file is opened.", global.statename, names->init_state);
    nrerror(ErrStr);
  }
  if ( strcmp( global.statename, "" ) == 0 ) {
    fprintf(stderr,"WARNING: Model compiled with SAVE_STATE activated, but no state files were defined.\n");
  }
#endif // SAVE_STATE
#if LAKE_MODEL
  if ( global.resolution == 0 && options.LAKES ) {
    sprintf(ErrStr, "The model grid cell resolution (RESOLUTION) must be defined in the global control file when the lake model is active.");
    nrerror(ErrStr);
  }
  if ( options.LAKES && global.resolution > 360 && !options.EQUAL_AREA ) {
    sprintf(ErrStr, "For EQUAL_AREA=FALSE, the model grid cell resolution (RESOLUTION) must be set to the number of lat or lon degrees per grid cell.  This cannot exceed 360.");
    nrerror(ErrStr);
  }
#endif  // LAKE_MODEL

#endif  // !OUTPUT_FORCE

#if OUTPUT_FORCE
  options.SNOW_STEP = global.dt;
#endif  // OUTPUT_FORCE

  /* set NR and NF */
  if (global.dt < 24 && global.dt != options.SNOW_STEP)
    nrerror("If the model step is smaller than daily, the snow model should run\nat the same time step as the rest of the model.");

  NF = global.dt/options.SNOW_STEP;
  if (global.dt % options.SNOW_STEP != 0 || options.SNOW_STEP > global.dt)
    nrerror("SNOW_STEP should be smaller than TIME_STEP and divide TIME_STEP evenly ");
  if (NF == 1)
    NR = 0;
  else
    NR = NF;

#if !OUTPUT_FORCE

  /*********************************
    Output major options to stderr
  *********************************/
#if VERBOSE
  display_current_settings(DISP_ALL,names,&global);
#else
  display_current_settings(DISP_VERSION,names,&global);
#endif

#if VERBOSE
  fprintf(stderr,"Time Step = %i hour(s)\n",global.dt);
  fprintf(stderr,"Simulation start date = %02i/%02i/%04i\n",
	  global.startday, global.startmonth, global.startyear);
  if ( global.nrecs > 0 )
    fprintf(stderr,"Number of Records = %i\n\n",global.nrecs);
  else 
    fprintf(stderr,"Simulation end date = %02i/%02i/%04i\n\n",
	    global.endday, global.endmonth, global.endyear);
  fprintf(stderr,"Full Energy...................(%i)\n",options.FULL_ENERGY);
  fprintf(stderr,"Use Distributed Precipitation.(%i)\n",options.DIST_PRCP);
  if(options.DIST_PRCP)
    fprintf(stderr,"..Using Precipitation Exponent of %f\n",options.PREC_EXPT);
  if ( options.GRND_FLUX ) {
    fprintf(stderr,"Ground heat flux will be estimated ");
    if ( options.QUICK_FLUX ) 
      fprintf(stderr,"using Liang, Wood and Lettenmaier (1999).\n");
    else 
      fprintf(stderr,"using Cherkauer and Lettenmaier (1999).\n");
  }
  else
    fprintf(stderr,"Ground heat flux not computed (no energy balance).\n");
  fprintf(stderr,"Use Frozen Soil Model.........(%i)\n",options.FROZEN_SOIL);
  if ( QUICK_FS )
    fprintf(stderr,".... Using linearized UFWC curve with %i temperatures.\n",
	    QUICK_FS_TEMPS);
  fprintf(stderr,"Run Snow Model Using a Time Step of %i hours\n", 
	  options.SNOW_STEP);
  fprintf(stderr,"Compress Output Files.........(%i)\n",options.COMPRESS);
  fprintf(stderr,"Correct Precipitation.........(%i)\n",options.CORRPREC);
  fprintf(stderr,"\n");
  fprintf(stderr,"Using %i Snow Bands\n",options.SNOW_BAND);
  fprintf(stderr,"Using %i Root Zones\n",options.ROOT_ZONES);
#if SAVE_STATE
  if ( global.stateyear != MISSING )
    fprintf(stderr,"Model state will be saved on = %02i/%02i/%04i\n\n",
	    global.stateday, global.statemonth, global.stateyear);
#endif
  if ( OPTIMIZE )
    fprintf(stderr,"Model is using optimized output (runoff and baseflow only).\n");
  else if ( LDAS_OUTPUT )
    fprintf(stderr,"Model output is in LDAS binary short int format.\n");
  else if ( options.BINARY_OUTPUT ) 
    fprintf(stderr,"Model output is in standard BINARY format.\n");
  else 
    fprintf(stderr,"Model output is in standard ASCII format.\n");
  if ( LINK_DEBUG ) 
    fprintf(stderr,"Debugging code has been included in the executable.\n");
  else 
    fprintf(stderr,"Debugging code has not been compiled.\n");
#endif

#endif // !OUTPUT_FORCE

  return global;

}
