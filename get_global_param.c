#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id$";

global_param_struct get_global_param(filenames_struct *names,
                                     FILE             *gp,
				     int              *force_dt)
/**********************************************************************
  get_global_param	Keith Cherkauer	            March 1998

  This routine reads the VIC model global control file, getting
  values for global parameters, model options, and debugging controls.

  Modifications:
  7-19-96 Modified to read time step		        KAC
  4-5-98  Modified to read model options and debugging
          controls from a single file                   KAC
**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;
  
  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  char flgstr[MAXSTRING];
  char ErrStr[MAXSTRING];
  global_param_struct temp;

  /** Initialize non-global parameters **/
  temp.Nnodes = 3;

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
        sscanf(cmdstr,"%*s %i",&temp.dt);
      }
      else if(strcasecmp("RESOLUTION",optstr)==0) {
        sscanf(cmdstr,"%*s %f",&temp.resolution);
      }
      else if(strcasecmp("STARTYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.startyear);
      }
      else if(strcasecmp("STARTMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.startmonth);
      }
      else if(strcasecmp("STARTDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.startday);
      }
      else if(strcasecmp("STARTHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.starthour);
      }
      else if(strcasecmp("ENDYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.endyear);
      }
      else if(strcasecmp("ENDMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.endmonth);
      }
      else if(strcasecmp("ENDDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.endday);
      }
      else if(strcasecmp("NRECS",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.nrecs);
      }
      else if(strcasecmp("WIND_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&temp.wind_h);
      }
      else if(strcasecmp("MEASURE_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&temp.measure_h);
      }
      else if(strcasecmp("NODES",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.Nnodes);
      }
      else if(strcasecmp("MIN_RAIN_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&temp.MIN_RAIN_TEMP);
      }
      else if(strcasecmp("MAX_SNOW_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&temp.MAX_SNOW_TEMP);
      }

      /********************
        Get Model Options
	******************/
      else if(strcasecmp("FULL_ENERGY",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.FULL_ENERGY=TRUE;
        else options.FULL_ENERGY = FALSE;
      }
      else if(strcasecmp("FROZEN_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  options.FROZEN_SOIL=TRUE;
	  options.FS_FLUXES=TRUE;
	}
        else options.FROZEN_SOIL = FALSE;
      }
      else if(strcasecmp("SNOW_MODEL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.SNOW_MODEL=TRUE;
        else options.SNOW_MODEL = FALSE;
      }
      else if(strcasecmp("DIST_PRCP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.DIST_PRCP=TRUE;
        else options.DIST_PRCP = FALSE;
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
      else if(strcasecmp("ARC_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ARC_SOIL=TRUE;
        else options.ARC_SOIL = FALSE;
      }
      else if(strcasecmp("SNOW_STEP",optstr)==0) {
	sscanf(cmdstr,"%*s %i",&options.SNOW_STEP);
      }
      else if(strcasecmp("FORCE_DT",optstr)==0) {
	sscanf(cmdstr,"%*s %i %i",&force_dt[0],&force_dt[1]);
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

      /************************************
        Get Frocing Data File Information
	**********************************/
      else if(strcasecmp("FORCE_TYPE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",options.FORCE_TYPE);
      }
      else if(strcasecmp("FORCING1",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->forcing[0]);
      }
      else if(strcasecmp("FORCING2",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->forcing[1]);
      }
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
      else if(strcasecmp("RESULT_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->result_dir);
      }
      else if(strcasecmp("INIT_SNOW",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->init_snow);
	options.INIT_SNOW = TRUE;
      }
      else if(strcasecmp("INIT_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->init_soil);
	options.INIT_SOIL = TRUE;
      }


      /******************************
        Get Model Debugging Options
	****************************/
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
      else if(strcasecmp("PRT_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_TEMP=TRUE;
        else debug.PRT_TEMP = FALSE;
      }
      else if(strcasecmp("DEBUG_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",debug.debug_dir);
      }

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
/*   if(!options.FULL_ENERGY && options.FROZEN_SOIL)  */
/*     options.FULL_ENERGY = TRUE; */
  if((!options.FULL_ENERGY && !options.FROZEN_SOIL) && options.CALC_SNOW_FLUX) 
    options.CALC_SNOW_FLUX = FALSE;
  if(force_dt[0]<0 || force_dt[1]<0)
    nrerror("Must define forcing file time steps (FORCE_DT <dt 1> <dt 2>) in control file");
  if(options.ROOT_ZONES<0)
    nrerror("ROOT_ZONES must be defined to a positive integer greater than 0, in the global control file.");
  if(options.Nlayer > MAX_LAYERS) {
    sprintf(ErrStr,"Global file wants more soil moisture layers (%i) than are defined by MAX_LAYERS (%i).  Edit vicNl_def.h and recompile.",options.Nlayer,MAX_LAYERS);
    nrerror(ErrStr);
  }
  if(temp.Nnodes > MAX_NODES) {
    sprintf(ErrStr,"Global file wants more soil thermal nodes (%i) than are defined by MAX_NODES (%i).  Edit vicNl_def.h and recompile.",temp.Nnodes,MAX_NODES);
    nrerror(ErrStr);
  }
  if(options.SNOW_BAND > MAX_BANDS) {
    sprintf(ErrStr,"Global file wants more snow bands (%i) than are defined by MAX_BANDS (%i).  Edit vicNl_def.h and recompile.",options.SNOW_BAND,MAX_BANDS);
    nrerror(ErrStr);
  }

  /*********************************
    Output major options to stderr
  *********************************/
  fprintf(stderr,"Time Step = %i hour(s)\n",temp.dt);
  fprintf(stderr,"Number of Records = %i\n\n",temp.nrecs);
  fprintf(stderr,"Full Energy...................(%i)\n",options.FULL_ENERGY);
  fprintf(stderr,"Use Distributed Precipitation.(%i)\n",options.DIST_PRCP);
  if(options.DIST_PRCP)
    fprintf(stderr,"..Using Precipitation Exponent of %lf\n",options.PREC_EXPT);
  fprintf(stderr,"Use Frozen Soil Model.........(%i)\n",options.FROZEN_SOIL);
  fprintf(stderr,"Run Snow Model................(%i)\n",options.SNOW_MODEL);
  if(options.SNOW_MODEL && !options.FULL_ENERGY)
    fprintf(stderr,"..Using Time Step of %i hours\n",options.SNOW_STEP);
  fprintf(stderr,"Compress Output Files.........(%i)\n",options.COMPRESS);
  fprintf(stderr,"Correct Precipitation.........(%i)\n",options.CORRPREC);

  fprintf(stderr,"\n");
  fprintf(stderr,"Using %i Snow Bands\n",options.SNOW_BAND);
  fprintf(stderr,"Using %i Root Zones\n",options.ROOT_ZONES);

  return temp;

}
