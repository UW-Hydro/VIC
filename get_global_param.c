#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
global_param_struct get_global_param(filenames_struct *names,
                                     FILE             *gp)
/**********************************************************************
	global_param_struct	Dag Lohmann	January 1996

  This routine reads the global parameters, model options, and
  debugging controls from an input file.

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
  global_param_struct temp;
  
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
      else if(strcasecmp("ULAYER",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.Ulayer);
      }
      else if(strcasecmp("LLAYER",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&temp.Llayer);
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
        if(strcasecmp("TRUE",flgstr)==0) options.FROZEN_SOIL=TRUE;
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
      else if(strcasecmp("CALC_SNOW_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.CALC_SNOW_FLUX=TRUE;
        else options.CALC_SNOW_FLUX = FALSE;
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
      else if(strcasecmp("GRID_DECIMAL",optstr)==0) {
        sscanf(cmdstr,"%*s %i",&options.GRID_DECIMAL);
      }
      else if(strcasecmp("BINARY_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.BINARY_OUTPUT=TRUE;
        else options.BINARY_OUTPUT = FALSE;
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
      else if(strcasecmp("VEGPARAM",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veg);
      }
      else if(strcasecmp("VEGLIB",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veglib);
      }
      else if(strcasecmp("RESULT_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->result_dir);
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
      else if(strcasecmp("DEBUG_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",debug.debug_dir);
      }
    }
    fgets(cmdstr,MAXSTRING,gp);
  }

  return temp;
}
