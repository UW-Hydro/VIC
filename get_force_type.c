#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

void get_force_type(char   *cmdstr, 
		    int     file_num,
		    int    *field) {
/*************************************************************
  get_force_type.c      Keith Cherkauer     January 20, 2000

  This routine determines the current forcing file data type
  and stores its location in the description of the current 
  forcing file.

  2006-Sep-14 (Port from 4.1.0) Added support for ALMA input variables. TJB

*************************************************************/

  extern param_set_struct param_set;

  char optstr[50];
  char flgstr[10];
  char ErrStr[MAXSTRING];
  int  type;

  if((*field) >= param_set.N_TYPES[file_num]) {
    sprintf(ErrStr,"Too many variables defined for forcing file %i.",file_num);
    nrerror(ErrStr);
  }

  sscanf(cmdstr,"%*s %s",optstr);

  /***************************************
    Get meteorological data forcing info
  ***************************************/

  /* type 0: air temperature [C] */
  if(strcasecmp("AIR_TEMP",optstr)==0){
    type = AIR_TEMP;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 1: albedo [fraction] */
  else if(strcasecmp("ALBEDO",optstr)==0){
    type = ALBEDO;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 2: convective rainfall [mm/s] */
  else if(strcasecmp("CRAINF",optstr)==0){
    type = CRAINF;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 3: air density [kg/m3] */
  else if(strcasecmp("DENSITY",optstr)==0){
    type = DENSITY;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 4: incoming longwave radiation [W/m2] */
  else if(strcasecmp("LONGWAVE",optstr)==0){
    type = LONGWAVE;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 5: large-scale rainfall [mm/s] */
  else if(strcasecmp("LSRAINF",optstr)==0){
    type = LSRAINF;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 6: precipitation [mm] */
  else if(strcasecmp("PREC",optstr)==0){
    type = PREC;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 7: air pressure [kPa] */
  else if(strcasecmp("PRESSURE",optstr)==0){
    type = PRESSURE;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 8: air pressure [Pa] */
  else if(strcasecmp("PSURF",optstr)==0){
    type = PSURF;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 9: specific humidity [kg/kg] */
  else if(strcasecmp("QAIR",optstr)==0){
    type = QAIR;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 10: rainfall rate [mm/s] */
  else if(strcasecmp("RAINF",optstr)==0){
    type = RAINF;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 11: shortwave radiation [W/m2] */
  else if(strcasecmp("SHORTWAVE",optstr)==0){
    type = SHORTWAVE;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 12: snowfall rate [mm/s] */
  else if(strcasecmp("SNOWF",optstr)==0){
    type = SNOWF;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 13: surface air temperature [K] */
  else if(strcasecmp("TAIR",optstr)==0){
    type = TAIR;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 14: maximum daily temperature [C] */
  else if(strcasecmp("TMAX",optstr)==0){
    type = TMAX;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 15: minimum daily temperature [C] */
  else if(strcasecmp("TMIN",optstr)==0){
    type = TMIN;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 16: sky cover */
  else if(strcasecmp("TSKC",optstr)==0){
    type = TSKC;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 17: vapor pressure [kPa] */
  else if(strcasecmp("VP",optstr)==0){
    type = VP;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 18: wind speed [m/s] */
  else if(strcasecmp("WIND",optstr)==0){
    type = WIND;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 19: zonal wind speed [m/s] */
  else if(strcasecmp("WIND_E",optstr)==0){
    type = WIND_E;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 20: meridional wind speed [m/s] */
  else if(strcasecmp("WIND_N",optstr)==0){
    type = WIND_N;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 21: unused (blank) data */
  else if(strcasecmp("SKIP",optstr)==0){
    type = SKIP;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      param_set.TYPE[type].multiplier = 1;
      param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /** Undefined variable type **/
  else {
    sprintf(ErrStr,"Undefined forcing variable type %s in file %i.",
	    optstr, file_num);
    nrerror(ErrStr);
  }

  (*field)++;

}
