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

  /* type 0: air temperature */
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

  /* type 1: albedo */
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

  /* type 2: air density */
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

  /* type 3: incomgin longwave radiation */
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

  /* type 4: precipitation */
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

  /* type 5: air pressure */
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

  /* type 6: shortwave radiation */
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

  /* type 7: maximum daily temperature */
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

  /* type 8: minimum daily temperature */
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

  /* type 9: sky cover */
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

  /* type 10: vapor pressure */
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

  /* type 11: wind speed */
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

  /* type 12: wind speed */
  else if(strcasecmp("WINDE",optstr)==0){
    type = WINDE;
    param_set.TYPE[type].SUPPLIED=file_num+1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;
    if(BINARY) {
      sscanf(cmdstr,"%*s %*s %s %lf",flgstr,
	     &param_set.TYPE[type].multiplier);
      if(strcasecmp("SIGNED",flgstr)==0) param_set.TYPE[type].SIGNED=TRUE;
      else param_set.TYPE[type].SIGNED=FALSE;
    }
  }

  /* type 13: unusead (blank) data */
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
