#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

filenames_struct make_in_and_outfiles(infiles_struct   *infp, 
				      filenames_struct *filenames,
				      soil_con_struct  *soil,
				      outfiles_struct  *outfp)
/**********************************************************************
	make_in_and_outfile	Dag Lohman	January 1996

  This program builds the files names for input and output of grided
  data files.

  Modifications:
  5/20/96	The routine was modified to accept a variable
		number of layers, as well as to work with 
		frozen soils					KAC

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern FILE *open_file(char string[], char type[]);

  char             latchar[10], lngchar[10], junk[5];
  filenames_struct fnames;

  fnames = *filenames;

  sprintf(junk, "%%.%if", options.GRID_DECIMAL);
  sprintf(latchar, junk, soil->lat);
  sprintf(lngchar, junk, soil->lng);
 
  strcat(fnames.forcing[0], latchar);
  strcat(fnames.forcing[0], "_");
  strcat(fnames.forcing[0], lngchar);
  if(param_set.FORCE_FORMAT[0] == BINARY)
    infp->forcing[0] = open_file(fnames.forcing[0], "rb");
  else
    infp->forcing[0] = open_file(fnames.forcing[0], "r");

  infp->forcing[1] = NULL;
  if(strcasecmp(fnames.forcing[1],"FALSE")!=0) {
    strcat(fnames.forcing[1], latchar);
    strcat(fnames.forcing[1], "_");
    strcat(fnames.forcing[1], lngchar);
    if(param_set.FORCE_FORMAT[0] == BINARY) 
      infp->forcing[1] = open_file(fnames.forcing[1], "rb");
    else 
      infp->forcing[1] = open_file(fnames.forcing[1], "r");
  }

  /** If running frozen soils model **/
#if !LDAS_OUTPUT && !OPTIMIZE
   if(options.FROZEN_SOIL) {
     strcpy(fnames.fdepth, fnames.result_dir);
     strcat(fnames.fdepth, "fdepth");
     strcat(fnames.fdepth, "_");
     strcat(fnames.fdepth, latchar);
     strcat(fnames.fdepth, "_");
     strcat(fnames.fdepth, lngchar);
     if(options.BINARY_OUTPUT) 
       outfp->fdepth = open_file(fnames.fdepth, "wb");
     else outfp->fdepth = open_file(fnames.fdepth, "w");
   }
#endif

  strcpy(fnames.fluxes, fnames.result_dir);
  strcat(fnames.fluxes, "fluxes");
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, latchar);
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, lngchar);
  if(options.BINARY_OUTPUT) 
    outfp->fluxes = open_file(fnames.fluxes, "wb");
  else outfp->fluxes = open_file(fnames.fluxes, "w");

#if !LDAS_OUTPUT && !OPTIMIZE
  strcpy(fnames.snow, fnames.result_dir);
  strcat(fnames.snow, "snow");
  strcat(fnames.snow, "_");
  strcat(fnames.snow, latchar);
  strcat(fnames.snow, "_");
  strcat(fnames.snow, lngchar);
  if(options.BINARY_OUTPUT) 
    outfp->snow = open_file(fnames.snow, "wb");
  else outfp->snow = open_file(fnames.snow, "w");

   if(options.PRT_SNOW_BAND) {
     strcpy(fnames.snowband, fnames.result_dir);
     strcat(fnames.snowband, "snowband");
     strcat(fnames.snowband, "_");
     strcat(fnames.snowband, latchar);
     strcat(fnames.snowband, "_");
     strcat(fnames.snowband, lngchar);
     if(options.BINARY_OUTPUT) 
       outfp->snowband = open_file(fnames.snowband, "wb");
     else outfp->snowband = open_file(fnames.snowband, "w");
   }
#endif

  return (fnames);

} 
