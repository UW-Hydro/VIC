#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

filenames_struct make_in_and_outfiles(infiles_struct *infp, 
				      filenames_struct filenames,
				      soil_con_struct soil,
				      outfiles_struct *outfp)
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
  extern option_struct options;
  extern debug_struct debug;

  char latchar[10], lngchar[10], junk[5];
  int layer;
  filenames_struct fnames;
  extern FILE *open_file(char string[], char type[]);

  sprintf(junk, "%%.%if", options.GRID_DECIMAL);
  sprintf(latchar, junk, soil.lat);
  sprintf(lngchar, junk, soil.lng);
 
  strcat(filenames.forcing[0], latchar);
  strcat(filenames.forcing[0], "_");
  strcat(filenames.forcing[0], lngchar);
  infp->forcing[0] = open_file(filenames.forcing[0], "r");

  infp->forcing[1] = NULL;
  if(strcasecmp(filenames.forcing[1],"FALSE")!=0) {
    strcat(filenames.forcing[1], latchar);
    strcat(filenames.forcing[1], "_");
    strcat(filenames.forcing[1], lngchar);
    infp->forcing[1] = open_file(filenames.forcing[1], "r");
    if((strcasecmp(options.FORCE_TYPE,"SAWD"))==0) options.HP = FALSE;
  }
  else if((strcasecmp(options.FORCE_TYPE,"SAWD"))==0) options.HP = TRUE;

  /** If using radar precipitation data **/
  if(options.RADAR) {
    strcpy(filenames.radar_prcp, "../radar_prec/");
    strcat(filenames.radar_prcp, "radar_prec");
    strcat(filenames.radar_prcp, "_");
    strcat(filenames.radar_prcp, latchar);
    strcat(filenames.radar_prcp, "_");
    strcat(filenames.radar_prcp, lngchar);
    infp->radar_prcp = open_file(filenames.radar_prcp, "r");
  }

  /** Initial Soil Variables Read from File **/
  if((options.FULL_ENERGY || options.FROZEN_SOIL)&& options.INIT_SOIL) {
    infp->init_soil = open_file(filenames.init_soil, "r");
  }

  /** Initial Snow Variables Read from File **/
  if((options.FULL_ENERGY || options.FROZEN_SOIL)&& options.INIT_SNOW) {
    infp->init_snow = open_file(filenames.init_snow, "r");
  }

  /** If running frozen soils model **/
  if(options.FROZEN_SOIL) {
    strcpy(filenames.fdepth, filenames.result_dir);
    strcat(filenames.fdepth, "fdepth");
    strcat(filenames.fdepth, "_");
    strcat(filenames.fdepth, latchar);
    strcat(filenames.fdepth, "_");
    strcat(filenames.fdepth, lngchar);
    outfp->fdepth = open_file(filenames.fdepth, "w");
  }

  strcpy(filenames.fluxes, filenames.result_dir);
  strcat(filenames.fluxes, "fluxes");
  strcat(filenames.fluxes, "_");
  strcat(filenames.fluxes, latchar);
  strcat(filenames.fluxes, "_");
  strcat(filenames.fluxes, lngchar);
  outfp->fluxes = open_file(filenames.fluxes, "w");

  if(options.FULL_ENERGY || options.SNOW_MODEL) {
    strcpy(filenames.snow, filenames.result_dir);
    strcat(filenames.snow, "snow");
    strcat(filenames.snow, "_");
    strcat(filenames.snow, latchar);
    strcat(filenames.snow, "_");
    strcat(filenames.snow, lngchar);
    outfp->snow = open_file(filenames.snow, "w");
  }

  fnames = filenames;
  return (fnames);

} 
