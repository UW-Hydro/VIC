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
  11-18-02 Modified to print notification that the output fluxes file
           will be in a binary format.                          LCB
  29-Oct-03 Distinguishing between input lakeparam file and output
	    lake file.						TJB
  2005-Mar-24 Modified to handle ALMA output files.		TJB

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
 
  /********************************
  Input Forcing Files
  ********************************/

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

  /********************************
  Output Files
  ********************************/

#if OUTPUT_FORCE

  strcpy(fnames.fluxes, fnames.result_dir);
  strcat(fnames.fluxes, "full_data");
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, latchar);
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, lngchar);
  if(options.BINARY_OUTPUT) {
    outfp->fluxes = open_file(fnames.fluxes, "wb");
    fprintf(stderr,"flux file opened as binary!\n");
  }
  else outfp->fluxes = open_file(fnames.fluxes, "w");

#else /* OUTPUT_FORCE */

#if LDAS_OUTPUT || OPTIMIZE

  strcpy(fnames.fluxes, fnames.result_dir);
  strcat(fnames.fluxes, "fluxes");
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, latchar);
  strcat(fnames.fluxes, "_");
  strcat(fnames.fluxes, lngchar);
  if(options.BINARY_OUTPUT) 
    outfp->fluxes = open_file(fnames.fluxes, "wb");
  else outfp->fluxes = open_file(fnames.fluxes, "w");

#else /* LDAS_OUTPUT || OPTIMIZE */

  if ( options.ALMA_OUTPUT ) {

    strcpy(fnames.eb, fnames.result_dir);
    strcat(fnames.eb, "eb");
    strcat(fnames.eb, "_");
    strcat(fnames.eb, latchar);
    strcat(fnames.eb, "_");
    strcat(fnames.eb, lngchar);
    outfp->eb = open_file(fnames.eb, "w");

    strcpy(fnames.wb, fnames.result_dir);
    strcat(fnames.wb, "wb");
    strcat(fnames.wb, "_");
    strcat(fnames.wb, latchar);
    strcat(fnames.wb, "_");
    strcat(fnames.wb, lngchar);
    outfp->wb = open_file(fnames.wb, "w");

    strcpy(fnames.sur, fnames.result_dir);
    strcat(fnames.sur, "sur");
    strcat(fnames.sur, "_");
    strcat(fnames.sur, latchar);
    strcat(fnames.sur, "_");
    strcat(fnames.sur, lngchar);
    outfp->sur = open_file(fnames.sur, "w");

    strcpy(fnames.sub, fnames.result_dir);
    strcat(fnames.sub, "sub");
    strcat(fnames.sub, "_");
    strcat(fnames.sub, latchar);
    strcat(fnames.sub, "_");
    strcat(fnames.sub, lngchar);
    outfp->sub = open_file(fnames.sub, "w");

    strcpy(fnames.eva, fnames.result_dir);
    strcat(fnames.eva, "eva");
    strcat(fnames.eva, "_");
    strcat(fnames.eva, latchar);
    strcat(fnames.eva, "_");
    strcat(fnames.eva, lngchar);
    outfp->eva = open_file(fnames.eva, "w");

    strcpy(fnames.csp, fnames.result_dir);
    strcat(fnames.csp, "csp");
    strcat(fnames.csp, "_");
    strcat(fnames.csp, latchar);
    strcat(fnames.csp, "_");
    strcat(fnames.csp, lngchar);
    outfp->csp = open_file(fnames.csp, "w");

  }
  else {

    strcpy(fnames.fluxes, fnames.result_dir);
    strcat(fnames.fluxes, "fluxes");
    strcat(fnames.fluxes, "_");
    strcat(fnames.fluxes, latchar);
    strcat(fnames.fluxes, "_");
    strcat(fnames.fluxes, lngchar);
    if(options.BINARY_OUTPUT) 
      outfp->fluxes = open_file(fnames.fluxes, "wb");
    else outfp->fluxes = open_file(fnames.fluxes, "w");

    strcpy(fnames.snow, fnames.result_dir);
    strcat(fnames.snow, "snow");
    strcat(fnames.snow, "_");
    strcat(fnames.snow, latchar);
    strcat(fnames.snow, "_");
    strcat(fnames.snow, lngchar);
    if(options.BINARY_OUTPUT) 
      outfp->snow = open_file(fnames.snow, "wb");
    else outfp->snow = open_file(fnames.snow, "w");

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

#if LAKE_MODEL
    if ( options.LAKES ) {
      strcpy(fnames.lake, fnames.result_dir);
      strcat(fnames.lake, "lake");
      strcat(fnames.lake, "_");
      strcat(fnames.lake, latchar);
      strcat(fnames.lake, "_");
      strcat(fnames.lake, lngchar);
      if(options.BINARY_OUTPUT) 
        outfp->lake = open_file(fnames.lake, "wb");
      else outfp->lake = open_file(fnames.lake, "w");
    }
#endif // LAKE_MODEL

  }

#endif /* LDAS_OUTPUT || OPTIMIZE */

#endif /* OUTPUT_FORCE - else */

  return (fnames);

} 
