#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
 
void close_files(infiles_struct   inf,
                 outfiles_struct  outf,
                 filenames_struct fnames)
/**********************************************************************
	close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  /**********************
    Close All Input Files
    **********************/

  fclose(inf.forcing[0]);
  if(options.COMPRESS) compress_files(fnames.forcing[0]);
  if(strcasecmp(fnames.forcing[1],"FALSE")!=0) {
    fclose(inf.forcing[1]);
    if(options.COMPRESS) compress_files(fnames.forcing[1]);
  }

  /*******************
    Close Output Files
    *******************/

  /** Energy and Moisture Fluxes Output File **/
  fclose(outf.fluxes);
  if(options.COMPRESS) compress_files(fnames.fluxes);

  /** Frozen Soils Output File **/
  if(options.FROZEN_SOIL) {
    fclose(outf.fdepth);
    if(options.COMPRESS) compress_files(fnames.fdepth);
  }

  /** Snow Data Output File **/
  if(options.FULL_ENERGY || options.SNOW_MODEL) {
    fclose(outf.snow);
    if(options.COMPRESS) compress_files(fnames.snow);
  }

  /*******************************
    Close All Used Debugging Files
    *******************************/ 

  if(debug.DEBUG || debug.PRT_TEMP) {
    fclose(debug.fg_temp);
  }
  if(debug.DEBUG || debug.PRT_MOIST) {
    fclose(debug.fg_moist);
  }
  if(debug.DEBUG || debug.PRT_KAPPA) {
    fclose(debug.fg_kappa);
  }
  if(debug.DEBUG || debug.PRT_BALANCE) {
    fclose(debug.fg_balance);
  }
  if(debug.DEBUG || debug.PRT_FLUX) {
    fclose(debug.fg_energy);
  }
  if(debug.DEBUG || debug.PRT_SNOW) {
    fclose(debug.fg_snow);
  }
  if(debug.DEBUG || debug.PRT_GRID) {
    fclose(debug.fg_grid);
  }

}
