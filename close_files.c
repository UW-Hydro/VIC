#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void close_files(infiles_struct   *inf,
                 outfiles_struct  *outf,
                 filenames_struct *fnames)
/**********************************************************************
	close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC
  02-27-01 Now closes files opened for lake model applications  KAC
  11-18-02 Now closes lake debugging file.                      LCB
  29-Oct-03 Distinguishing between input lakeparam file and output
	    lake file.						TJB
  2005-Mar-24 Added support for ALMA output files.		TJB

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  /**********************
    Close All Input Files
    **********************/

  fclose(inf->forcing[0]);
  if(options.COMPRESS) compress_files(fnames->forcing[0]);
  if(inf->forcing[1]!=NULL) {
    fclose(inf->forcing[1]);
    if(options.COMPRESS) compress_files(fnames->forcing[1]);
  }

  /*******************
    Close Output Files
    *******************/

#if LDAS_OUTPUT || OPTIMIZE

  /** Energy and Moisture Fluxes Output File **/
  fclose(outf->fluxes);
  if(options.COMPRESS) compress_files(fnames->fluxes);

#else /* LDAS_OUTPUT || OPTIMIZE */

  if ( options.ALMA_OUTPUT ) {

    /** Energy Balance Output File **/
    fclose(outf->eb);
    if(options.COMPRESS) compress_files(fnames->eb);

    /** Water Balance Output File **/
    fclose(outf->wb);
    if(options.COMPRESS) compress_files(fnames->wb);

    /** Surface Output File **/
    fclose(outf->sur);
    if(options.COMPRESS) compress_files(fnames->sur);

    /** Subsurace Output File **/
    fclose(outf->sub);
    if(options.COMPRESS) compress_files(fnames->sub);

    /** Evaporation Output File **/
    fclose(outf->eva);
    if(options.COMPRESS) compress_files(fnames->eva);

    /** Cold-season Processes Output File **/
    fclose(outf->csp);
    if(options.COMPRESS) compress_files(fnames->csp);

  }
  else {

    /** Energy and Moisture Fluxes Output File **/
    fclose(outf->fluxes);
    if(options.COMPRESS) compress_files(fnames->fluxes);

    /** These output files are not used when using LDAS binary format **/

    /** Frozen Soils Output File **/
    if(options.FROZEN_SOIL) {
      fclose(outf->fdepth);
      if(options.COMPRESS) compress_files(fnames->fdepth);
    }

    /** Snow Data Output File **/
    fclose(outf->snow);
    if(options.COMPRESS) compress_files(fnames->snow);

    if(options.PRT_SNOW_BAND) {
      fclose(outf->snowband);
      if(options.COMPRESS) compress_files(fnames->snowband);
    }

#if LAKE_MODEL
    if ( options.LAKES ) {
      /** Lake Data Output File **/
      fclose(outf->lake);
      if(options.COMPRESS) compress_files(fnames->lake);
    }
#endif /* LAKE_MODEL */

  }

#endif /* LDAS_OUTPUT || OPTIMIZE */

  /*******************************
    Close All Used Debugging Files
    *******************************/ 

#if LINK_DEBUG
  if(debug.DEBUG || debug.PRT_TEMP) {
    fclose(debug.fg_temp);
  }
  if(debug.DEBUG || debug.PRT_MOIST) {
    fclose(debug.fg_moist);
  }
  if(debug.DEBUG || debug.PRT_KAPPA) {
    fclose(debug.fg_kappa);
  }
  if(debug.DEBUG || debug.PRT_LAKE) {
    fclose(debug.fg_lake);
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
#endif

}
