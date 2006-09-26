#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void close_files(infiles_struct       *inf,
                 out_data_file_struct *out_data_files,
                 filenames_struct     *fnames)
/**********************************************************************
	close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC
  2006-Sep-01 (Port from 4.1.0) Added logic for OUTPUT_FORCE option. TJB
  2006-Sep-11 Implemented flexible output configuration, using
              new out_data_files structure. TJB

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int filenum;

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

  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    fclose(out_data_files[filenum].fh);
    if(options.COMPRESS) compress_files(out_data_files[filenum].filename);
  }

#if !OUTPUT_FORCE

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
#endif /* LINK_DEBUG */
#endif /* !OUTPUT_FORCE */

}
