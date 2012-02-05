#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void close_files(filep_struct         *filep,
                 out_data_file_struct *out_data_files,
                 filenames_struct     *fnames)
/**********************************************************************
	close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC
  02-27-01 Now closes files opened for lake model applications  KAC
  11-18-02 Now closes lake debugging file.                      LCB
  2003-Oct-29 Distinguishing between input lakeparam file and output
	      lake file.						TJB
  2005-Mar-24 Added support for ALMA output files.			TJB
  2005-Apr-10 Added logic for OUTPUT_FORCE option.			TJB
  2006-Sep-23 Implemented flexible output configuration; uses new
	      out_data_files structure.					TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN

**********************************************************************/
{
  extern option_struct options;
  int filenum;

  /**********************
    Close All Input Files
    **********************/

  fclose(filep->forcing[0]);
  if(options.COMPRESS) compress_files(fnames->forcing[0]);
  if(filep->forcing[1]!=NULL) {
    fclose(filep->forcing[1]);
    if(options.COMPRESS) compress_files(fnames->forcing[1]);
  }

  /*******************
    Close Output Files
    *******************/
  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    fclose(out_data_files[filenum].fh);
    if(options.COMPRESS) compress_files(out_data_files[filenum].filename);
  }

}
