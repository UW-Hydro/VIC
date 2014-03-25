#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void vicerror(char error_text[])
/**********************************************************************
	vicerror.c	Keith Cherkauer		April 23, 1997

  This subroutine was written to handle numerical errors within the
  VIC model.  This will flush all file buffers so that all records 
  that have been run will be written to disk before the model is exited.

  Modifications:
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures.xi			TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
**********************************************************************/
{
        extern option_struct options;
	extern Error_struct Error;
        filenames_struct fnames;
	void _exit();

        options.COMPRESS=FALSE;	/* turn off compression of last set of files */

	fprintf(stderr,"VIC model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now writing output files...\n");
        close_files(&(Error.filep), Error.out_data_files, &fnames);
	fprintf(stderr,"...now exiting to system...\n");
        fflush(stdout);
        fflush(stderr);
	_exit(1);
}
