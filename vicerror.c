#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

void vicerror(char error_text[])
/**********************************************************************
	vicerror.c	Keith Cherkauer		April 23, 1997

  This subroutine was written to handle numerical errors within the
  VIC model.  This will dump all records currently stored, as well
  current values of all data structures.

**********************************************************************/
{
        extern option_struct options;
	extern Error_struct Error;
        extern debug_struct debug;

        filenames_struct fnames;
	void _exit();
        int i;

        options.COMPRESS=FALSE;	/* turn off compression of last set of files */

	fprintf(stderr,"VIC model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now writing output files...\n");
/*****
        write_data(Error.out_data,Error.outfp,Error.rec,Error.dt);
*****/
        close_files(Error.infp,Error.outfp,fnames);
/*****
        for(i=0;i<Error.veg_con[0].vegetat_type_num;i++) {
	  write_vegvar(Error.veg_var[i],i);
        }
*****/
	fprintf(stderr,"...now exiting to system...\n");
        fflush(stdout);
        fflush(stderr);
	_exit(1);
}
