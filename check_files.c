#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void check_files(infiles_struct *infp, 
		 filenames_struct fnames)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegitation, and global parameters.

**********************************************************************/
{
  extern FILE *open_file(char string[], char type[]);

  infp->soilparam   = open_file(fnames.soil, "r");
  infp->veglib      = open_file(fnames.veglib, "r");
  infp->vegparam    = open_file(fnames.veg, "r");
}



