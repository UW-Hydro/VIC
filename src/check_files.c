#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void check_files(filep_struct     *filep, 
		 filenames_struct *fnames)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegetation, and global parameters.

  Modifcations:
  02-27-01 Added controls for lake model parameter file    KAC
  2005-Apr-13 Added logic for OUTPUT_FORCE option.			TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct.			TJB
**********************************************************************/
{
  extern option_struct  options;
  extern FILE          *open_file(char string[], char type[]);

  filep->soilparam   = open_file(fnames->soil, "r");
  if (!options.OUTPUT_FORCE) {
    filep->veglib      = open_file(fnames->veglib, "r");
    filep->vegparam    = open_file(fnames->veg, "r");
    if(options.SNOW_BAND>1)
      filep->snowband    = open_file(fnames->snowband, "r");
    if ( options.LAKES )
      filep->lakeparam = open_file(fnames->lakeparam,"r");
  }

}


