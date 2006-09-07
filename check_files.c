#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void check_files(infiles_struct   *infp, 
		 filenames_struct *fnames)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegetation, and global parameters.

  Modifications:
  16-Jun-04 Added logic to count number of fields in first line of
	    soil file.  If there are 54 fields, assume that the
	    final field is avgJulyAirTemp, and set JULY_TAVG_SUPPLIED
	    to TRUE.						TJB
  2006-Sep-01 (Port from 4.1.0) Added logic for OUTPUT_FORCE option. TJB
**********************************************************************/
{
  extern option_struct  options;
  extern FILE          *open_file(char string[], char type[]);

  int in_word = 0;
  int num_words = 0;
  int nextchar;

  infp->soilparam   = open_file(fnames->soil, "r");
#if !OUTPUT_FORCE
  if (!options.ARC_SOIL) {
    /* count number of fields in first line of file */
    while ( (nextchar = fgetc(infp->soilparam)) != EOF
      && (char)nextchar != '\n' ) {
      if ((char)nextchar != ' ' && (char)nextchar != '\t') {
        if (!in_word) {
          /* we've reached the beginning of a word */
          in_word = 1;
          num_words++;
        }
      }
      else {
        if (in_word) {
          /* we've reached the beginning of whitespace */
          in_word = 0;
        }
      }
    }
    if (num_words == 54) {
      options.JULY_TAVG_SUPPLIED = TRUE;
    }
    /* rewind file back to beginning */
    fseek(infp->soilparam, 0, 0);
  }
  infp->veglib      = open_file(fnames->veglib, "r");
  infp->vegparam    = open_file(fnames->veg, "r");
  if(options.SNOW_BAND>1)
    infp->snowband    = open_file(fnames->snow_band, "r");
#endif /* !OUTPUT_FORCE */
}



