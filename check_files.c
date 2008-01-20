#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void check_files(filep_struct     *filep, 
		 filenames_struct *fnames)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegetation, and global parameters.

  Modifications:
  16-Jun-04 Added logic to count number of fields in first line of
	    soil file.  If there are 54 fields, assume that the
	    final field is avgJulyAirTemp, and set JULY_TAVG_SUPPLIED
	    to TRUE.								TJB
  2006-Sep-01 (Port from 4.1.0) Added logic for OUTPUT_FORCE option.		TJB
  2006-Oct-26 Merged infiles and outfiles structs into filep_struct.		TJB
  2008-Jan-19 Moved check on JULY_TAVG_SUPPLIED into OUTPUT_FORCE section.	TJB

**********************************************************************/
{
  extern option_struct  options;
  extern FILE          *open_file(char string[], char type[]);

  int in_word = 0;
  int num_words = 0;
  int nextchar;

  filep->soilparam   = open_file(fnames->soil, "r");
  if (!options.ARC_SOIL) {
    /* count number of fields in first line of file */
    while ( (nextchar = fgetc(filep->soilparam)) != EOF
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
    fseek(filep->soilparam, 0, 0);
  }
#if !OUTPUT_FORCE
  filep->veglib      = open_file(fnames->veglib, "r");
  filep->vegparam    = open_file(fnames->veg, "r");
  if(options.SNOW_BAND>1)
    filep->snowband    = open_file(fnames->snowband, "r");
#endif /* !OUTPUT_FORCE */
}



