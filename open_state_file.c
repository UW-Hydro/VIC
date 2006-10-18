#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";


FILE *open_state_file(global_param_struct *global,
		      filenames_struct     filenames,
		      int                  Nlayer,
		      int                  Nnodes) 
/*********************************************************************
  open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file for output.

  Modifications:
  04-10-03 Modified to open and write to a binary state file.    KAC
  06-03-03 modified to handle both ASCII and BINARY state files.  KAC
  2005-11-29 SAVE_STATE is set in global param file, not in user_def.h GCT
  2005-12-06 Moved setting of statename to get_global_param     GCT
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to filenames->statefile. TJB

*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
  sprintf(filename,"%s", filenames.statefile);
  if ( options.BINARY_STATE_FILE )
    statefile = open_file(filename,"wb");
  else
    statefile = open_file(filename,"w");

  /* Write save state date information */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &global->stateyear, sizeof(int), 1, statefile );
    fwrite( &global->statemonth, sizeof(int), 1, statefile );
    fwrite( &global->stateday, sizeof(int), 1, statefile );
  }
  else {
    fprintf(statefile,"%i %i %i\n", global->stateyear, 
	    global->statemonth, global->stateday);
  }

  /* Write simulation flags */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &Nlayer, sizeof(int), 1, statefile );
    fwrite( &Nnodes, sizeof(int), 1, statefile );
  }
  else {
    fprintf(statefile,"%i %i\n", Nlayer, Nnodes);
  }

  return(statefile);

}

