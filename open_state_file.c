#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id $";

FILE *open_state_file(global_param_struct *global,
		      int                  Nlayer,
		      int                  Nnodes) 
/*********************************************************************
  open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file for output.

  Modifications:
  04-10-03 Modified to open and write to a binary state file.    KAC
  06-03-03 modified to handle both ASCII and BINARY state files.  KAC
  10-02-03 Modified to write statefile as year, month, day rather than
           day, month, year.  This makes it consistent with how the
           file is read by the model.                            KAC
  11-May-04 Modified the statefile name to contain year, month, day
	    rather than day, month, year.  This makes it consistent
	    with the planned release of 4.1.0.			TJB
  2005-11-09 Removed '#if SAVE_STATE'                           GCT
  2005-11-10 Moved setting of statename to get_global_param     GCT
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... GCT
*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
  sprintf(filename,"%s", global->statename);
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

