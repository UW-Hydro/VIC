#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id $";

#if SAVE_STATE

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

*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
  sprintf(filename,"%s_%04i%02i%02i", global->statename, 
	  global->stateyear, global->statemonth, global->stateday);
  if ( options.BINARY_STATE_FILE )
    statefile = open_file(filename,"wb");
  else
    statefile = open_file(filename,"w");

  /* Write save state date information */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &global->stateyear, 1, sizeof(int), statefile );
    fwrite( &global->statemonth, 1, sizeof(int), statefile );
    fwrite( &global->stateday, 1, sizeof(int), statefile );
  }
  else {
    fprintf(statefile,"%i %i %i\n", global->stateyear, 
	    global->statemonth, global->stateday);
  }

  /* Write simulation flags */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &Nlayer, 1, sizeof(int), statefile );
    fwrite( &Nnodes, 1, sizeof(int), statefile );
  }
  else {
    fprintf(statefile,"%i %i\n", Nlayer, Nnodes);
  }

  return(statefile);

}

#endif
