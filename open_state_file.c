#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#if SAVE_STATE

FILE *open_state_file(global_param_struct *global,
		      int                  Nlayer,
		      int                  Nnodes) 
/*********************************************************************
  open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file for output.

  Modifications:
  04-10-03 Modified to open and write to a binary state file.    KAC

*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  double  Nsum;

  /* open state file */
  sprintf(filename,"%s_%02i%02i%04i", global->statename, 
	  global->stateday, global->statemonth, global->stateyear);
  statefile = open_file(filename,"wb");

  /* Write save state date information */
  /*fprintf(statefile,"%i %i %i\n", global->stateday, 
	  global->statemonth, global->stateyear);*/
  fwrite( &global->stateyear, 1, sizeof(int), statefile );
  fwrite( &global->statemonth, 1, sizeof(int), statefile );
  fwrite( &global->stateday, 1, sizeof(int), statefile );

  /* Write simulation flags */
  //fprintf(statefile,"%i %i\n", Nlayer, Nnodes);
  fwrite( &Nlayer, 1, sizeof(int), statefile );
  fwrite( &Nnodes, 1, sizeof(int), statefile );

  return(statefile);

}

#endif
