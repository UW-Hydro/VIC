#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

FILE *check_state_file(char                *init_state,
		       dmy_struct           dmy,
		       global_param_struct *global,
		       int                  Nlayer,
		       int                  Nnodes) 
/*********************************************************************
  check_state_file      Keith Cherkauer           April 17, 2000

  This subroutine opens a model state file and verifys that the 
  starting date, number of layers and number of thermal nodes in the 
  file agrees with what was defined in the model global control file.

*********************************************************************/
{
  extern option_struct options;

  FILE   *statefile;
  char    filename[MAXSTRING];
  char    ErrStr[MAXSTRING];
  double  Nsum;
  int     tmp_Nlayer;
  int     tmp_Nnodes;
  int     startday, startmonth, startyear;

  /* open state file */
  statefile = open_file(init_state,"r");

  /* Check state date information */
  fscanf(statefile,"%i %i %i\n", &startday, &startmonth, &startyear);
  if ( startday != dmy.day || startmonth != dmy.month || 
       startyear != dmy.year || dmy.hour != 0 ) {
    sprintf(ErrStr,"Starting date of the state file (%i/%i/%i) does not match that of the global file (%i/%i/%i).  Model cannot be restarted, check your files.  NOTE: Model has to be restarted at hour 0 (currently %i) of the day indicated at the top of the state file.", startmonth, startday, startyear, dmy.month, dmy.day, dmy.year, dmy.hour);
    nrerror(ErrStr);
  }

  /* Check simulation options */
  fscanf(statefile,"%i %i\n", &tmp_Nlayer, &tmp_Nnodes);
  if ( tmp_Nlayer != Nlayer ) {
    sprintf(ErrStr,"The number of soil moisture layers in the model state file (%i) does not equal that defined in the global control file (%i).  Check your input files.", tmp_Nlayer, Nlayer);
    nrerror(ErrStr);
  }
  if ( tmp_Nnodes != Nnodes ) {
    sprintf(ErrStr,"The number of soil thermal nodes in the model state file (%i) does not equal that defined in the global control file (%i).  Check your input files.", tmp_Nnodes, Nnodes);
    nrerror(ErrStr);
  }

  return(statefile);

}
