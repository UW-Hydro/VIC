#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id$";

double **read_forcing_data(FILE                **infile,
			   global_param_struct   global_param)
/**********************************************************************
  read_forcing_data    Keith Cherkauer      January 10, 2000

  This subroutine controls the order and number of forcing variables
  read from the forcing data files.  Two forcing files are allowed, 
  variables, time step and file format must be defined in the global
  control file.

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NR, NF;

  char                 errorstr[MAXSTRING];
  int                  i;
  double             **forcing_data;

  /** Allocate data arrays for input forcing data **/
  forcing_data = (double **)calloc(N_FORCING_TYPES,sizeof(double*));
  for(i=0;i<N_FORCING_TYPES;i++) 
    if (param_set.TYPE[i].SUPPLIED) 
      forcing_data[i] = (double *)calloc((global_param.nrecs * NF),
			   sizeof(double));

  /** Read First Forcing Data File **/
  if(param_set.FORCE_DT[0] > 0) {
    read_atmos_data(infile[0], global_param, 0, global_param.forceskip[0],
		    forcing_data);
  }
  else {
    sprintf(errorstr,"ERROR: File time step must be defined for at least the first forcing file (FILE_DT).\n");
    vicerror(errorstr);
  }

  /** Read Second Forcing Data File **/
  if(param_set.FORCE_DT[1] > 0) {
    read_atmos_data(infile[1], global_param, 1, global_param.forceskip[1], 
		    forcing_data);
  }

  return(forcing_data);

}
