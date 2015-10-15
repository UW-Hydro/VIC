#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id$";

double **read_forcing_data(FILE                **infile,
			   global_param_struct   global_param,
			   double            ****veg_hist_data)
/**********************************************************************
  read_forcing_data    Keith Cherkauer      January 10, 2000

  This subroutine controls the order and number of forcing variables
  read from the forcing data files.  Two forcing files are allowed, 
  variables, time step and file format must be defined in the global
  control file.

  Modifications:
  2014-Apr-25 Added non-climatological veg parameters (as forcing
	      variables).						TJB
  2014-Apr-25 Added partial vegcover fraction.				TJB
**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern int              NR, NF;

  char                 errorstr[MAXSTRING];
  int                  i,j;
  double             **forcing_data;

  /** Allocate data arrays for input forcing data **/
  forcing_data = (double **)calloc(N_FORCING_TYPES,sizeof(double*));
  (*veg_hist_data) = (double ***)calloc(N_FORCING_TYPES,sizeof(double**));
  for(i=0;i<N_FORCING_TYPES;i++) {
    if (param_set.TYPE[i].SUPPLIED) {
      if (i != ALBEDO && i != LAI_IN && i != VEGCOVER) {
        forcing_data[i] = (double *)calloc((global_param.nrecs * NF), sizeof(double));
      }
      else {
        (*veg_hist_data)[i] = (double **)calloc(param_set.TYPE[i].N_ELEM, sizeof(double*));
        for(j=0;j<param_set.TYPE[i].N_ELEM;j++) {
          (*veg_hist_data)[i][j] = (double *)calloc((global_param.nrecs * NF), sizeof(double));
        }
      }
    }
  }

  /** Read First Forcing Data File **/
  if(param_set.FORCE_DT[0] > 0) {
    read_atmos_data(infile[0], global_param, 0, global_param.forceskip[0],
		    forcing_data, (*veg_hist_data));
  }
  else {
    sprintf(errorstr,"ERROR: File time step must be defined for at least the first forcing file (FILE_DT).\n");
    vicerror(errorstr);
  }

  /** Read Second Forcing Data File **/
  if(param_set.FORCE_DT[1] > 0) {
    read_atmos_data(infile[1], global_param, 1, global_param.forceskip[1], 
		    forcing_data, (*veg_hist_data));
  }
  if(param_set.FORCE_DT[2] > 0) {
    read_atmos_data(infile[2], global_param, 2, global_param.forceskip[2], 
		    forcing_data, (*veg_hist_data));
  }

  return(forcing_data);

}
