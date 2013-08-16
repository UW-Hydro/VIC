#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: open_debug.c,v 4.1.2.1 2004/05/10 18:42:55 tbohn Exp $";

#if LINK_DEBUG
void open_debug() {
/**********************************************************************
  open_debug		Keith Cherkauer		October 10, 1997

  This subroutine opens all requested debugging output files.

  Modifications:
    07-May-04 Initialize debug_store_moist array when debug.PRT_MOIST
	      is true (as well as under the other previously-defined
	      conditions).					TJB

**********************************************************************/

  extern debug_struct debug;
  extern option_struct options;

  char tempname[MAXSTRING];
  int  i, j;
  
  if(debug.DEBUG || debug.PRT_TEMP) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_temp.out");
    if((debug.fg_temp=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_temp.out");
  }
  if(debug.DEBUG || debug.PRT_MOIST) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_moist.out");
    if((debug.fg_moist=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_moist.out");
  }
  if(debug.DEBUG || debug.PRT_KAPPA) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_kappa.out");
    if((debug.fg_kappa=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_kappa.out");
  }
  if(debug.DEBUG || debug.PRT_BALANCE || debug.PRT_MOIST) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_balance.out");
    if((debug.fg_balance=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_balance.out");
    for(i=0;i<2;i++) {
      debug.inflow[i]      = (double **)calloc(options.SNOW_BAND,
					       sizeof(double*));
      debug.outflow[i]     = (double **)calloc(options.SNOW_BAND,
					       sizeof(double*));
      debug.store_moist[i] = (double **)calloc(options.SNOW_BAND,
					       sizeof(double*));
      for(j=0;j<options.SNOW_BAND;j++) {
	debug.inflow[i][j]      = (double *)calloc(options.Nlayer+3,
						   sizeof(double));
	debug.outflow[i][j]     = (double *)calloc(options.Nlayer+3,
						   sizeof(double));
	debug.store_moist[i][j] = (double *)calloc(options.Nlayer+3,
						   sizeof(double));
      }
    }
  }
  if(debug.DEBUG || debug.PRT_FLUX) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_energy.out");
    if((debug.fg_energy=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_energy.out");
  }
  if(debug.DEBUG || debug.PRT_SNOW) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_snow.out");
    if((debug.fg_snow=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_snow.out");
  }
  if(debug.DEBUG || debug.PRT_GRID) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_grid.out");
    if((debug.fg_grid=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_grid.out");
  }
  if(debug.DEBUG || debug.PRT_ATMOS) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_snowstep_atmos.out");
    if((debug.fg_snowstep_atmos=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_snowstep_atmos.out");
  }
  if(debug.DEBUG || debug.PRT_ATMOS) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_modelstep_atmos.out");
    if((debug.fg_modelstep_atmos=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_modelstep_atmos.out");
  }
}
#endif
