#include <stdio.h>
#include <string.h>
#include <vicNl.h>

void open_debug() {
/**********************************************************************
  open_debug		Keith Cherkauer		October 10, 1997

  This subroutine opens all requested debugging output files.

**********************************************************************/

  extern debug_struct debug;
  extern option_struct options;

  char tempname[MAXSTRING];

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
  if(debug.DEBUG || debug.PRT_BALANCE) {
    strcpy(tempname,debug.debug_dir);
    strcat(tempname,"/VIC_balance.out");
    if((debug.fg_balance=fopen(tempname,"w"))==NULL)
      nrerror("ERROR: Unable to open VIC_balance.out");
    debug.inflow = (double *)calloc(options.Nlayer+3,sizeof(double));
    debug.outflow = (double *)calloc(options.Nlayer+3,sizeof(double));
    debug.store_moist = (double *)calloc(options.Nlayer+3,sizeof(double));
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
      nrerror("ERROR: Unable to open VIC_grid.xyz");
  }

}
