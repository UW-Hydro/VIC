#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id: set_output_defaults.c,v 1.1.2.2 2006/09/26 20:24:27 vicadmin Exp $";

out_data_file_struct *set_output_defaults(out_data_struct *out_data) {
/*************************************************************
  set_output_defaults.c      Ted Bohn     September 08, 2006

  This routine sets the out_data_files and out_data structures to default values.
  These can be overridden by the user in the global control file.

  Modifications:
  2006-Sep-23 Simplified logic for fdepth and snowbands files.  TJB

*************************************************************/

  extern option_struct options;
  out_data_file_struct *out_data_files;
  int v, i;
  int filenum;
  int varnum;

#if OUTPUT_FORCE

  // Output files
  options.Noutfiles = 1;
  out_data_files = (out_data_file_struct *)calloc(options.Noutfiles,sizeof(out_data_file_struct));
  strcpy(out_data_files[0].prefix,"full_data");
  out_data_files[0].nvars = 8;
  out_data_files[0].varid = (int *)calloc(out_data_files[0].nvars, sizeof(int));

  // Variables in first file
  filenum = 0;
  varnum = 0;
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_PREC", varnum++, "%.4f", OUT_TYPE_USINT, 40);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_AIR_TEMP", varnum++, "%.4f", OUT_TYPE_SINT, 100);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SHORTWAVE", varnum++, "%.4f", OUT_TYPE_USINT, 50);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LONGWAVE", varnum++, "%.4f", OUT_TYPE_USINT, 80);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_DENSITY", varnum++, "%.4f", OUT_TYPE_USINT, 100);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_PRESSURE", varnum++, "%.4f", OUT_TYPE_USINT, 100);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_VP", varnum++, "%.4f", OUT_TYPE_SINT, 100);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_WIND", varnum++, "%.4f", OUT_TYPE_USINT, 100);

#else

  // Output files
  options.Noutfiles = 2;
  if (options.FROZEN_SOIL) {
    options.Noutfiles++;
  }
  if (options.PRT_SNOW_BAND) {
    options.Noutfiles++;
  }
  out_data_files = (out_data_file_struct *)calloc(options.Noutfiles,sizeof(out_data_file_struct));
  filenum = 0;
  strcpy(out_data_files[filenum].prefix,"fluxes");
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    out_data_files[filenum].nvars = 20;
  }
  else {
    out_data_files[filenum].nvars = 16;
  }
  filenum++;
  strcpy(out_data_files[filenum].prefix,"snow");
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    out_data_files[filenum].nvars = 7;
  }
  else {
    out_data_files[filenum].nvars = 3;
  }
  if (options.FROZEN_SOIL) {
    filenum++;
    strcpy(out_data_files[filenum].prefix,"fdepth");
    out_data_files[filenum].nvars = 3;
  }
  if (options.PRT_SNOW_BAND) {
    filenum++;
    strcpy(out_data_files[filenum].prefix,"snowband");
    if (options.FULL_ENERGY) {
      out_data_files[filenum].nvars = 13;
    }
    else {
      out_data_files[filenum].nvars = 9;
    }
  }
  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    out_data_files[filenum].varid = (int *)calloc(out_data_files[filenum].nvars, sizeof(int));
  }

  // Variables in first file
  filenum = 0;
  varnum = 0;
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_PREC", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_EVAP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_RUNOFF", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_BASEFLOW", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_WDEW", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SOIL_LIQ", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_RAD_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_NET_SHORT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_R_NET", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LATENT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_EVAP_CANOP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_TRANSP_VEG", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_EVAP_BARE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SUB_CANOP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SUB_SNOW", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SENSIBLE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_GRND_FLUX", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_AERO_RESIST", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SURF_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ALBEDO", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    
  // Variables in second file
  filenum++;
  varnum = 0;
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SWE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_DEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_CANOPY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  if (options.FULL_ENERGY || options.FROZEN_SOIL)  {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ADVECTION", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_DELTACC", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_FLUX", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_REFREEZE_ENERGY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }

  // Variables in other files
  if (options.FROZEN_SOIL) { 
    filenum++;
    varnum = 0;
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_FDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_TDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SOIL_MOIST", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  if (options.PRT_SNOW_BAND) {
    filenum++;
    varnum = 0;
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SWE_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_DEPTH_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_CANOPY_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    if (options.FULL_ENERGY) {
      set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ADVECTION_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
      set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_DELTACC_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
      set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_FLUX_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
      set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_REFREEZE_ENERGY_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_NET_SHORT_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_NET_LONG_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ALBEDO_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LATENT_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SENSIBLE_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_GRND_FLUX_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }

#endif //OUTPUT_FORCE

  return out_data_files;

}
