#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

out_data_file_struct *set_output_defaults(out_data_struct *out_data) {
/*************************************************************
  set_output_defaults.c      Ted Bohn     September 08, 2006

  This routine sets the out_data_files and out_data structures to default values.
  These can be overridden by the user in the global control file.

  Modifications:
  2006-Oct-10 Shortened the names of output variables whose names were
              too long.							TJB
  2007-Oct-09 Updated to reflect variables present in traditional 4.1.0
	      output files.  Previously the defaults matched the traditional
	      4.0.6 output files.					TJB
  2008-Apr-11 Added OUT_SUB_BLOWING, OUT_SUB_SURFACE, and OUT_SUB_SNOW to
	      default snow output file for case of options.BLOWING == TRUE.
	      This makes it almost the same as previous versions of 4.1.0,
	      (r3 and earlier) with the exception that previous versions
	      of 4.1.0 multiplied these terms by 100 when saving to the
	      snow file.						TJB
  2010-Sep-24 Renamed RUNOFF_IN and OUT_RUNOFF_IN to CHANNEL_IN and
	      OUT_LAKE_CHAN_IN, respectively.  Renamed OUT_EVAP_LAKE
	      to OUT_LAKE_EVAP.  Added other lake water balance terms
	      to set of output variables.  Added volumetric versions
	      of these too.						TJB
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct.			TJB
*************************************************************/

  extern option_struct options;
  out_data_file_struct *out_data_files;
  int v, i;
  int filenum;
  int varnum;

if(options.OUTPUT_FORCE) {

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

}
else {

  // Output files
  options.Noutfiles = 2;
  if (options.FROZEN_SOIL) {
    options.Noutfiles++;
  }
  if (options.PRT_SNOW_BAND) {
    options.Noutfiles++;
  }
  if (options.LAKES) {
    options.Noutfiles++;
  }
  out_data_files = (out_data_file_struct *)calloc(options.Noutfiles,sizeof(out_data_file_struct));
  filenum = 0;
  strcpy(out_data_files[filenum].prefix,"fluxes");
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    out_data_files[filenum].nvars = 26;
  }
  else {
    out_data_files[filenum].nvars = 20;
  }
  filenum++;
  strcpy(out_data_files[filenum].prefix,"snow");
  if (options.FULL_ENERGY || options.FROZEN_SOIL) {
    out_data_files[filenum].nvars = 14;
  }
  else {
    out_data_files[filenum].nvars = 4;
  }
  if (options.BLOWING) {
    out_data_files[filenum].nvars+= 3;
  }
  if (options.FROZEN_SOIL) {
    filenum++;
    strcpy(out_data_files[filenum].prefix,"fdepth");
    out_data_files[filenum].nvars = 4;
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
  if (options.LAKES) {
    filenum++;
    strcpy(out_data_files[filenum].prefix,"lake");
    out_data_files[filenum].nvars = 8;
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
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_DELTAH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_FUSION", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_AERO_RESIST", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SURF_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ALBEDO", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_REL_HUMID", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_IN_LONG", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_AIR_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_WIND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    
  // Variables in second file
  filenum++;
  varnum = 0;
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SWE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_DEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_CANOPY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_COVER", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  if (options.FULL_ENERGY || options.FROZEN_SOIL)  {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ADVECTION", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_DELTACC", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_FLUX", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_RFRZ_ENERGY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_MELT_ENERGY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ADV_SENS", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LATENT_SUB", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_SURF_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_PACK_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SNOW_MELT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  if (options.BLOWING)  {
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SUB_BLOWING", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SUB_SURFACE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SUB_SNOW", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }

  // Variables in other files
  if (options.FROZEN_SOIL) { 
    filenum++;
    varnum = 0;
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_FDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_TDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SOIL_MOIST", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SURF_FROST_FRAC", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
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
      set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_RFRZ_ENERGY_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_NET_SHORT_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_NET_LONG_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_ALBEDO_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LATENT_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_SENSIBLE_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_GRND_FLUX_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }
  if (options.LAKES) {
    filenum++;
    varnum = 0;
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_ICE_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_ICE_HEIGHT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_ICE_FRACT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_DEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_SURF_AREA", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_VOLUME", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_SURF_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, TRUE, filenum, out_data, "OUT_LAKE_EVAP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
  }

} // !OUTPUT_FORCE

  return out_data_files;

}
