#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

void initialize_global() {

  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;

  options.FULL_ENERGY = FALSE;
  options.FROZEN_SOIL = FALSE;
  options.SNOW_MODEL = FALSE;
  options.CALC_SNOW_FLUX = FALSE;
  options.DIST_PRCP = FALSE;
  options.RADAR = FALSE;
  options.INIT_SOIL = FALSE;
  options.Nlayer = 2;
  options.CORRPREC = FALSE;
  options.MOISTFRACT = FALSE;
  options.GRID_DECIMAL = 2;
 
  debug.DEBUG = FALSE;
  debug.PRT_SOIL = FALSE;
  debug.PRT_VEGE = FALSE;
  debug.PRT_GLOBAL = FALSE;
  debug.PRT_ATMOS = FALSE;
  debug.PRT_SNOW = FALSE;
  debug.PRT_FLUX = FALSE;
  debug.PRT_VAR = FALSE;
  debug.PRT_TEMP = FALSE;
  debug.PRT_MOIST = FALSE;
  debug.PRT_KAPPA = FALSE;
  debug.PRT_BALANCE = FALSE;
  debug.PRT_GRID = FALSE;
  strcpy(debug.debug_dir,"./");

  param_set.SHORTWAVE = FALSE;
  param_set.LONGWAVE = FALSE;
  param_set.PRESSURE = FALSE;
  param_set.TSKC = FALSE;
  param_set.SVP = FALSE;
  param_set.VP = FALSE;
  param_set.VPD = FALSE;
  param_set.REL_HUMID = FALSE;
  param_set.SPEC_HUMID = FALSE;
  param_set.ALBEDO = FALSE;
  param_set.AIR_TEMP = FALSE;
  param_set.TMAX = FALSE;
  param_set.TMIN = FALSE;
  param_set.PREC = FALSE;
  param_set.WIND = FALSE;
  param_set.DENSITY = FALSE;

}
