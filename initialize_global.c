#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_global() {
/*********************************************************************
  initialize_global              Keith Cherkauer       March 1998

  This subroutine initalizes all global parameters before they are 
  called by the model.

  option_strudt:           structure containing all global model options
  options.FULL_ENERGY    = TRUE - compute full energy balance
  options.FROZEN_SOIL    = TRUE - compute frozen soils
  options.SNOW_MODEL     = TRUE - compute full energy balance snow model
  options.CALC_SNOW_FLUX = TRUE - compute thermal fluxes through the snowpack
  options.DIST_PRCP      = TRUE - use distributed precipitation
  options.INIT_SOIL      = TRUE - use file to initialize soil column
  options.CORRPREC       = TRUE - correct precipitation measurements to
                           account for gauge loss due to wind
  options.MOISTFRACT     = TRUE - output moisture as fractional moisture
                           content instead of total water in mm
  options.BINARY_OUTPUT  = TRUE - create binary output fles
  options.Nlayer         = Number of soil layers to use in model (at least
                           3 must be used for the energy balance model
  options.GRID_DECIMAL   = Number of decimal places used in the gridded
                           input and output file names
  options.SNOW_BAND      = Number of elevation bands over which to solve the
                           enery balance snow model
 
  debug_struct:            Structure cantains all debugging flags
  debug.DEBUG            = TRUE - turn on all debugging
  debug.PRT_SOIL         = TRUE - print soil parameter debugging files
  debug.PRT_VEGE         = TRUE - print vegetation parameter debugging files
  debug.PRT_GLOBAL       = TRUE - print global parameter debugging files
  debug.PRT_ATMOS        = TRUE - print forcing data debugging files
  debug.PRT_SNOW         = TRUE - print snow pack debugging files
  debug.PRT_FLUX         = TRUE - print energy flux debugging files
  debug.PRT_VAR          = TRUE - print variable debugging files
  debug.PRT_TEMP         = TRUE - 
  debug.PRT_MOIST        = TRUE - 
  debug.PRT_KAPPA        = TRUE - 
  debug.PRT_BALANCE      = TRUE - 
  debug.PRT_GRID         = TRUE - 
  debug.debug_dir        = 

  param_set.SHORTWAVE    = FALSE;
  param_set.LONGWAVE     = FALSE;
  param_set.PRESSURE     = FALSE;
  param_set.TSKC         = FALSE;
  param_set.VP           = FALSE;
  param_set.VPD          = FALSE;
  param_set.REL_HUMID    = FALSE;
  param_set.SPEC_HUMID   = FALSE;
  param_set.ALBEDO       = FALSE;
  param_set.AIR_TEMP     = FALSE;
  param_set.TMAX         = FALSE;
  param_set.TMIN         = FALSE;
  param_set.PREC         = FALSE;
  param_set.WIND         = FALSE;
  param_set.DENSITY      = FALSE;

*********************************************************************/

  extern option_struct options;
  extern debug_struct debug;
  extern param_set_struct param_set;

  options.FULL_ENERGY           = FALSE;
  options.FROZEN_SOIL           = FALSE;
  options.SNOW_MODEL            = FALSE;
  options.CALC_SNOW_FLUX        = TRUE;
  options.DIST_PRCP             = FALSE;
  options.RADAR                 = FALSE;
  options.INIT_SOIL             = FALSE;
  options.CORRPREC              = FALSE;
  options.MOISTFRACT            = FALSE;
  options.BINARY_OUTPUT         = FALSE;
  options.PRT_SNOW_BAND         = FALSE;
  options.Nlayer                = 2;
  options.GRID_DECIMAL          = 2;
  options.SNOW_BAND             = 1;
  options.SNOW_STEP             = 1;
  options.OVERSTORY_ATTENUATION = 0.5;
  options.TRUNK_RATIO           = 0.2;
 
  debug.DEBUG       = FALSE;
  debug.PRT_SOIL    = FALSE;
  debug.PRT_VEGE    = FALSE;
  debug.PRT_GLOBAL  = FALSE;
  debug.PRT_ATMOS   = FALSE;
  debug.PRT_SNOW    = FALSE;
  debug.PRT_FLUX    = FALSE;
  debug.PRT_VAR     = FALSE;
  debug.PRT_TEMP    = FALSE;
  debug.PRT_MOIST   = FALSE;
  debug.PRT_KAPPA   = FALSE;
  debug.PRT_BALANCE = FALSE;
  debug.PRT_GRID    = FALSE;
  strcpy(debug.debug_dir,"./");

  param_set.SHORTWAVE  = FALSE;
  param_set.LONGWAVE   = FALSE;
  param_set.PRESSURE   = FALSE;
  param_set.TSKC       = FALSE;
  param_set.VP         = FALSE;
  param_set.VPD        = FALSE;
  param_set.REL_HUMID  = FALSE;
  param_set.SPEC_HUMID = FALSE;
  param_set.ALBEDO     = FALSE;
  param_set.AIR_TEMP   = FALSE;
  param_set.TMAX       = FALSE;
  param_set.TMIN       = FALSE;
  param_set.PREC       = FALSE;
  param_set.WIND       = FALSE;
  param_set.DENSITY    = FALSE;

}
