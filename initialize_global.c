#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_global() {
/*********************************************************************
  initialize_global              Keith Cherkauer       March 1998

  This subroutine initalizes all global parameters before they are 
  called by the model.

  option_struct:           structure containing all global model options
  options.FULL_ENERGY    = TRUE - compute full energy balance
  options.FROZEN_SOIL    = TRUE - compute frozen soils
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

  param_set.ALBEDO       = FALSE;
  param_set.AIR_TEMP     = FALSE;
  param_set.CRAINF       = FALSE;
  param_set.CSNOWF       = FALSE;
  param_set.DENSITY      = FALSE;
  param_set.LONGWAVE     = FALSE;
  param_set.LSRAINF      = FALSE;
  param_set.LSSNOWF      = FALSE;
  param_set.PREC         = FALSE;
  param_set.PRESSURE     = FALSE;
  param_set.QAIR         = FALSE;
  param_set.RAINF        = FALSE;
  param_set.REL_HUMID    = FALSE;
  param_set.SHORTWAVE    = FALSE;
  param_set.SNOWF        = FALSE;
  param_set.TMAX         = FALSE;
  param_set.TMIN         = FALSE;
  param_set.TSKC         = FALSE;
  param_set.VP           = FALSE;
  param_set.WIND         = FALSE;
  param_set.WIND_E       = FALSE;
  param_set.WIND_N       = FALSE;

  Modifications:
    06-03-2003 modified to handle both ASCII and BINARY state files.  KAC
    09-02-2003 Moved COMPUTE_TREELINE flag from user_def.h to the 
               options structure.  Now when not set to FALSE, the 
               value indicates the default above treeline vegetation
               if no usable vegetation types are in the grid cell 
               (i.e. everything has a canopy).  A negative value  
               will cause the model to use bare soil.  Make sure that 
               positive index value refer to a non-canopied vegetation
               type in the vegetation library.                   KAC
    07-May-04 Initialize ARC_SOIL, COMPRESS, and ARNO_PARAMS to FALSE.
	      Also changed limit on loop over forcing types from
	      hard-coded 17 to variable N_FORCING_TYPES.	TJB
    16-Jun-04 Added JULY_TAVG_SUPPLIED.				TJB
   2005-11-21 (Port from 4.1.0) Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW. GCT 	      
  2006-Jan-22    Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option. TJB
  2006-Sep-11 Implemented flexible output configuration; added new
              options.Noutfiles and organized the options according
              to function. TJB
  2006-Sep-14 (Port from 4.1.0) Added support for ALMA-compliant input and output. TJB
  2006-Dec-29 Added REL_HUMID to the list of supported met input variables. TJB
  2006-Dec-29 Added initialization of some uninitialized met input variables. TJB
  2007-Jan-02 Added CSNOWF and LSSNOWF to list of supported met input variables. TJB
  2007-Jan-02 Added ALMA_INPUT option; removed TAIR and PSURF from list of supported met input variables. TJB

*********************************************************************/

  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif
  extern param_set_struct param_set;

  int i, j;

  /** Initialize model option flags **/

  // simulation modes
  options.AboveTreelineVeg      = -1;
  options.COMPUTE_TREELINE      = FALSE;
  options.CORRPREC              = FALSE;
  options.DIST_PRCP             = FALSE;
  options.FROZEN_SOIL           = FALSE;
  options.FULL_ENERGY           = FALSE;
  options.GRND_FLUX             = FALSE;
  options.MIN_WIND_SPEED        = 0.0;
  options.MOISTFRACT            = FALSE;
  options.Nlayer                = 2;
  options.Nnode                 = 3;
  options.NOFLUX                = FALSE;
  options.PREC_EXPT             = 0.6;
  options.QUICK_FLUX            = TRUE;
  options.SNOW_BAND             = 1;
  options.SNOW_STEP             = 1;
  // input options
  options.ALMA_INPUT            = FALSE;
  options.ARC_SOIL              = FALSE;
  options.BASEFLOW              = ARNO;
  options.GRID_DECIMAL          = 2;
  options.GLOBAL_LAI            = FALSE;
  options.JULY_TAVG_SUPPLIED    = FALSE;
  options.ROOT_ZONES            = MISSING;
  // state options
  options.INIT_STATE            = FALSE;
  options.SAVE_STATE            = FALSE;
  options.BINARY_STATE_FILE     = TRUE;
  // output options
  options.ALMA_OUTPUT           = FALSE;
  options.BINARY_OUTPUT         = FALSE;
  options.COMPRESS              = FALSE;
  options.Noutfiles             = 2;
  options.PRT_SNOW_BAND         = FALSE;

#if LINK_DEBUG 

  /** Initialize debugging control flags **/

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
#endif

  /** Initialize forcing file input controls **/

  param_set.TYPE[AIR_TEMP].SUPPLIED   = FALSE;
  param_set.TYPE[ALBEDO].SUPPLIED     = FALSE;
  param_set.TYPE[CRAINF].SUPPLIED     = FALSE;
  param_set.TYPE[CSNOWF].SUPPLIED     = FALSE;
  param_set.TYPE[DENSITY].SUPPLIED    = FALSE;
  param_set.TYPE[LONGWAVE].SUPPLIED   = FALSE;
  param_set.TYPE[LSRAINF].SUPPLIED    = FALSE;
  param_set.TYPE[LSSNOWF].SUPPLIED    = FALSE;
  param_set.TYPE[PREC].SUPPLIED       = FALSE;
  param_set.TYPE[PRESSURE].SUPPLIED   = FALSE;
  param_set.TYPE[QAIR].SUPPLIED       = FALSE;
  param_set.TYPE[RAINF].SUPPLIED      = FALSE;
  param_set.TYPE[REL_HUMID].SUPPLIED  = FALSE;
  param_set.TYPE[SHORTWAVE].SUPPLIED  = FALSE;
  param_set.TYPE[SNOWF].SUPPLIED      = FALSE;
  param_set.TYPE[TMAX].SUPPLIED       = FALSE;
  param_set.TYPE[TMIN].SUPPLIED       = FALSE;
  param_set.TYPE[TSKC].SUPPLIED       = FALSE;
  param_set.TYPE[VP].SUPPLIED         = FALSE;
  param_set.TYPE[WIND].SUPPLIED       = FALSE;
  param_set.TYPE[WIND_E].SUPPLIED     = FALSE;
  param_set.TYPE[WIND_N].SUPPLIED     = FALSE;
  param_set.TYPE[SKIP].SUPPLIED       = FALSE;
  for(i=0;i<2;i++) {
    param_set.FORCE_DT[i] = MISSING;
    param_set.N_TYPES[i] = MISSING;
    param_set.FORCE_FORMAT[i] = MISSING;
    for(j=0;j<N_FORCING_TYPES;j++) param_set.FORCE_INDEX[i][j] = MISSING;
  }

}
