/**********************************************************************
  This header file contains model parameters that can be modified by
  the user to control model performance.  When this file is modified 
  the model needs to be recompiled for the changes to take effect.

  NOTE (to those who add or remove parameters from VIC):
  Any additions or removals of parameters in this file must also
  be made in display_current_settings.c.

  $Id$

  Modifications:
  2006-Sep-23 Implemented flexible output configuration; removed the
              OPTIMIZE and LDAS_OUTPUT options.				TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2008-Feb-17 Moved constants related to snow albedo and trace snow to
	      snow.h.							TJB
  2009-Jan-12 Removed COMPUTE_TREELINE option (moved into options
	      struct).							TJB
  2011-Jan-04 Added MAX_ZWTVMOIST for storing array of soil moisture vs
	      water table position.					TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-26 Removed LWAVE_COR option.					TJB
  2013-Dec-26 Moved MAX_VEG and other array lengths to vicNl_def.h.	TJB
  2013-Dec-26 Removed OUTPUT_FORCE_STATS option.			TJB
  2013-Dec-26 Replaced LOW_RES_MOIST compile-time option (user_def.h) with
	      LOG_MATRIC run-time option (vicNl_def.h).			TJB
  2013-Dec-26 Removed CLOSE_ENERGY option (moved into options struct).	TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_SNOW to options_struct (vicNl_def.h).	TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct (vicNl_def.h).	TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct (vicNl_def.h).	TJB
  2013-Dec-28 Removed NO_REWIND option.					TJB
**********************************************************************/

/***** If TRUE include all model messages to stdout, and stderr *****/
#define VERBOSE TRUE

