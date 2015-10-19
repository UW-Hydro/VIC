/**********************************************************************
                        Global Variables

  NOTE: This file exists because global variables that are shared among
        files via the "extern" statement must be initially declared
        (without the word "extern") ONLY once.  Currently, vicNl_def.h
        is included (via vicNl.h) in every .c file, meaning that any
        declarations in vicNl_def.h end up happening multiple times
        (once per .c file).  Thus, these "extern" variables cannot be
        declared in vicNl_def.h.  This is not a problem for #define
        statements and typedef statements, which is what vicNl_def.h
        is primarily composed of.

  $Id$

  2003-Oct-29 Added version string and removed unused options from
	      optstring.						TJB
  2009-Jun-09 Added definitions of reference landcover types, used
	      mainly for pot_evap computations but also defines the
	      characteristics of bare soil.				TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2014-May-20 Added ref_veg_vegcover.					TJB
**********************************************************************/
char *version = "4.2.a.irrigation";
char *optstring = "g:vo";
int flag;

global_param_struct global_param;
veg_lib_struct *veg_lib;
option_struct options;
Error_struct Error;
param_set_struct param_set;

  /**************************************************************************
    Define some reference landcover types that always exist regardless
    of the contents of the library (mainly for potential evap calculations):
    Non-natural:
      satsoil = saturated bare soil
      h2osurf = open water surface (deep enough to have albedo of 0.08)
      short   = short reference crop (grass)
      tall    = tall reference crop (alfalfa)
    Natural:
      natveg  = current vegetation
      vegnocr = current vegetation with canopy resistance set to 0
    NOTE: these are external variables, declared in vicNl_def.h.
    NOTE2: bare soil roughness and displacement will be overwritten by the
           values found in the soil parameter file; bare soil wind_h will
	   be overwritten by the value specified in the global param file.
  **************************************************************************/

  /* One element for each non-natural PET type */
  char   ref_veg_over[]        = { 0, 0, 0, 0 };
  double ref_veg_rarc[]        = { 0.0, 0.0, 25, 25 };
  double ref_veg_rmin[]        = { 0.0, 0.0, 100, 100 };
  double ref_veg_lai[]         = { 1.0, 1.0, 2.88, 4.45 };
  double ref_veg_albedo[]      = { BARE_SOIL_ALBEDO, H2O_SURF_ALBEDO, 0.23, 0.23 };
  double ref_veg_vegcover[]      = { MIN_VEGCOVER, MIN_VEGCOVER, 1.00, 1.00 };
  double ref_veg_crop_frac[]   = { 0.0, 0.0, 1.00, 1.00 };
  double ref_veg_rough[]       = { 0.001, 0.001, 0.0148, 0.0615 };
  double ref_veg_displ[]       = { 0.0054, 0.0054, 0.08, 0.3333 };
  double ref_veg_wind_h[]      = { 10.0, 10.0, 10.0, 10.0 };
  double ref_veg_RGL[]         = { 0.0, 0.0, 100, 100 };
  double ref_veg_rad_atten[]   = { 0.0, 0.0, 0.0, 0.0 };
  double ref_veg_wind_atten[]  = { 0.0, 0.0, 0.0, 0.0 };
  double ref_veg_trunk_ratio[] = { 0.0, 0.0, 0.0, 0.0 };
  /* One element for each PET type (non-natural or natural) */
  char ref_veg_ref_crop[] = { FALSE, FALSE, TRUE, TRUE, FALSE, FALSE };

