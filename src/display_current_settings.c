#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void display_current_settings(int                 mode,
                              filenames_struct    *names,
                              global_param_struct *global) 
/**********************************************************************
  display_current_settings	Ted Bohn			2003

  This routine displays the current settings of options defined in
  user_def.h and the global parameter file.

  NOTE: This file must be kept in sync with any additions, removals,
  or modifications to names of parameters in user_def.h or get_global_param.c.

  Modifications:
  2005-03-08 Added EQUAL_AREA option.				TJB
  2005-03-24 Added ALMA_OUTPUT option.				TJB
  2005-04-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW.	TJB
  2005-11-29 SAVE_STATE is now set in global param file         GCT
  2006-09-13 Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option. TJB/GCT
  2006-Sep-23 Implemented flexible output configuration
              and aggregation of output variables.		TJB
  2006-Oct-10 Moved printing of soil_dir inside if{} block.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to names->statefile. TJB
  2006-Nov-07 Removed LAKE_MODEL option.			TJB
  2007-Jan-03 Added ALMA_INPUT option.				TJB
  2007-Jan-15 Added PRT_HEADER option.				TJB
  2007-Apr-24 Added IMPLICIT option.				JCA
  2007-Apr-24 Added EXP_TRANS option.				JCA
  2007-Aug-08 Added EXCESS_ICE option.				JCA
  2008-Apr-21 Added SNOW_ALBEDO option.				KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.			TJB
  2009-Jan-12 Added COMPUTE_TREELINE and JULY_TAVG_SUPPLIED options.	TJB
  2009-Jan-16 Added AERO_RESIST_CANSNOW option.			TJB
  2009-May-17 Added AR_406_LS to AERO_RESIST_CANSNOW.		TJB
  2009-May-18 Added PLAPSE option.				TJB
  2009-May-20 Added GRND_FLUX_TYPE option.			TJB
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.	TJB
  2009-Sep-19 Moved TFALLBACK to its own separate option.	TJB
  2009-Nov-15 Redirected output to stderr.			TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.TJB
  2011-May-31 Removed GRND_FLUX option.				TJB
  2011-Jun-03 Added ORGANIC_FRACT option.			TJB
  2011-Nov-04 Added options for accessing new forcing estimation
	      features.						TJB
  2012-Jan-16 Removed LINK_DEBUG code				BN
  2012-Jan-28 Removed AR_COMBO and GF_FULL.			TJB
  2014-Mar-24 Removed ARC_SOIL option       BN
**********************************************************************/
{

  extern char *version;
  extern option_struct options;
  extern param_set_struct param_set;

  int file_num;

  if (mode == DISP_VERSION) {
    fprintf(stderr,"***** VIC Version %s *****\n",version);
    return;
  }
  else {
    fprintf(stderr,"\n***** VIC Version %s - Current Model Settings *****\n",version);
  }

  fprintf(stderr,"\n");
  fprintf(stderr,"COMPILE-TIME OPTIONS (set in user_def.h)\n");
  fprintf(stderr,"----------------------------------------\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Output to Screen:\n");
#if OUTPUT_FORCE_STATS
  fprintf(stderr,"OUTPUT_FORCE_STATS\tTRUE\n");
#else
  fprintf(stderr,"OUTPUT_FORCE_STATS\tFALSE\n");
#endif
#if VERBOSE
  fprintf(stderr,"VERBOSE\t\t\tTRUE\n");
#else
  fprintf(stderr,"VERBOSE\t\t\tFALSE\n");
#endif

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Files:\n");
#if NO_REWIND
  fprintf(stderr,"NO_REWIND\t\tTRUE\n");
#else
  fprintf(stderr,"NO_REWIND\t\tFALSE\n");
#endif

  fprintf(stderr,"\n");
  fprintf(stderr,"Output Files:\n");
#if OUTPUT_FORCE
  fprintf(stderr,"OUTPUT_FORCE\t\tTRUE\n");
#else
  fprintf(stderr,"OUTPUT_FORCE\t\tFALSE\n");
#endif

  fprintf(stderr,"\n");
  fprintf(stderr,"Simulation Parameters:\n");
#if CLOSE_ENERGY
  fprintf(stderr,"CLOSE_ENERGY\t\tTRUE\n");
#else
  fprintf(stderr,"CLOSE_ENERGY\t\tFALSE\n");
#endif
#if LOW_RES_MOIST
  fprintf(stderr,"LOW_RES_MOIST\t\tTRUE\n");
#else
  fprintf(stderr,"LOW_RES_MOIST\t\tFALSE\n");
#endif
#if QUICK_FS
  fprintf(stderr,"QUICK_FS\t\tTRUE\n");
  fprintf(stderr,"QUICK_FS_TEMPS\t%d\n",QUICK_FS_TEMPS);
#else
  fprintf(stderr,"QUICK_FS\t\tFALSE\n");
#endif
#if SPATIAL_FROST
  fprintf(stderr,"SPATIAL_FROST\t\tTRUE\n");
  fprintf(stderr,"FROST_SUBAREAS\t\t%d\n",FROST_SUBAREAS);
#else
  fprintf(stderr,"SPATIAL_FROST\t\tFALSE\n");
#endif
#if SPATIAL_SNOW
  fprintf(stderr,"SPATIAL_SNOW\t\tTRUE\n");
#else
  fprintf(stderr,"SPATIAL_SNOW\t\tFALSE\n");
#endif
#if EXCESS_ICE
  fprintf(stderr,"EXCESS_ICE\t\tTRUE\n");
#else
  fprintf(stderr,"EXCESS_ICE\t\tFALSE\n");
#endif

  fprintf(stderr,"\n");
  fprintf(stderr,"Maximum Array Sizes:\n");
  fprintf(stderr,"MAX_BANDS\t\t%2d\n",MAX_BANDS);
  fprintf(stderr,"MAX_FRONTS\t\t%2d\n",MAX_FRONTS);
  fprintf(stderr,"MAX_LAKE_NODES\t\t%2d\n",MAX_LAKE_NODES);
  fprintf(stderr,"MAX_LAYERS\t\t%2d\n",MAX_LAYERS);
  fprintf(stderr,"MAX_NODES\t\t%2d\n",MAX_NODES);
  fprintf(stderr,"MAX_VEG\t\t\t%2d\n",MAX_VEG);
  fprintf(stderr,"\n");
  fprintf(stderr,"Snow Constants:\n");
  fprintf(stderr,"NEW_SNOW_ALB\t\t%f\n",NEW_SNOW_ALB);
  fprintf(stderr,"SNOW_ALB_ACCUM_A\t%f\n",SNOW_ALB_ACCUM_A);
  fprintf(stderr,"SNOW_ALB_ACCUM_B\t%f\n",SNOW_ALB_ACCUM_B);
  fprintf(stderr,"SNOW_ALB_THAW_A\t\t%f\n",SNOW_ALB_THAW_A);
  fprintf(stderr,"SNOW_ALB_THAW_B\t\t%f\n",SNOW_ALB_THAW_B);
  fprintf(stderr,"TraceSnow\t\t%f\n",TraceSnow);
  fprintf(stderr,"\n");
  fprintf(stderr,"Other Constants:\n");
  fprintf(stderr,"LAI_WATER_FACTOR\t%f\n",LAI_WATER_FACTOR);
  fprintf(stderr,"LWAVE_COR\t\t%f\n",LWAVE_COR);
  fprintf(stderr,"MAXIT_FE\t\t%2d\n",MAXIT_FE);

  if (mode == DISP_COMPILE_TIME) {
    return;
  }

  fprintf(stderr,"\n");
  fprintf(stderr,"RUN-TIME OPTIONS (set in global parameter file)\n");
  fprintf(stderr,"-----------------------------------------------\n");

  fprintf(stderr,"Simulation Dimensions:\n");
  fprintf(stderr,"NLAYER\t\t\t%d\n",options.Nlayer);
  if ( options.EQUAL_AREA ) {
    fprintf(stderr,"EQUAL_AREA\t\tTRUE\n");
  }
  else {
    fprintf(stderr,"EQUAL_AREA\t\tFALSE\n");
  }
  fprintf(stderr,"RESOLUTION\t\t%f\n",global->resolution);
  fprintf(stderr,"TIME_STEP\t\t%d\n",global->dt);
  fprintf(stderr,"SNOW_STEP\t\t%d\n",options.SNOW_STEP);
  fprintf(stderr,"STARTYEAR\t\t%d\n",global->startyear);
  fprintf(stderr,"STARTMONTH\t\t%d\n",global->startmonth);
  fprintf(stderr,"STARTDAY\t\t%d\n",global->startday);
  fprintf(stderr,"STARTHOUR\t\t%d\n",global->starthour);
  if ( global->nrecs > 0 )
    fprintf(stderr,"NRECS\t\t%d\n",global->nrecs);
  else {
    fprintf(stderr,"ENDYEAR\t\t\t%d\n",global->endyear);
    fprintf(stderr,"ENDMONTH\t\t%d\n",global->endmonth);
    fprintf(stderr,"ENDDAY\t\t\t%d\n",global->endday);
  }

  fprintf(stderr,"\n");
  fprintf(stderr,"Simulation Parameters:\n");
  if (options.AERO_RESIST_CANSNOW == AR_406)
    fprintf(stderr,"AERO_RESIST_CANSNOW\t\tAR_406\n");
  else if (options.AERO_RESIST_CANSNOW == AR_406_LS)
    fprintf(stderr,"AERO_RESIST_CANSNOW\t\tAR_406_LS\n");
  else if (options.AERO_RESIST_CANSNOW == AR_406_FULL)
    fprintf(stderr,"AERO_RESIST_CANSNOW\t\tAR_406_FULL\n");
  else if (options.AERO_RESIST_CANSNOW == AR_410)
    fprintf(stderr,"AERO_RESIST_CANSNOW\t\tAR_410\n");
  if (options.BLOWING)
    fprintf(stderr,"BLOWING\t\t\tTRUE\n");
  else
    fprintf(stderr,"BLOWING\t\t\tFALSE\n");
  if (options.COMPUTE_TREELINE)
    fprintf(stderr,"COMPUTE_TREELINE\t\tTRUE\n");
  else
    fprintf(stderr,"COMPUTE_TREELINE\t\tFALSE\n");
  if (options.CONTINUEONERROR == TRUE)
    fprintf(stderr,"CONTINUEONERROR\t\tTRUE\n");
  else
    fprintf(stderr,"CONTINUEONERROR\t\tFALSE\n");
  if (options.CORRPREC)
    fprintf(stderr,"CORRPREC\t\tTRUE\n");
  else
    fprintf(stderr,"CORRPREC\t\tFALSE\n");
  if (options.DIST_PRCP)
    fprintf(stderr,"DIST_PRCP\t\tTRUE\n");
  else
    fprintf(stderr,"DIST_PRCP\t\tFALSE\n");
  if (options.EXP_TRANS)
    fprintf(stderr,"EXP_TRANS\t\tTRUE\n");
  else
    fprintf(stderr,"EXP_TRANS\t\tFALSE\n");
  if (options.FROZEN_SOIL)
    fprintf(stderr,"FROZEN_SOIL\t\tTRUE\n");
  else
    fprintf(stderr,"FROZEN_SOIL\t\tFALSE\n");
  if (options.FULL_ENERGY)
    fprintf(stderr,"FULL_ENERGY\t\tTRUE\n");
  else
    fprintf(stderr,"FULL_ENERGY\t\tFALSE\n");
  if (options.GRND_FLUX_TYPE == GF_406)
    fprintf(stderr,"GRND_FLUX_TYPE\t\tGF_406\n");
  else if (options.GRND_FLUX_TYPE == GF_410)
    fprintf(stderr,"GRND_FLUX_TYPE\t\tGF_410\n");
  if (options.LW_TYPE == LW_TVA)
    fprintf(stderr,"LW_TYPE\t\tLW_TVA\n");
  else if (options.LW_TYPE == LW_ANDERSON)
    fprintf(stderr,"LW_TYPE\t\tLW_ANDERSON\n");
  else if (options.LW_TYPE == LW_BRUTSAERT)
    fprintf(stderr,"LW_TYPE\t\tLW_BRUTSAERT\n");
  else if (options.LW_TYPE == LW_SATTERLUND)
    fprintf(stderr,"LW_TYPE\t\tLW_SATTERLUND\n");
  else if (options.LW_TYPE == LW_IDSO)
    fprintf(stderr,"LW_TYPE\t\tLW_IDSO\n");
  else if (options.LW_TYPE == LW_PRATA)
    fprintf(stderr,"LW_TYPE\t\tLW_PRATA\n");
  if (options.LW_CLOUD == LW_CLOUD_DEARDORFF)
    fprintf(stderr,"LW_CLOUD\t\tLW_CLOUD_DEARDORFF\n");
  else
    fprintf(stderr,"LW_CLOUD\t\tLW_CLOUD_BRAS\n");
  if (options.IMPLICIT)
    fprintf(stderr,"IMPLICIT\t\tTRUE\n");
  else
    fprintf(stderr,"IMPLICIT\t\tFALSE\n");
  if (options.NOFLUX)
    fprintf(stderr,"NOFLUX\t\t\tTRUE\n");
  else
    fprintf(stderr,"NOFLUX\t\t\tFALSE\n");
  if (options.MTCLIM_SWE_CORR)
    fprintf(stderr,"MTCLIM_SWE_CORR\t\tTRUE\n");
  else
    fprintf(stderr,"MTCLIM_SWE_CORR\t\tFALSE\n");
  if (options.PLAPSE)
    fprintf(stderr,"PLAPSE\t\tTRUE\n");
  else
    fprintf(stderr,"PLAPSE\t\tFALSE\n");
  if (options.QUICK_FLUX)
    fprintf(stderr,"QUICK_FLUX\t\tTRUE\n");
  else
    fprintf(stderr,"QUICK_FLUX\t\tFALSE\n");
  if (options.QUICK_SOLVE)
    fprintf(stderr,"QUICK_SOLVE\t\tTRUE\n");
  else
    fprintf(stderr,"QUICK_SOLVE\t\tFALSE\n");
  if (options.SNOW_ALBEDO == USACE)
    fprintf(stderr,"SNOW_ALBEDO\t\tUSACE\n");
  else if (options.SNOW_ALBEDO == SUN1999)
    fprintf(stderr,"SNOW_ALBEDO\t\tSUN1999\n");
  if (options.SNOW_DENSITY == DENS_BRAS)
    fprintf(stderr,"SNOW_DENSITY\t\tDENS_BRAS\n");
  else if (options.SNOW_DENSITY == DENS_SNTHRM)
    fprintf(stderr,"SNOW_DENSITY\t\tDENS_SNTHRM\n");
  fprintf(stderr,"SW_PREC_THRESH\t\t%f\n",options.SW_PREC_THRESH);
  if (options.TFALLBACK == TRUE)
    fprintf(stderr,"TFALLBACK\t\tTRUE\n");
  else
    fprintf(stderr,"TFALLBACK\t\tFALSE\n");
  if (options.VP_INTERP == TRUE)
    fprintf(stderr,"VP_INTERP\t\tTRUE\n");
  else
    fprintf(stderr,"VP_INTERP\t\tFALSE\n");
  if (options.VP_ITER == VP_ITER_NONE)
    fprintf(stderr,"VP_ITER\t\tVP_ITER_NONE\n");
  else if (options.VP_ITER == VP_ITER_ALWAYS)
    fprintf(stderr,"VP_ITER\t\tVP_ITER_ALWAYS\n");
  else if (options.VP_ITER == VP_ITER_ANNUAL)
    fprintf(stderr,"VP_ITER\t\tVP_ITER_ANNUAL\n");
  else if (options.VP_ITER == VP_ITER_CONVERGE)
    fprintf(stderr,"VP_ITER\t\tVP_ITER_CONVERGE\n");
  fprintf(stderr,"PREC_EXPT\t\t%f\n",options.PREC_EXPT);
  fprintf(stderr,"WIND_H\t\t\t%f\n",global->wind_h);
  fprintf(stderr,"MEASURE_H\t\t%f\n",global->measure_h);
  fprintf(stderr,"NODES\t\t\t%d\n",options.Nnode);
  fprintf(stderr,"MIN_RAIN_TEMP\t\t%f\n",global->MIN_RAIN_TEMP);
  fprintf(stderr,"MAX_SNOW_TEMP\t\t%f\n",global->MAX_SNOW_TEMP);
  fprintf(stderr,"MIN_WIND_SPEED\t\t%f\n",options.MIN_WIND_SPEED);

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Forcing Data:\n");
  for (file_num=0; file_num<2; file_num++) {
    if (global->forceyear[file_num] > 0) {
      fprintf(stderr,"Forcing File %d:\t\t%s*\n",file_num+1,names->f_path_pfx[file_num]);
      fprintf(stderr,"FORCEYEAR\t\t%d\n",global->forceyear[file_num]);
      fprintf(stderr,"FORCEMONTH\t\t%d\n",global->forcemonth[file_num]);
      fprintf(stderr,"FORCEDAY\t\t%d\n",global->forceday[file_num]);
      fprintf(stderr,"FORCEHOUR\t\t%d\n",global->forcehour[file_num]);
      fprintf(stderr,"N_TYPES\t\t\t%d\n",param_set.N_TYPES[file_num]);
      fprintf(stderr,"FORCE_DT\t\t%d\n",param_set.FORCE_DT[file_num]);
      if (param_set.FORCE_ENDIAN[file_num] == LITTLE)
        fprintf(stderr,"FORCE_ENDIAN\t\tLITTLE\n");
      else
        fprintf(stderr,"FORCE_ENDIAN\t\tBIG\n");
      if (param_set.FORCE_FORMAT[file_num] == BINARY)
        fprintf(stderr,"FORCE_FORMAT\t\tBINARY\n");
      else
        fprintf(stderr,"FORCE_FORMAT\t\tASCII\n");
    }
  }
  fprintf(stderr,"GRID_DECIMAL\t\t%d\n",options.GRID_DECIMAL);
  if (options.ALMA_INPUT)
    fprintf(stderr,"ALMA_INPUT\t\tTRUE\n");
  else
    fprintf(stderr,"ALMA_INPUT\t\tFALSE\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Soil Data:\n");
  fprintf(stderr,"Soil file\t\t%s\n",names->soil);
  if (options.BASEFLOW == ARNO)
    fprintf(stderr,"BASEFLOW\t\tARNO\n");
  else if (options.BASEFLOW == NIJSSEN2001)
    fprintf(stderr,"BASEFLOW\t\tNIJSSEN2001\n");
  if (options.JULY_TAVG_SUPPLIED)
    fprintf(stderr,"JULY_TAVG_SUPPLIED\t\tTRUE\n");
  else
    fprintf(stderr,"JULY_TAVG_SUPPLIED\t\tFALSE\n");
  if (options.ORGANIC_FRACT)
    fprintf(stderr,"ORGANIC_FRACT\t\tTRUE\n");
  else
    fprintf(stderr,"ORGANIC_FRACT\t\tFALSE\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Veg Data:\n");
  fprintf(stderr,"Veg library file\t%s\n",names->veglib);
  fprintf(stderr,"Veg param file\t\t%s\n",names->veg);
  fprintf(stderr,"ROOT_ZONES\t\t%d\n",options.ROOT_ZONES);
  if (options.VEGPARAM_LAI)
    fprintf(stderr,"VEGPARAM_LAI\t\tTRUE\n");
  else
    fprintf(stderr,"VEGPARAM_LAI\t\tFALSE\n");
  if (options.LAI_SRC == LAI_FROM_VEGPARAM)
    fprintf(stderr,"LAI_SRC\t\tLAI_FROM_VEGPARAM\n");
  else if (options.LAI_SRC == LAI_FROM_VEGLIB)
    fprintf(stderr,"LAI_SRC\t\tLAI_FROM_VEGLIB\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Elevation Data:\n");
  if (options.SNOW_BAND > 1)
    fprintf(stderr,"SNOW_BAND\t\t%d\t%s\n",options.SNOW_BAND,names->snowband);
  else if (options.SNOW_BAND == 1)
    fprintf(stderr,"SNOW_BAND\t\t%d\t(no input file needed for SNOW_BAND=1)\n",options.SNOW_BAND);
  else
    fprintf(stderr,"SNOW_BAND\t\t%d\n",options.SNOW_BAND);

  fprintf(stderr,"\n");
  fprintf(stderr,"Input Lake Data:\n");
  if (options.LAKES)
    fprintf(stderr,"LAKES\t\tTRUE\t%s\n",names->lakeparam);
  else
    fprintf(stderr,"LAKES\t\tFALSE\n");
  if (options.LAKE_PROFILE)
    fprintf(stderr,"LAKE_PROFILE\t\tTRUE\n");
  else
    fprintf(stderr,"LAKE_PROFILE\t\tFALSE\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Input State File:\n");
  if (options.INIT_STATE) {
    fprintf(stderr,"INIT_STATE\t\tTRUE\t%s\n",names->init_state);
    if (options.BINARY_STATE_FILE)
      fprintf(stderr,"BINARY_STATE_FILE\tTRUE\n");
    else
      fprintf(stderr,"BINARY_STATE_FILE\tFALSE\n");
  }
  else
    fprintf(stderr,"INIT_STATE\t\tFALSE\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"Output State File:\n");
  if (options.SAVE_STATE) {
    fprintf(stderr,"SAVE_STATE\t\tTRUE\n");
    fprintf(stderr,"STATENAME\t\t%s\n",names->statefile);
    fprintf(stderr,"STATEYEAR\t\t%d\n",global->stateyear);
    fprintf(stderr,"STATEMONTH\t\t%d\n",global->statemonth);
    fprintf(stderr,"STATEDAY\t\t%d\n",global->stateday);
    if (options.BINARY_STATE_FILE)
      fprintf(stderr,"BINARY_STATE_FILE\tTRUE\n");
    else
      fprintf(stderr,"BINARY_STATE_FILE\tFALSE\n");
  }
  else {
    fprintf(stderr,"SAVE_STATE\t\tFALSE\n");
  }

  fprintf(stderr,"\n");
  fprintf(stderr,"Output Data:\n");
  fprintf(stderr,"Result dir:\t\t%s\n",names->result_dir);
  fprintf(stderr,"OUT_STEP\t\t%d\n",global->out_dt);
  if (options.ALMA_OUTPUT)
    fprintf(stderr,"ALMA_OUTPUT\t\tTRUE\n");
  else
    fprintf(stderr,"ALMA_OUTPUT\t\tFALSE\n");
  if (options.BINARY_OUTPUT)
    fprintf(stderr,"BINARY_OUTPUT\t\tTRUE\n");
  else
    fprintf(stderr,"BINARY_OUTPUT\t\tFALSE\n");
  if (options.COMPRESS)
    fprintf(stderr,"COMPRESS\t\tTRUE\n");
  else
    fprintf(stderr,"COMPRESS\t\tFALSE\n");
  if (options.MOISTFRACT)
    fprintf(stderr,"MOISTFRACT\t\tTRUE\n");
  else
    fprintf(stderr,"MOISTFRACT\t\tFALSE\n");
  if (options.PRT_HEADER)
    fprintf(stderr,"PRT_HEADER\t\tTRUE\n");
  else
    fprintf(stderr,"PRT_HEADER\t\tFALSE\n");
  if (options.PRT_SNOW_BAND)
    fprintf(stderr,"PRT_SNOW_BAND\t\tTRUE\n");
  else
    fprintf(stderr,"PRT_SNOW_BAND\t\tFALSE\n");
  fprintf(stderr,"SKIPYEAR\t\t%d\n",global->skipyear);
  fprintf(stderr,"\n");

}
