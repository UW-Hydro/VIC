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

  07-May-04 Replaced "Release" with "Version".  Moved COMPUTE_TREELINE
	    lower down in output display.			TJB
  12-May-04 Changed formatting from %s to %d in STATE* print
	    statements.						TJB
  2005-11-21 (Port from 4.1.0) Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW. GCT
  2006-Jan-22    Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option. TJB
  2006-Sep-11 Implemented flexible output configuration; removed
              OPTIMIZE and LDAS_OUTPUT options.			TJB
  2006-Sep-14 Implemented ALMA-compliant input and output.	TJB
  2006-Sep-18 Implemented aggregation of output variables.	TJB
  2006-Oct-26 Moved printing of soil_dir inside if{} block.	TJB
  2006-Oct-26 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to names->statefile. TJB

**********************************************************************/
{

  extern char *version;
  extern option_struct options;
  extern param_set_struct param_set;

  int file_num;
  char typestr[20];

  if (mode == DISP_VERSION) {
    fprintf(stdout,"***** VIC Version %s *****\n",version);
    return;
  }
  else {
    fprintf(stdout,"\n***** VIC Version %s - Current Model Settings *****\n",version);
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"COMPILE-TIME OPTIONS (set in user_def.h)\n");
  fprintf(stdout,"----------------------------------------\n");

  fprintf(stdout,"\n");
  fprintf(stdout,"Output to Screen:\n");
#if VERBOSE
  fprintf(stdout,"VERBOSE\t\t\tTRUE\n");
#else
  fprintf(stdout,"VERBOSE\t\t\tFALSE\n");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Input Files:\n");
#if NO_REWIND
  fprintf(stdout,"NO_REWIND\t\tTRUE\n");
#else
  fprintf(stdout,"NO_REWIND\t\tFALSE\n");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Output Files:\n");
#if LINK_DEBUG
  fprintf(stdout,"LINK_DEBUG\t\tTRUE\n");
#else
  fprintf(stdout,"LINK_DEBUG\t\tFALSE\n");
#endif
#if OUTPUT_FORCE
  fprintf(stdout,"OUTPUT_FORCE\t\tTRUE\n");
#else
  fprintf(stdout,"OUTPUT_FORCE\t\tFALSE\n");
#endif
  fprintf(stdout,"\n");
  fprintf(stdout,"Simulation Parameters:\n");
#if LOW_RES_MOIST
  fprintf(stdout,"LOW_RES_MOIST\t\tTRUE\n");
#else
  fprintf(stdout,"LOW_RES_MOIST\t\tFALSE\n");
#endif
#if QUICK_FS
  fprintf(stdout,"QUICK_FS\t\tTRUE\n");
  fprintf(stdout,"QUICK_FS_TEMPS\t%d\n",QUICK_FS_TEMPS);
#else
  fprintf(stdout,"QUICK_FS\t\tFALSE\n");
#endif

  fprintf(stdout,"\n");
  fprintf(stdout,"Maximum Array Sizes:\n");
  fprintf(stdout,"MAX_BANDS\t\t%2d\n",MAX_BANDS);
  fprintf(stdout,"MAX_FRONTS\t\t%2d\n",MAX_FRONTS);
  fprintf(stdout,"MAX_LAYERS\t\t%2d\n",MAX_LAYERS);
  fprintf(stdout,"MAX_NODES\t\t%2d\n",MAX_NODES);
  fprintf(stdout,"MAX_VEG\t\t\t%2d\n",MAX_VEG);
  fprintf(stdout,"\n");
  fprintf(stdout,"Snow Constants:\n");
  fprintf(stdout,"NEW_SNOW_ALB\t\t%f\n",NEW_SNOW_ALB);
  fprintf(stdout,"SNOW_ALB_ACCUM_A\t%f\n",SNOW_ALB_ACCUM_A);
  fprintf(stdout,"SNOW_ALB_ACCUM_B\t%f\n",SNOW_ALB_ACCUM_B);
  fprintf(stdout,"SNOW_ALB_THAW_A\t\t%f\n",SNOW_ALB_THAW_A);
  fprintf(stdout,"SNOW_ALB_THAW_B\t\t%f\n",SNOW_ALB_THAW_B);
  fprintf(stdout,"TraceSnow\t\t%f\n",TraceSnow);
  fprintf(stdout,"\n");
  fprintf(stdout,"Other Constants:\n");
  fprintf(stdout,"LAI_WATER_FACTOR\t%f\n",LAI_WATER_FACTOR);
  fprintf(stdout,"LWAVE_COR\t\t%f\n",LWAVE_COR);
  fprintf(stdout,"MAXIT_FE\t\t%2d\n",MAXIT_FE);

  if (mode == DISP_COMPILE_TIME) {
    return;
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"RUN-TIME OPTIONS (set in global parameter file)\n");
  fprintf(stdout,"-----------------------------------------------\n");

  fprintf(stdout,"Simulation Dimensions:\n");
  fprintf(stdout,"NLAYER\t\t\t%d\n",options.Nlayer);
  fprintf(stdout,"RESOLUTION\t\t%f\n",global->resolution);
  fprintf(stdout,"TIME_STEP\t\t%d\n",global->dt);
  fprintf(stdout,"SNOW_STEP\t\t%d\n",options.SNOW_STEP);
  fprintf(stdout,"STARTYEAR\t\t%d\n",global->startyear);
  fprintf(stdout,"STARTMONTH\t\t%d\n",global->startmonth);
  fprintf(stdout,"STARTDAY\t\t%d\n",global->startday);
  fprintf(stdout,"STARTHOUR\t\t%d\n",global->starthour);
  if ( global->nrecs > 0 )
    fprintf(stdout,"NRECS\t\t%d\n",global->nrecs);
  else {
    fprintf(stdout,"ENDYEAR\t\t\t%d\n",global->endyear);
    fprintf(stdout,"ENDMONTH\t\t%d\n",global->endmonth);
    fprintf(stdout,"ENDDAY\t\t\t%d\n",global->endday);
  }

  fprintf(stdout,"\n");
  fprintf(stdout,"Simulation Parameters:\n");
  if (options.FULL_ENERGY)
    fprintf(stdout,"FULL_ENERGY\t\tTRUE\n");
  else
    fprintf(stdout,"FULL_ENERGY\t\tFALSE\n");
  if (options.FROZEN_SOIL)
    fprintf(stdout,"FROZEN_SOIL\t\tTRUE\n");
  else
    fprintf(stdout,"FROZEN_SOIL\t\tFALSE\n");
  if (options.QUICK_FLUX)
    fprintf(stdout,"QUICK_FLUX\t\tTRUE\n");
  else
    fprintf(stdout,"QUICK_FLUX\t\tFALSE\n");
  if (options.GRND_FLUX)
    fprintf(stdout,"GRND_FLUX\t\tTRUE\n");
  else
    fprintf(stdout,"GRND_FLUX\t\tFALSE\n");
  if (options.NOFLUX)
    fprintf(stdout,"NOFLUX\t\t\tTRUE\n");
  else
    fprintf(stdout,"NOFLUX\t\t\tFALSE\n");
  if (options.CORRPREC)
    fprintf(stdout,"CORRPREC\t\tTRUE\n");
  else
    fprintf(stdout,"CORRPREC\t\tFALSE\n");
  if (options.DIST_PRCP)
    fprintf(stdout,"DIST_PRCP\t\tTRUE\n");
  else
    fprintf(stdout,"DIST_PRCP\t\tFALSE\n");
  if (options.MOISTFRACT)
    fprintf(stdout,"MOISTFRACT\t\tTRUE\n");
  else
    fprintf(stdout,"MOISTFRACT\t\tFALSE\n");
  if (options.COMPUTE_TREELINE)
    fprintf(stdout,"COMPUTE_TREELINE\t\tTRUE\n");
  else
    fprintf(stdout,"COMPUTE_TREELINE\t\tFALSE\n");
  fprintf(stdout,"PREC_EXPT\t\t%f\n",options.PREC_EXPT);
  fprintf(stdout,"WIND_H\t\t\t%f\n",global->wind_h);
  fprintf(stdout,"MEASURE_H\t\t%f\n",global->measure_h);
  fprintf(stdout,"NODES\t\t\t%d\n",options.Nnode);
  fprintf(stdout,"MIN_RAIN_TEMP\t\t%f\n",global->MIN_RAIN_TEMP);
  fprintf(stdout,"MAX_SNOW_TEMP\t\t%f\n",global->MAX_SNOW_TEMP);
  fprintf(stdout,"MIN_WIND_SPEED\t\t%f\n",options.MIN_WIND_SPEED);

  fprintf(stdout,"\n");
  fprintf(stdout,"Input Forcing Data:\n");
  for (file_num=0; file_num<2; file_num++) {
    if (global->forceyear[file_num] > 0) {
      fprintf(stdout,"Forcing File %d:\t\t%s*\n",file_num+1,names->f_path_pfx[file_num]);
      fprintf(stdout,"FORCEYEAR\t\t%d\n",global->forceyear[file_num]);
      fprintf(stdout,"FORCEMONTH\t\t%d\n",global->forcemonth[file_num]);
      fprintf(stdout,"FORCEDAY\t\t%d\n",global->forceday[file_num]);
      fprintf(stdout,"FORCEHOUR\t\t%d\n",global->forcehour[file_num]);
      fprintf(stdout,"N_TYPES\t\t\t%d\n",param_set.N_TYPES[file_num]);
      fprintf(stdout,"FORCE_DT\t\t%d\n",param_set.FORCE_DT[file_num]);
      if (param_set.FORCE_ENDIAN[file_num] == LITTLE)
        fprintf(stdout,"FORCE_ENDIAN\t\tLITTLE\n");
      else
        fprintf(stdout,"FORCE_ENDIAN\t\tBIG\n");
      if (param_set.FORCE_FORMAT[file_num] == BINARY)
        fprintf(stdout,"FORCE_FORMAT\t\tBINARY\n");
      else
        fprintf(stdout,"FORCE_FORMAT\t\tASCII\n");
    }
  }
  fprintf(stdout,"GRID_DECIMAL\t\t%d\n",options.GRID_DECIMAL);

  fprintf(stdout,"\n");
  fprintf(stdout,"Input Soil Data:\n");
  fprintf(stdout,"Soil file\t\t%s\n",names->soil);
  if (options.ARC_SOIL) {
    fprintf(stdout,"ARC_SOIL\t\tTRUE\n");
    fprintf(stdout,"Soil dir\t\t%s\n",names->soil_dir);
  }
  else
    fprintf(stdout,"ARC_SOIL\t\tFALSE\n");
  if (options.BASEFLOW == ARNO)
    fprintf(stdout,"BASEFLOW\t\tARNO\n");
  else if (options.BASEFLOW == NIJSSEN2001)
    fprintf(stdout,"BASEFLOW\t\tNIJSSEN2001\n");

  fprintf(stdout,"\n");
  fprintf(stdout,"Input Veg Data:\n");
  fprintf(stdout,"Veg param file\t\t%s\n",names->veg);
  fprintf(stdout,"ROOT_ZONES\t\t%d\n",options.ROOT_ZONES);
  fprintf(stdout,"Veg library file\t%s\n",names->veglib);
  if (options.GLOBAL_LAI)
    fprintf(stdout,"GLOBAL_LAI\t\tTRUE\n");
  else
    fprintf(stdout,"GLOBAL_LAI\t\tFALSE\n");

  fprintf(stdout,"\n");
  fprintf(stdout,"Input Elevation Data:\n");
  if (options.SNOW_BAND > 1)
    fprintf(stdout,"SNOW_BAND\t\t%d\t%s\n",options.SNOW_BAND,names->snowband);
  else if (options.SNOW_BAND == 1)
    fprintf(stdout,"SNOW_BAND\t\t%d\t(no input file needed for SNOW_BAND=1)\n",options.SNOW_BAND);
  else
    fprintf(stdout,"SNOW_BAND\t\t%d\n",options.SNOW_BAND);

  fprintf(stdout,"\n");
  fprintf(stdout,"Input State File:\n");
  if (options.INIT_STATE) {
    fprintf(stdout,"INIT_STATE\t\tTRUE\t%s\n",names->init_state);
    if (options.BINARY_STATE_FILE)
      fprintf(stdout,"BINARY_STATE_FILE\tTRUE\n");
    else
      fprintf(stdout,"BINARY_STATE_FILE\tFALSE\n");
  }
  else
    fprintf(stdout,"INIT_STATE\t\tFALSE\n");

  fprintf(stdout,"\n");
  fprintf(stdout,"Output State File:\n");
  if (options.SAVE_STATE ) {
    fprintf(stdout, "SAVE_STATE\t\tTRUE\n");
    fprintf(stdout,"STATENAME\t\t%s\n",names->statefile);
    fprintf(stdout,"STATEYEAR\t\t%d\n",global->stateyear);
    fprintf(stdout,"STATEMONTH\t\t%d\n",global->statemonth);
    fprintf(stdout,"STATEDAY\t\t%d\n",global->stateday);
    if (options.BINARY_STATE_FILE)
      fprintf(stdout,"BINARY_STATE_FILE\tTRUE\n");
    else
      fprintf(stdout,"BINARY_STATE_FILE\tFALSE\n");
  }
  else {
    fprintf(stdout, "SAVE_STATE\t\tFALSE\n");
  }
  fprintf(stdout,"\n");
  fprintf(stdout,"Output Data:\n");
  fprintf(stdout,"Result dir:\t\t%s\n",names->result_dir);
  fprintf(stdout,"OUT_STEP\t\t%d\n",global->out_dt);
  if (options.BINARY_OUTPUT)
    fprintf(stdout,"BINARY_OUTPUT\t\tTRUE\n");
  else
    fprintf(stdout,"BINARY_OUTPUT\t\tFALSE\n");
  if (options.COMPRESS)
    fprintf(stdout,"COMPRESS\t\tTRUE\n");
  else
    fprintf(stdout,"COMPRESS\t\tFALSE\n");
  fprintf(stdout,"SKIPYEAR\t\t%d\n",global->skipyear);
  if (options.PRT_SNOW_BAND)
    fprintf(stdout,"PRT_SNOW_BAND\t\tTRUE\n");
  else
    fprintf(stdout,"PRT_SNOW_BAND\t\tFALSE\n");
  if (options.ALMA_OUTPUT)
    fprintf(stdout,"ALMA_OUTPUT\t\tTRUE\n");
  else
    fprintf(stdout,"ALMA_OUTPUT\t\tFALSE\n");
  fprintf(stdout,"\n");

}
