#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_atmos(atmos_data_struct        *atmos,
                      dmy_struct               *dmy,
		      FILE                    **infile,
#if OUTPUT_FORCE
		      soil_con_struct          *soil_con,
                      out_data_file_struct     *out_data_files,
                      out_data_struct          *out_data)
#else /* OUTPUT_FORCE */
		      soil_con_struct          *soil_con)
#endif /* OUTPUT_FORCE */
/**********************************************************************
  initialize_atmos	Keith Cherkauer		February 3, 1997

  This routine initializes atmospheric variables for both the model
  time step, and the time step used by the snow algorithm (if different).
  Air temperature is estimated using MTCLIM (see routine for reference),
  atmospheric moisture is estimated using Kimball's algorithm (see 
  routine for reference), and radiation is estimated using Bras's algorithms
  (see routines for reference).

  WARNING: This subroutine is site specific.  Location parameters
    must be changed before compilation.

  UNITS: mks
	energy - W/m^2

  Modifications:
  11-18-98  Removed variable array yearly_epot, since yearly potential
            evaporation is no longer used for estimating the dew
            point temperature from daily minimum temperature.   KAC
  11-25-98  Added second check to make sure that the difference 
            between tmax and tmin is positive, after being reset
            when it was equal to 0.                        DAG, EFW
  12-1-98   Changed relative humidity computations so that they 
            use air temperature for the time step, instead of average
            daily temperature.  This allows relative humidity to
            change during the day, when the time step is less than
            daily.                                              KAC
  8-19-99   MIN_TDEW was added to prevent the dew point temperature
            estimated by Kimball's equations from becoming so low
            that svp() fails.							Bart
  9-4-99    Code was largely rewritten to change make use of the MTCLIM
            meteorological preprocessor which estimates sub-daily 
	    met forcings for all time steps.  The atmos_data_struct was
	    also reconfigured so that it has a new record for each
	    model time step, but stores sub-time step forcing data
	    (that might be needed for the snow model) within each
	    record, eliminating the on the fly estimations used in
	    previous versions of the model.					Bart and Greg
  01-17-01  Pressure and vapor pressure read from a forcing file are
            converted from kPa to Pa.  This preserves the original
            format of the forcing files (where pressure was supposed 
            to be in kPa, but allows VIC to use Pa internally, eliminating
            the need to convert to Pa every time it is used.			KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.			KAC
  04-Oct-04 Changed logic to allow VP to be supplied without
	    SHORTWAVE.								TJB
  2005-Mar-24 Modified to handle ALMA forcing variables.			TJB
  2005-Apr-30 Fixed typo in QAIR calculation.					TJB
  2005-May-01 Added logic for CSNOWF and LSSNOWF.				TJB
  2005-May-02 Added logic for WIND_E and WIND_N.				TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
	      out_data and out_data_files structures.				TJB
  2006-Dec-20 Replaced 1000.0 with kPa2Pa in pressure conversion.		TJB
  2006-Dec-29 Added REL_HUMID to the list of supported met input variables.	TJB
  2007-Jan-02 Added ALMA_INPUT option; removed TAIR and PSURF from list of
	      supported met input variables.					TJB
  2008-Jan-25 Fixed conditions under which net longwave replaces incoming
	      longwave in atmos[rec].longwave[NR].  Previously, net longwave
	      was stored if SNOW_STEP != global.dt.  Now, net longwave is
	      stored if options.FULL_ENERGY and options.FROZEN_SOIL are both
	      FALSE, i.e. for a water balance mode run.				TJB
  2009-Jan-12 Modified to pass avgJulyAirTemp argument to
	      compute_treeline(). 						TJB
  2009-May-18 Added options.PLAPSE, which when TRUE changes pressure
	      calculation to be a function of elevation and air temperature
	      (as opposed to a constant 95.5 kPa, as it was previously).
	      Made similar change to density calculation.			TJB
  2009-Jun-10 Fixed incorrect handling of cases when incoming longwave and
	      shortwave radiation are supplied.					TJB
  2009-Jul-26 Removed the special logic for the water balance mode, in
	      which net longwave is stored in the "longwave" variable.		TJB
  2009-Oct-13 Removed condition if(options.SNOW_BAND) for call to
	      compute_treeline(), since options.SNOW_BAND is always > 0.	TJB
  2010-Mar-31 Added RUNOFF_IN.							TJB
  2010-Apr-28 Removed individual soil_con variables from argument list and
	      replaced with *soil_con.						TJB
  2010-Sep-24 Renamed RUNOFF_IN to CHANNEL_IN.					TJB

**********************************************************************/
{
  extern option_struct       options;
  extern param_set_struct    param_set;
  extern global_param_struct global_param;
  extern int                 NR, NF;

  int     i;
  int     j;
  int     band;
  int     day;
  int     hour;
  int     rec;
  int     step;
  int     idx;
  int    *tmaxhour;
  int    *tminhour;
  double  deltat;
  double  cell_area;
  double  theta_l;
  double  theta_s;
  double  phi;
  double  elevation;
  double  annual_prec;
  double  wind_h;
  double  roughness;
  double  avgJulyAirTemp;
  double *Tfactor;
  char   *AboveTreeLine;
  double  min_Tfactor;
  double  shortwave;
  double  svp_tair;
  double *hourlyrad;
  double *prec;
  double *tmax;
  double *tmin;
  double *tair;
  double *tskc;
  double *vp;
  double  min, max;
  double  rainonly;
  int     Ndays;
  int     stepspday;
  double  sum;
  double **forcing_data;
  int     type;
  double  air_temp;
  double  factor;

  wind_h = global_param.wind_h;
  theta_l = (double)soil_con->time_zone_lng;
  theta_s = (double)soil_con->lng;
  phi = soil_con->lat;
  elevation = soil_con->elevation;
  annual_prec = soil_con->annual_prec;
  roughness = soil_con->rough;
  cell_area = soil_con->cell_area;
  avgJulyAirTemp = soil_con->avgJulyAirTemp;
  Tfactor = soil_con->Tfactor;
  AboveTreeLine = soil_con->AboveTreeLine;

  /* compute number of simulation days */
  Ndays = ( global_param.nrecs * global_param.dt) / 24;

  /* compute number of full model time steps per day */
  stepspday = 24/global_param.dt;
  
  if ( !param_set.TYPE[PREC].SUPPLIED
    && ( ( !param_set.TYPE[RAINF].SUPPLIED && ( !param_set.TYPE[LSRAINF].SUPPLIED || !param_set.TYPE[CRAINF].SUPPLIED ) )
      || ( ( !param_set.TYPE[SNOWF].SUPPLIED && ( !param_set.TYPE[LSSNOWF].SUPPLIED || !param_set.TYPE[CSNOWF].SUPPLIED ) ) ) ) )
    nrerror("Precipitation (PREC, or { {RAINF or {LSRAINF and CRAINF}} and {SNOWF or {LSSNOWF and CSNOWF}} }) must be given to the model, check input files\n");

  if ((!param_set.TYPE[TMAX].SUPPLIED || !param_set.TYPE[TMIN].SUPPLIED) && !param_set.TYPE[AIR_TEMP].SUPPLIED )
    nrerror("Daily maximum and minimum air temperature or sub-daily air temperature must be given to the model, check input files\n");
  
  if ( param_set.TYPE[AIR_TEMP].SUPPLIED && param_set.FORCE_DT[param_set.TYPE[AIR_TEMP].SUPPLIED-1] == 24 )
    nrerror("Model cannot use daily average temperature, must provide daily maximum and minimum or sub-daily temperatures.");
  
  if ( param_set.TYPE[SHORTWAVE].SUPPLIED
        && !(param_set.TYPE[VP].SUPPLIED || param_set.TYPE[QAIR].SUPPLIED || param_set.TYPE[REL_HUMID].SUPPLIED) )
    nrerror("Sub-daily shortwave and vapor pressure forcing data must be supplied together.");

  /* mtclim routine memory allocations */

  hourlyrad  = (double *) calloc(Ndays*24, sizeof(double));
  prec       = (double *) calloc(Ndays*24, sizeof(double));
  tair       = (double *) calloc(Ndays*24, sizeof(double));
  tmax       = (double *) calloc(Ndays, sizeof(double));
  tmaxhour   = (int *)    calloc(Ndays, sizeof(double));
  tmin       = (double *) calloc(Ndays, sizeof(double));
  tminhour   = (int *)    calloc(Ndays, sizeof(double));
  tskc       = (double *) calloc(Ndays*24, sizeof(double));
  vp         = (double *) calloc(Ndays*24, sizeof(double));
  
  if (hourlyrad == NULL || prec == NULL || tair == NULL || tmax == NULL ||
      tmaxhour == NULL ||  tmin == NULL || tminhour == NULL || tskc == NULL ||
      vp == NULL)
    nrerror("Memory allocation failure in initialize_atmos()");
  
  /*******************************
    read in meteorological data 
  *******************************/

  forcing_data = read_forcing_data(infile, global_param);
  
  fprintf(stderr,"\nRead meteorological forcing file\n");

  /*************************************************
    Pre-processing
  *************************************************/

  /*************************************************
    Convert units from ALMA to VIC standard, if necessary
  *************************************************/
  if (options.ALMA_INPUT) {
    for (type=0; type<N_FORCING_TYPES; type++) {
      if (param_set.TYPE[type].SUPPLIED) {
        /* Convert moisture flux rates to accumulated moisture flux per time step */
        if (   type == PREC
            || type == RAINF
            || type == CRAINF
            || type == LSRAINF
            || type == SNOWF
            || type == CSNOWF
            || type == LSSNOWF
            || type == CHANNEL_IN
           ) {
          for (idx=0; idx<(global_param.nrecs*NF); idx++) {
            forcing_data[type][idx] *= global_param.dt * 3600;
          }
        }
        /* Convert temperatures from K to C */
        else if (   type == AIR_TEMP
                 || type == TMIN
                 || type == TMAX
                ) {
          for (idx=0; idx<(global_param.nrecs*NF); idx++) {
            forcing_data[type][idx] -= KELVIN;
          }
        }
      }
    }
  }
  else {
    for (type=0; type<N_FORCING_TYPES; type++) {
      if (param_set.TYPE[type].SUPPLIED) {
        /* Convert pressures from kPa to Pa */
        if (   type == PRESSURE
            || type == VP
           ) {
          for (idx=0; idx<(global_param.nrecs*NF); idx++) {
            forcing_data[type][idx] *= kPa2Pa;
          }
        }
      }
    }
  }

  /*************************************************
    Incoming Upslope Runoff
  *************************************************/

  /*************************************************
    Create sub-daily channel_in if not provided
  *************************************************/

  if(param_set.TYPE[CHANNEL_IN].SUPPLIED) {
    if(param_set.FORCE_DT[param_set.TYPE[CHANNEL_IN].SUPPLIED-1] == 24) {
      /* daily channel_in provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
        for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < NF; j++) {
	    atmos[rec].channel_in[j] = forcing_data[CHANNEL_IN][day] 
	      / (float)(NF * stepspday);
	    atmos[rec].channel_in[j] *= 1000/cell_area; // convert to mm over grid cell 
	    sum += atmos[rec].channel_in[j];
	  }
	  if(NF>1) atmos[rec].channel_in[NR] = sum;
	  if(global_param.dt == 24) atmos[rec].channel_in[NR] = forcing_data[CHANNEL_IN][day];
	  rec++;
        }
      }
    }
    else {
      /* sub-daily channel_in provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        for(i = 0; i < NF; i++) {
	  atmos[rec].channel_in[i] = forcing_data[CHANNEL_IN][idx];
	  atmos[rec].channel_in[i] *= 1000/cell_area; // convert to mm over grid cell 
	  sum += atmos[rec].channel_in[i];
	  idx++;
        }
        if(NF>1) atmos[rec].channel_in[NR] = sum;
      }
    }
  }
  else {
    idx = 0;
    for(rec = 0; rec < global_param.nrecs; rec++) {
      sum = 0;
      for(i = 0; i < NF; i++) {
        atmos[rec].channel_in[i] = 0;
        sum += atmos[rec].channel_in[i];
        idx++;
      }
      if(NF>1) atmos[rec].channel_in[NR] = sum;
    }
  }

  /*************************************************
    Precipitation
  *************************************************/

  /*************************************************
    If provided, translate rainfall and snowfall
    into total precipitation
    NOTE: this overwrites any PREC data that was supplied
  *************************************************/

  if(param_set.TYPE[RAINF].SUPPLIED && param_set.TYPE[SNOWF].SUPPLIED) {
    /* rainfall and snowfall supplied */
    if (forcing_data[PREC] == NULL) {
      forcing_data[PREC] = (double *)calloc((global_param.nrecs * NF),sizeof(double));
    }
    for (idx=0; idx<(global_param.nrecs*NF); idx++) {
      forcing_data[PREC][idx] = forcing_data[RAINF][idx] + forcing_data[SNOWF][idx];
    }
    param_set.TYPE[PREC].SUPPLIED = param_set.TYPE[RAINF].SUPPLIED;
  }
  else if(param_set.TYPE[CRAINF].SUPPLIED && param_set.TYPE[LSRAINF].SUPPLIED
    && param_set.TYPE[CSNOWF].SUPPLIED && param_set.TYPE[LSSNOWF].SUPPLIED) {
    /* convective and large-scale rainfall and snowfall supplied */
    if (forcing_data[PREC] == NULL) {
      forcing_data[PREC] = (double *)calloc((global_param.nrecs * NF),sizeof(double));
    }
    for (idx=0; idx<(global_param.nrecs*NF); idx++) {
      forcing_data[PREC][idx] = forcing_data[CRAINF][idx] + forcing_data[LSRAINF][idx]
                               + forcing_data[CSNOWF][idx] + forcing_data[LSSNOWF][idx];
    }
    param_set.TYPE[PREC].SUPPLIED = param_set.TYPE[LSRAINF].SUPPLIED;
  }

  /*************************************************
    Create sub-daily precipitation if not provided
  *************************************************/

  if(param_set.FORCE_DT[param_set.TYPE[PREC].SUPPLIED-1] == 24) {
    /* daily prec provided */
    rec = 0;
    for (day = 0; day < Ndays; day++) {
      for (i = 0; i < stepspday; i++) {
	sum = 0;
	for (j = 0; j < NF; j++) {
	  atmos[rec].prec[j] = forcing_data[PREC][day] 
	    / (float)(NF * stepspday);
	  sum += atmos[rec].prec[j];
	}
	if(NF>1) atmos[rec].prec[NR] = sum;
	if(global_param.dt == 24) atmos[rec].prec[NR] = forcing_data[PREC][day];
	rec++;
      }
    }
  }
  else {
    /* sub-daily prec provided */
    idx = 0;
    for(rec = 0; rec < global_param.nrecs; rec++) {
      sum = 0;
      for(i = 0; i < NF; i++) {
	atmos[rec].prec[i] = forcing_data[PREC][idx];
	sum += atmos[rec].prec[i];
	idx++;
      }
      if(NF>1) atmos[rec].prec[NR] = sum;
    }
  }

  /*************************************************
    Air Temperature
  *************************************************/

  /************************************************
    Set maximum daily air temperature if provided 
  ************************************************/

  if(param_set.TYPE[TMAX].SUPPLIED) {
    if(param_set.FORCE_DT[param_set.TYPE[TMAX].SUPPLIED-1] == 24) {
      /* daily tmax provided */
      for (day = 0; day < Ndays; day++) {
	tmax[day] = forcing_data[TMAX][day];
      }
    }
    else {
      /* sub-daily tmax provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	tmax[rec/stepspday] = forcing_data[TMAX][idx];
	for(i = 0; i < NF; i++) idx++;
      }
    }
  }

  /************************************************
    Set minimum daily air temperature if provided 
  ************************************************/

  if(param_set.TYPE[TMIN].SUPPLIED) {
    if(param_set.FORCE_DT[param_set.TYPE[TMIN].SUPPLIED-1] == 24) {
      /* daily tmin provided */
      for (day = 0; day < Ndays; day++) {
	tmin[day] = forcing_data[TMIN][day];
      }
    }
    else {
      /* sub-daily tmin provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	tmin[rec/stepspday] = forcing_data[TMIN][idx];
	for(i = 0; i < NF; i++) idx++;
      }
    }
  }

  /*************************************************
    Store sub-daily air temperature if provided
  *************************************************/

  if(param_set.TYPE[AIR_TEMP].SUPPLIED) {
    /* forcing data defined as equal to or less than SNOW_STEP */
    idx = 0;
    for (rec = 0; rec < global_param.nrecs; rec++) {
      sum = 0;
      for (i = 0; i < NF; i++, step++) {
	atmos[rec].air_temp[i] = forcing_data[AIR_TEMP][idx];
	sum += atmos[rec].air_temp[i];
	idx++;
      }
      if(NF > 1) atmos[rec].air_temp[NR] = sum / (float)NF;
    }
  }

  /******************************************************
    Determine Tmax and Tmin from sub-daily temperatures
  ******************************************************/

  if(!(param_set.TYPE[TMAX].SUPPLIED && param_set.TYPE[TMIN].SUPPLIED)) {
    rec = 0;
    while(rec < global_param.nrecs) {
      min = max = atmos[rec].air_temp[0];
      for (j = 0; j < stepspday; j++) {
	for (i = 0; i < NF; i++, step++) {
	  if ( atmos[rec].air_temp[i] > max ) max = atmos[rec].air_temp[i];
	  if ( atmos[rec].air_temp[i] < min ) min = atmos[rec].air_temp[i];
	}
	rec++;
      }
      tmax[(rec-1)/stepspday] = max;
      tmin[(rec-1)/stepspday] = min;
    }
  }


  /*************************************************
    Pressures
  *************************************************/

  /*************************************************
    If provided, translate specific humidity and atm. pressure
    into vapor pressure
    NOTE: this overwrites any VP data that was supplied
  *************************************************/

  if(param_set.TYPE[QAIR].SUPPLIED && param_set.TYPE[PRESSURE].SUPPLIED) {
    /* specific humidity and atm. pressure supplied */
    if (forcing_data[VP] == NULL) {
      forcing_data[VP] = (double *)calloc((global_param.nrecs * NF),sizeof(double));
    }
    for (idx=0; idx<(global_param.nrecs*NF); idx++) {
      forcing_data[VP][idx] = forcing_data[QAIR][idx] * forcing_data[PRESSURE][idx] / EPS;
    }
    param_set.TYPE[VP].SUPPLIED = param_set.TYPE[QAIR].SUPPLIED;
  }

  /*************************************************
    If provided, translate relative humidity and atm. pressure
    into vapor pressure
    NOTE: this overwrites any VP data that was supplied
  *************************************************/

  if(param_set.TYPE[REL_HUMID].SUPPLIED && param_set.TYPE[PRESSURE].SUPPLIED) {
    /* relative humidity and atm. pressure supplied */
    if (forcing_data[VP] == NULL) {
      forcing_data[VP] = (double *)calloc((global_param.nrecs * NF),sizeof(double));
    }
    for (idx=0; idx<(global_param.nrecs*NF); idx++) {
      forcing_data[VP][idx] = forcing_data[REL_HUMID][idx] * svp(forcing_data[AIR_TEMP][idx]) / 100.;
    }
    param_set.TYPE[VP].SUPPLIED = param_set.TYPE[REL_HUMID].SUPPLIED;
  }


  /*************************************************
    Shortwave, vp, air temp, pressure, and density
  *************************************************/

  /**************************************************
    use the mtclim code to get the hourly shortwave 
    and the daily dew point temperature 

    requires prec, tmax, and tmin
  **************************************************/
  for (i = 0; i < Ndays; i++)
    prec[i] = 0;
  for (rec = 0; rec < global_param.nrecs; rec++) {
    prec[rec/stepspday] += atmos[rec].prec[NR];
  }
  mtclim42_wrapper(0, 0, (theta_l-theta_s)*24./360., elevation, annual_prec,
		   phi, &global_param, dmy, prec, tmax, tmin, tskc, vp,
		   hourlyrad);

  /***********************************************************
    reaggregate the hourly shortwave to the larger timesteps 
  ***********************************************************/

  if(!param_set.TYPE[SHORTWAVE].SUPPLIED) {
    for (rec = 0, hour = 0; rec < global_param.nrecs; rec++) {
      for (i = 0; i < NF; i++) {
	atmos[rec].shortwave[i] = 0;
	for (j = 0; j < options.SNOW_STEP; j++, hour++) {
	  atmos[rec].shortwave[i] += hourlyrad[hour];
	}
	atmos[rec].shortwave[i] /= options.SNOW_STEP;
      }
      if (NF > 1) {
	atmos[rec].shortwave[NR] = 0;
	for (i = 0; i < NF; i++) {
	  atmos[rec].shortwave[NR] += atmos[rec].shortwave[i];
	}
	atmos[rec].shortwave[NR] /= NF;
      }
    }
  }
  else {
    if(param_set.FORCE_DT[param_set.TYPE[SHORTWAVE].SUPPLIED-1] == 24) {
      /* daily shortwave provided; to get sub-daily, we will take the
         mtclim estimates and scale them to match the supplied daily totals */
      for (rec = 0, hour = 0; rec < global_param.nrecs; rec++) {
        for (i = 0; i < NF; i++) {
	  atmos[rec].shortwave[i] = 0;
	  for (j = 0; j < options.SNOW_STEP; j++, hour++) {
	    atmos[rec].shortwave[i] += hourlyrad[hour];
	  }
	  atmos[rec].shortwave[i] /= options.SNOW_STEP;
        }
        if (NF > 1) {
	  atmos[rec].shortwave[NR] = 0;
	  for (i = 0; i < NF; i++) {
	    atmos[rec].shortwave[NR] += atmos[rec].shortwave[i];
	  }
	  atmos[rec].shortwave[NR] /= NF;
        }
      }
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	if (forcing_data[SHORTWAVE][day] > 0 && atmos[rec].shortwave[NR] > 0)
	  factor = forcing_data[SHORTWAVE][day]/atmos[rec].shortwave[NR];
	else
	  factor = 0;
	for (i = 0; i < stepspday; i++) {
	  for (j = 0; j < NF; j++)
	    atmos[rec].shortwave[j] *= factor;
	  if(NF > 1)
	    atmos[rec].shortwave[NR] *= factor;
	  rec++;
        }
      }
    }
    else {
      /* sub-daily shortwave provided, so it will be used instead
         of the mtclim estimates */
      idx = 0;
      for (rec = 0, hour = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        for (i = 0; i < NF; i++, step++) {
	  atmos[rec].shortwave[i] = ( forcing_data[SHORTWAVE][idx] < 0 ) ? 0 : forcing_data[SHORTWAVE][idx];
	  sum += atmos[rec].shortwave[i];
	  idx++;
        }
        if (NF > 1) atmos[rec].shortwave[NR] = sum / (float)NF;
      }
    }
  }

  /**************************************************************************
    Calculate the hours at which the minimum and maximum temperatures occur
  **************************************************************************/
  if(!param_set.TYPE[AIR_TEMP].SUPPLIED) {
    set_max_min_hour(hourlyrad, Ndays, tmaxhour, tminhour);

  /**********************************************************************
    Calculate the subdaily and daily temperature based on tmax and tmin 
  **********************************************************************/

    HourlyT(options.SNOW_STEP, Ndays, tmaxhour, tmax, tminhour, tmin, tair);
    for (rec = 0, step = 0; rec < global_param.nrecs; rec++) {
      for (i = 0; i < NF; i++, step++) {
	atmos[rec].air_temp[i] = tair[step];
      }
      if (NF > 1) {
	atmos[rec].air_temp[NR] = 0;
	for (i = 0; i < NF; i++) {
	  atmos[rec].air_temp[NR] += atmos[rec].air_temp[i];
	}
	atmos[rec].air_temp[NR] /= NF;
      }
    }
  }

  /**************************************************
    calculate the subdaily and daily vapor pressure 
    and vapor pressure deficit
  **************************************************/

  if(!param_set.TYPE[VP].SUPPLIED) {
    for (rec = 0; rec < global_param.nrecs; rec++) {
      atmos[rec].vp[NR] = vp[rec/stepspday];
      atmos[rec].vpd[NR] = svp(atmos[rec].air_temp[NR]) - atmos[rec].vp[NR];

      if(atmos[rec].vpd[NR]<0) {
	atmos[rec].vpd[NR]=0;
	atmos[rec].vp[NR]=svp(atmos[rec].air_temp[NR]);
      }
      
      for (i = 0; i < NF; i++) {
	atmos[rec].vp[i]  = atmos[rec].vp[NR];
	atmos[rec].vpd[i]  = (svp(atmos[rec].air_temp[i]) - atmos[rec].vp[i]);

	if(atmos[rec].vpd[i]<0) {
	  atmos[rec].vpd[i]=0;
	  atmos[rec].vp[i]=svp(atmos[rec].air_temp[i]);
	}
      
      }
    }
  }
  else {
    if(param_set.FORCE_DT[param_set.TYPE[VP].SUPPLIED-1] == 24) {
      /* daily vp provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < NF; j++) {
	    atmos[rec].vp[j] = forcing_data[VP][day];
	    atmos[rec].vpd[j] = (svp(atmos[rec].air_temp[j]) 
				 - atmos[rec].vp[j]);
	    sum += atmos[rec].vp[j];
	  }
	  if(NF > 1) {
	    atmos[rec].vp[NR] = sum / (float)NF;
	    atmos[rec].vpd[NR] = (svp(atmos[rec].air_temp[NR]) 
				  - atmos[rec].vp[NR]);
	  }
	  rec++;
	}
      }
    }
    else {
      /* sub-daily vp provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < NF; i++) {
	  atmos[rec].vp[i] = forcing_data[VP][idx];
	  atmos[rec].vpd[i] = (svp(atmos[rec].air_temp[i]) 
			       - atmos[rec].vp[i]);
	  sum += atmos[rec].vp[i];
	  idx++;
	}
	if(NF > 1) {
	  atmos[rec].vp[NR] = sum / (float)NF;
	  atmos[rec].vpd[NR] = (svp(atmos[rec].air_temp[NR]) 
				    - atmos[rec].vp[NR]);
	}
      }
    }
  }

  /*************************************************
    Longwave
  *************************************************/

  /****************************************************************************
    calculate the daily and sub-daily longwave.  There is a separate case for
    the full energy and the water balance modes.  For water balance mode we 
    need to calculate the net longwave for the daily timestep and the incoming
    longwave for the SNOW_STEPs, for the full energy balance mode we always
    want the incoming longwave. 
  ****************************************************************************/

  if ( !param_set.TYPE[LONGWAVE].SUPPLIED ) {
    /** Incoming longwave radiation not supplied **/
    for (rec = 0; rec < global_param.nrecs; rec++) {
      sum = 0;
      for (i = 0; i < NF; i++) {
	calc_longwave(&(atmos[rec].longwave[i]), tskc[rec/stepspday],
		      atmos[rec].air_temp[i], atmos[rec].vp[i]);
        sum += atmos[rec].longwave[i];
      }
      if(NF>1) atmos[rec].longwave[NR] = sum / (float)NF;
    }
  }
  else if(param_set.FORCE_DT[param_set.TYPE[LONGWAVE].SUPPLIED-1] == 24) {
    /* daily incoming longwave radiation provided */
    rec = 0;
    for (day = 0; day < Ndays; day++) {
      for (i = 0; i < stepspday; i++) {
	sum = 0;
	for (j = 0; j < NF; j++) {
	  atmos[rec].longwave[j] = forcing_data[LONGWAVE][day];
	  sum += atmos[rec].longwave[j];
	}
	if(NF>1) atmos[rec].longwave[NR] = sum / (float)NF; 
	rec++;
      }
    }
  }
  else {
    /* sub-daily incoming longwave radiation provided */
    idx = 0;
    for(rec = 0; rec < global_param.nrecs; rec++) {
      sum = 0;
      for(i = 0; i < NF; i++) {
	atmos[rec].longwave[i] = forcing_data[LONGWAVE][idx];
	sum += atmos[rec].longwave[i];
	idx++;
      }
      if(NF>1) atmos[rec].longwave[NR] = sum / (float)NF; 
    }
  }

  /*************************************************
    Wind Speed
  *************************************************/

  /*************************************************
    If provided, translate WIND_E and WIND_N into WIND
    NOTE: this overwrites any WIND data that was supplied
  *************************************************/

  if(param_set.TYPE[WIND_E].SUPPLIED && param_set.TYPE[WIND_N].SUPPLIED) {
    /* specific humidity and atm. pressure supplied */
    if (forcing_data[WIND] == NULL) {
      forcing_data[WIND] = (double *)calloc((global_param.nrecs * NF),sizeof(double));
    }
    for (idx=0; idx<(global_param.nrecs*NF); idx++) {
      forcing_data[WIND][idx] = sqrt( forcing_data[WIND_E][idx]*forcing_data[WIND_E][idx]
                                    + forcing_data[WIND_N][idx]*forcing_data[WIND_N][idx] );
    }
    param_set.TYPE[WIND].SUPPLIED = param_set.TYPE[WIND_E].SUPPLIED;
  }

  /********************
    set the windspeed 
  ********************/

  if (param_set.TYPE[WIND].SUPPLIED) {
    if(param_set.FORCE_DT[param_set.TYPE[WIND].SUPPLIED-1] == 24) {
      /* daily wind provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < NF; j++) {
	    if(forcing_data[WIND][day] < options.MIN_WIND_SPEED)
	      atmos[rec].wind[j] = options.MIN_WIND_SPEED;
	    else 
	      atmos[rec].wind[j] = forcing_data[WIND][day];
	    sum += atmos[rec].wind[j];
	  }
	  if(NF>1) atmos[rec].wind[NR] = sum / (float)NF;
	  if(global_param.dt == 24) {
	    if(forcing_data[WIND][day] < options.MIN_WIND_SPEED)
	      atmos[rec].wind[j] = options.MIN_WIND_SPEED;
	    else 
	      atmos[rec].wind[NR] = forcing_data[WIND][day];
	  }
	  rec++;
	}
      }
    }
    else {
      /* sub-daily wind speed provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < NF; i++) {
	  if(forcing_data[WIND][idx] <options.MIN_WIND_SPEED)
	    atmos[rec].wind[i] = options.MIN_WIND_SPEED;
	  else
	    atmos[rec].wind[i] = forcing_data[WIND][idx];
	  sum += atmos[rec].wind[i];
	  idx++;
	}
	if(NF>1) atmos[rec].wind[NR] = sum / (float)NF;
      }
    }
  }
  else {
    /* no wind data provided, use default constant */
    for (rec = 0; rec < global_param.nrecs; rec++) {
      for (i = 0; i < NF; i++) {
	atmos[rec].wind[i] = 1.5;
      }
      atmos[rec].wind[NR] = 1.5;	
    }
  }

  /*************************************************
    Store atmospheric density if provided (kg/m^3)
  *************************************************/

  if(param_set.TYPE[DENSITY].SUPPLIED) {
    if(param_set.FORCE_DT[param_set.TYPE[DENSITY].SUPPLIED-1] == 24) {
      /* daily density provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < NF; j++) {
	    atmos[rec].density[j] = forcing_data[DENSITY][day];
	    sum += atmos[rec].density[j];
	  }
	  if(NF>1) atmos[rec].density[NR] = sum / (float)NF;
	  rec++;
	}
      }
    }
    else {
      /* sub-daily density provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < NF; i++) {
	  atmos[rec].density[i] = forcing_data[DENSITY][idx];
	  sum += atmos[rec].density[i];
	  idx++;
	}
	if(NF>1) atmos[rec].density[NR] = sum / (float)NF;
      }
    }
  }

  /**************************************
    Estimate Atmospheric Pressure (Pa) 
  **************************************/

  if(!param_set.TYPE[PRESSURE].SUPPLIED) {
    if(!param_set.TYPE[DENSITY].SUPPLIED) {
      /* Estimate pressure */
      if (options.PLAPSE) {
        /* Assume average virtual temperature in air column
           between ground and sea level = KELVIN+atmos[rec].air_temp[NR] + 0.5*elevation*LAPSE_PM */
        for (rec = 0; rec < global_param.nrecs; rec++) {
          atmos[rec].pressure[NR] = PS_PM*exp(-elevation*G/(Rd*(KELVIN+atmos[rec].air_temp[NR]+0.5*elevation*LAPSE_PM)));
          for (i = 0; i < NF; i++) {
            atmos[rec].pressure[i] = PS_PM*exp(-elevation*G/(Rd*(KELVIN+atmos[rec].air_temp[i]+0.5*elevation*LAPSE_PM)));
          }
        }
      }
      else {
        /* set pressure to constant value */
        for (rec = 0; rec < global_param.nrecs; rec++) {
	  atmos[rec].pressure[NR] = 95500.;
	  for (i = 0; i < NF; i++) {
	    atmos[rec].pressure[i] = atmos[rec].pressure[NR];
	  }
        }
      }
    }
    else {
      /* use observed densities to estimate pressure */
      if (options.PLAPSE) {
        for (rec = 0; rec < global_param.nrecs; rec++) {
          atmos[rec].pressure[NR] = (KELVIN+atmos[rec].air_temp[NR])*atmos[rec].density[NR]*Rd;
          for (i = 0; i < NF; i++) {
            atmos[rec].pressure[i] = (KELVIN+atmos[rec].air_temp[i])*atmos[rec].density[i]*Rd;
          }
        }
      }
      else {
        for (rec = 0; rec < global_param.nrecs; rec++) {
	  atmos[rec].pressure[NR] = (275.0 + atmos[rec].air_temp[NR])
	    *atmos[rec].density[NR]/0.003486;
	  for (i = 0; i < NF; i++) {
	    atmos[rec].pressure[i] = (275.0 + atmos[rec].air_temp[i])
	      *atmos[rec].density[i]/0.003486;
	  }
        }
      }
    }
  }
  else {
    /* observed atmospheric pressure supplied */
    if(param_set.FORCE_DT[param_set.TYPE[PRESSURE].SUPPLIED-1] == 24) {
      /* daily pressure provided */
      rec = 0;
      for (day = 0; day < Ndays; day++) {
	for (i = 0; i < stepspday; i++) {
	  sum = 0;
	  for (j = 0; j < NF; j++) {
	    atmos[rec].pressure[j] = forcing_data[PRESSURE][day];
	    sum += atmos[rec].pressure[j];
	  }
	  if(NF>1) atmos[rec].pressure[NR] = sum / (float)NF;
	  rec++;
	}
      }
    }
    else {
      /* sub-daily pressure provided */
      idx = 0;
      for(rec = 0; rec < global_param.nrecs; rec++) {
	sum = 0;
	for(i = 0; i < NF; i++) {
	  atmos[rec].pressure[i] = forcing_data[PRESSURE][idx];
	  sum += atmos[rec].pressure[i];
	  idx++;
	}
	if(NF>1) atmos[rec].pressure[NR] = sum / (float)NF;
      }
    }
  }

  /********************************************************
    Estimate Atmospheric Density if not provided (kg/m^3)
  ********************************************************/

  if(!param_set.TYPE[DENSITY].SUPPLIED) {
    /* use pressure to estimate density */
    if (options.PLAPSE) {
      for (rec = 0; rec < global_param.nrecs; rec++) {
        atmos[rec].density[NR] = atmos[rec].pressure[NR]/(Rd*(KELVIN+atmos[rec].air_temp[NR]));
        for (i = 0; i < NF; i++) {
          atmos[rec].density[i] = atmos[rec].pressure[i]/(Rd*(KELVIN+atmos[rec].air_temp[i]));
        }
      }
    }
    else {
      for (rec = 0; rec < global_param.nrecs; rec++) {
        atmos[rec].density[NR] = 0.003486*atmos[rec].pressure[NR]/
	  (275.0 + atmos[rec].air_temp[NR]);
        for (i = 0; i < NF; i++) {
	  atmos[rec].density[i] = 0.003486*atmos[rec].pressure[i]/
	    (275.0 + atmos[rec].air_temp[i]);
        }
      }
    }
  }

  /****************************************************
    Determine if Snow will Fall During Each Time Step
  ****************************************************/

#if !OUTPUT_FORCE
  min_Tfactor = Tfactor[0];
  for (band = 1; band < options.SNOW_BAND; band++) {
    if (Tfactor[band] < min_Tfactor)
      min_Tfactor = Tfactor[band];
  }
  for (rec = 0; rec < global_param.nrecs; rec++) {
    atmos[rec].snowflag[NR] = FALSE;
    for (i = 0; i < NF; i++) {
      if ((atmos[rec].air_temp[i] + min_Tfactor) < global_param.MAX_SNOW_TEMP
	  &&  atmos[rec].prec[i] > 0) {
	atmos[rec].snowflag[i] = TRUE;
	atmos[rec].snowflag[NR] = TRUE;
      }
      else
	atmos[rec].snowflag[i] = FALSE;
    }
  }
#endif // OUTPUT_FORCE
 
  // Free temporary parameters
  free(hourlyrad);
  free(prec);
  free(tair);
  free(tmax);
  free(tmaxhour);
  free(tmin);
  free(tminhour);
  free(tskc);
  free(vp);

  for(i=0;i<N_FORCING_TYPES;i++) 
    if (param_set.TYPE[i].SUPPLIED) 
      free((char *)forcing_data[i]);
  free((char *)forcing_data);

#if OUTPUT_FORCE_STATS
  calc_forcing_stats(global_param.nrecs, atmos);
#endif // OUTPUT_FORCE_STATS

#if !OUTPUT_FORCE

  // If COMPUTE_TREELINE is TRUE and the treeline computation hasn't
  // specifically been turned off for this cell (by supplying avgJulyAirTemp
  // and setting it to -999), calculate which snowbands are above the
  // treeline, based on average July air temperature.
  if (options.COMPUTE_TREELINE) {
    if ( !(options.JULY_TAVG_SUPPLIED && avgJulyAirTemp == -999) ) {
      compute_treeline( atmos, dmy, avgJulyAirTemp, Tfactor, AboveTreeLine );
    }
  }

#else

  // If OUTPUT_FORCE is set to TRUE in user_def.h then the full
  // forcing data array is dumped into a new set of files.
  write_forcing_file(atmos, global_param.nrecs, out_data_files, out_data);

#endif // OUTPUT_FORCE 

}
