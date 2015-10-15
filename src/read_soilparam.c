#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

soil_con_struct read_soilparam(FILE *soilparam,
			       char *RUN_MODEL,
			       char *MODEL_DONE)
/**********************************************************************
	read_soilparam		Dag Lohmann		January 1996

  This routine reads soil parameters for each grid cell.

  Parameters Read from File:
  TYPE   NAME                    UNITS   DESCRIPTION
  int    gridcel;                N/A     grid cell number
  float  lat;		         degrees grid cell central latitude
  float  lng;		         degrees grid cell central longitude
  double b_infilt;  	         N/A     infiltration parameter
  double Ds;		         fract   fraction of maximum subsurface
                                         flow rate
  double Dsmax;  	         mm/day  maximum subsurface flow rate
  double Ws;		         fract   fraction of maximum soil moisture
  double c;                      N/A     exponent in ARNO baseflow curve
  double expt[MAX_LAYERS];         N/A     exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6
  double Ksat[MAX_LAYERS];         mm/day  saturated hydraulic conductivity
  double phi_s[MAX_LAYERS];        mm/mm   saturated matrix potential
  double init_moist[MAX_LAYERS];   mm      initial layer moisture level
  float  elevation;	         m       grid cell elevation
  double depth[MAX_LAYERS];        m       thickness of each layer
  double avg_temp;	         C       average soil temperature
  double dp;		         m       soil thermal damping depth
  double bubble;	         cm      bubbling pressure, HBH 5.15
  double quartz;	         fract   quartz content of soil
  double organic;	         fract   organic content of soil
  double bulk_density[MAX_LAYERS]; kg/m^3  soil bulk density
  double soil_density;		 kg/m^3  soil partical density
  double rough;		         m       soil surface roughness
  double snow_rough;             m       snow surface roughness

  Parameters Computed from those in the File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double max_moist[MAX_LAYERS];    mm      maximum moisture content per layer
  double max_infil;	         N/A     maximum infiltration rate
  double Wcr[MAX_LAYERS];	         mm      critical moisture level for soil
                                         layer, evaporation is no longer
                                         affected moisture stress in the soil
  double Wpwp[MAX_LAYERS];         mm      soil moisture content at permanent
                                         wilting point
  float  time_zone_lng;	         degrees central meridian of the time zone


  Modifications:
  7-19-96	Modified to read through variable layers, and
		read soil depth and average temperature for
		full energy and frozen soil versions of the
		model.						KAC
  4-12-98       Modified to read all parameters from a single
                standard input file.                            KAC
  3-13-00       Modified to read more parameters as separate
                layer values                                    KAC
  6-6-2000      Modified to skip individual parameter reads
                if model grid cell is not read.                 KAC
  xx-xx-01      Modified to read in spatial snow and frost
                parameters.                                     KAC
  11-18-02      Modified to read Bart's new Arno parameters.    IHA
  10-May-04     Replaced rint(something) with (float)(int)(something + 0.5)
		to handle rounding without resorting to rint().		TJB
  11-May-04	(fix by Chunmei Zhu and Alan Hamlet)
		Added check to make sure that wilting point is
		greater than residual moisture.				TJB
  07-Jul-04	Changed lower limit on initial soil moisture to be
		residual moisture instead of wilting point.  Also
		cleaned up validation statements.			TJB
  07-Jul-04	Removed extraneous tmp variable.			TJB
  07-Jul-04	Only validate initial soil moisture if INIT_STATE
		is FALSE.						TJB
  26-Oct-04	Added validation of depth_full_snow_cover and
		frost_slope.						TJB
  2005-Apr-13 Added logic for OUTPUT_FORCE option.			TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW.		TJB
  2006-Sep-13 Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option.	TJB/GCT
  2007-May-23 Replaced 'fscanf' statements with 'sscanf' statements
	      to trap missing fields.					GCT
  2007-Aug-08 Added EXCESS_ICE option.					JCA
  2007-Sep-14 Clarified description in comment before BASEFLOW check.	TJB
  2007-Nov-06 Moved computation of cell_area from read_lakeparam() to
	      here.							TJB
  2009-Jan-12 Added logic for JULY_TAVG_SUPPLIED.			TJB
  2009-May-22 Added validation of expt and bubble.			TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-17 Modified to understand both tabs and spaces as delimiters.TJB
  2009-Jul-31 Removed unused layer_node_fract array.			TJB
  2009-Sep-11 Added correct OUTPUT_FORCE logic around the new bare
	      soil/veg lib code.					TJB
  2009-Sep-28 Added initialization of snowband parameters.		TJB
  2011-Jan-04 Made read_soilparam_arc() a sub-function of
	      read_soilparam().						TJB
  2011-Jan-04 Added computation of relationship between soil moisture
	      and water table depth given by soil water retention curve.TJB
  2011-Mar-01 Updated soil water retention curve relationship.		TJB
  2011-Mar-05 Now does fgets whether or not RUN_MODEL is true - this
	      fixes bug introduced when read_soilparam_arc() was moved
	      to a sub-function of read_soilparam().			TJB
  2011-May-25 Expanded latchar, lngchar, and junk allocations to handle
	      GRID_DECIMAL > 4.						TJB
  2011-May-25 Moved fgets() so that it always gets called, even when
	      cells are skipped.					TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  If TRUE, VIC expects organic
	      fraction and organic bulk and soil densities to be supplied
	      for each grid cell.					TJB
  2011-Sep-28 Added validation of b_infilt.				TJB
  2011-Nov-04 Added hard-coding of slope, aspect, and horizons to 0.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
  2013-Jul-25 Added calculation of soil albedo in PAR range.		TJB
  2013-Dec-26 Removed EXCESS_ICE option.							TJB
  2013-Dec-27 Moved SPATIAL_SNOW from compile-time to run-time options.	TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct.			TJB
  2014-Mar-24 Removed ARC_SOIL option                               BN
  2014-Mar-28 Removed DIST_PRCP option.								TJB
**********************************************************************/
{
  void ttrim( char *string );
  extern option_struct options;
  extern global_param_struct global_param;
  extern veg_lib_struct *veg_lib;
  char            ErrStr[MAXSTRING];
  char            line[MAXSTRING];
  char            tmpline[MAXSTRING];
  const char      delimiters[] = " \t";
  char            *token;
  int             layer, i, tempint, j;
  double          Wcr_FRACT[MAX_LAYERS];
  double          Wpwp_FRACT[MAX_LAYERS];
  double          off_gmt;
  double          tempdbl;
  double          extra_depth;
  double          lat;
  double          lng;
  double          start_lat;
  double          right_lng;
  double          left_lng;
  double          delta;
  double          dist;
  size_t          length;
  int             Nbands,band;
  int             Ncells;
  int             flag;
  double          tmp_depth;
  double          tmp_depth2, tmp_depth2_save;
  double          b, b_save;
  double          bubble, bub_save;
  double          tmp_max_moist;
  double          tmp_resid_moist;
  double          zwt_prime, zwt_prime_eff;
  double          tmp_moist;
  double          w_avg;
  char   latchar[20], lngchar[20], junk[6];
  soil_con_struct temp;

    /** Read plain ASCII soil parameter file **/
  if ((fscanf(soilparam, "%d", &flag)) != EOF) {
    if (flag) {
      *RUN_MODEL = TRUE;
    }
    else {
      *RUN_MODEL = FALSE;
    }

    if (fgets( line, MAXSTRING, soilparam ) == NULL) {
      sprintf(ErrStr,"ERROR: Unexpected EOF while reading soil file\n");
      nrerror(ErrStr);
    }
  }
  else {
    *MODEL_DONE = TRUE;
    *RUN_MODEL = FALSE;
  }

    if(!(*MODEL_DONE) && (*RUN_MODEL)) {

      strcpy(tmpline, line);
      ttrim( tmpline );
      if( ( token = strtok (tmpline, delimiters)) == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for CELL NUMBER in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%d", &temp.gridcel);
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for CELL LATITUDE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%f", &temp.lat);
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for CELL LONGITUDE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%f", &temp.lng);
#if VERBOSE
      /* add print statements for grid cell number -- EDM */
      fprintf(stderr,"\ncell: %d,  lat: %.4f, long: %.4f\n",temp.gridcel,temp.lat,temp.lng);
#endif

      /* read infiltration parameter */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for INFILTRATION in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.b_infilt);
      if (temp.b_infilt <= 0) {
        sprintf(ErrStr,"ERROR: b_infilt (%f) in soil file is <= 0; b_infilt must be positive\n",temp.b_infilt);
        nrerror(ErrStr);
      }

      /* read fraction of baseflow rate */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for FRACTION OF BASEFLOW RATE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.Ds);

      /* read maximum baseflow rate */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for MAXIMUM BASEFLOW RATE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.Dsmax);

      /* read fraction of bottom soil layer moisture */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for FRACTION OF BOTTOM SOIL LAYER MOISTURE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.Ws);

      /* read exponential */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for EXPONENTIAL in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.c);

      /* read expt for each layer */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for EXPT for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.expt[layer]);
        if (!options.OUTPUT_FORCE) {
          if(temp.expt[layer] < 3.0) {
            fprintf(stderr,"ERROR: Exponent in layer %d is %f < 3.0; This must be > 3.0\n", layer, temp.expt[layer]);
            exit(0);
          }
        }
      }

      /* read layer saturated hydraulic conductivity */
      for(layer = 0; layer < options.Nlayer; layer++){
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for SATURATED HYDRAULIC CONDUCTIVITY for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.Ksat[layer]);
      }

      /* read layer phi_s */
      for(layer = 0; layer < options.Nlayer; layer++){
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for PHI_S for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.phi_s[layer]);
      }

      /* read layer initial moisture */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for INITIAL MOISTURE for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.init_moist[layer]);
        if (!options.OUTPUT_FORCE) {
          if(temp.init_moist[layer] < 0.) {
            sprintf(ErrStr,"ERROR: Initial moisture for layer %d cannot be negative (%f)",layer,temp.init_moist[layer]);
            nrerror(ErrStr);
          }
        }
      }

      /* read cell mean elevation */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for CELL MEAN ELEVATION in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%f", &temp.elevation);

      /* soil layer thicknesses */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for LAYER THICKNESS for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.depth[layer]);
      }
      if (!options.OUTPUT_FORCE) {
        /* round soil layer thicknesses to nearest mm */
        for(layer = 0; layer < options.Nlayer; layer++)
          temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
      }

      /* read average soil temperature */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for AVERAGE SOIL TEMPERATURE in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.avg_temp);
      if (!options.OUTPUT_FORCE) {
        if(options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
          fprintf(stderr,"Need valid average soil temperature in degrees C to run");
          fprintf(stderr," Full Energy model, %f is not acceptable.\n",
            temp.avg_temp);
          exit(0);
        }
      }

      /* read soil damping depth */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for SOIL DAMPING DEPTH in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.dp);

      /* read layer bubbling pressure */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for BUBBLING PRESSURE for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.bubble[layer]);
        if (!options.OUTPUT_FORCE) {
          if((options.FULL_ENERGY || options.FROZEN_SOIL) && temp.bubble[layer] < 0) {
            fprintf(stderr,"ERROR: Bubbling pressure in layer %d is %f < 0; This must be positive for FULL_ENERGY = TRUE or FROZEN_SOIL = TRUE\n", layer, temp.bubble[layer]);
            exit(0);
          }
        }
      }

      /* read layer quartz content */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for QUARTZ CONTENT for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.quartz[layer]);
        if (!options.OUTPUT_FORCE) {
          if(options.FULL_ENERGY && (temp.quartz[layer] > 1. || temp.quartz[layer] < 0)) {
            fprintf(stderr,"Need valid quartz content as a fraction to run");
            fprintf(stderr," Full Energy model, %f is not acceptable.\n", temp.quartz[layer]);
            exit(0);
          }
        }
      }

      /* read layer bulk density */
      for(layer = 0; layer < options.Nlayer; layer++){
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for mineral BULK DENSITY for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.bulk_dens_min[layer]);
        if (!options.OUTPUT_FORCE) {
          if(temp.bulk_dens_min[layer] <= 0) {
            sprintf(ErrStr,"ERROR: layer %d mineral bulk density (%f) must be > 0", layer, temp.bulk_dens_min[layer] );
            nrerror(ErrStr);
          }
        }
      }

      /* read layer soil density */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for mineral SOIL DENSITY for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.soil_dens_min[layer]);
        if (!options.OUTPUT_FORCE) {
          if(temp.soil_dens_min[layer] <= 0) {
            sprintf(ErrStr,"ERROR: layer %d mineral soil density (%f) must be > 0", layer, temp.soil_dens_min[layer] );
            nrerror(ErrStr);
          }
          if(temp.bulk_dens_min[layer]>=temp.soil_dens_min[layer]) {
            sprintf(ErrStr,"ERROR: layer %d mineral bulk density (%f) must be less than mineral soil density (%f)", layer, temp.bulk_dens_min[layer], temp.soil_dens_min[layer] );
            nrerror(ErrStr);
          }
        }
      }

      if (options.ORGANIC_FRACT) {
        /* read layer organic content */
        for(layer = 0; layer < options.Nlayer; layer++){
          token = strtok (NULL, delimiters);
          while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
          if( token == NULL ) {
            sprintf(ErrStr,"ERROR: Can't find values for ORGANIC CONTENT for layer %d in soil file\n", layer );
            nrerror(ErrStr);
          }
          sscanf(token, "%lf", &temp.organic[layer]);
          if (!options.OUTPUT_FORCE) {
            if(temp.organic[layer] > 1. || temp.organic[layer] < 0) {
              sprintf(ErrStr,"ERROR: Need valid volumetric organic soil fraction when options.ORGANIC_FRACT is set to TRUE.\n  %f is not acceptable.\n", temp.organic[layer]);
              nrerror(ErrStr);
            }
          }
        }

        /* read layer bulk density */
        for(layer = 0; layer < options.Nlayer; layer++){
          token = strtok (NULL, delimiters);
          while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
          if( token == NULL ) {
            sprintf(ErrStr,"ERROR: Can't find values for organic BULK DENSITY for layer %d in soil file\n", layer );
            nrerror(ErrStr);
          }
          sscanf(token, "%lf", &temp.bulk_dens_org[layer]);
          if (!options.OUTPUT_FORCE) {
            if(temp.bulk_dens_org[layer] <= 0 && temp.organic[layer] > 0) {
              fprintf(stderr,"WARNING: layer %d organic bulk density (%f) must be > 0; setting to mineral bulk density (%f)\n", layer, temp.bulk_dens_org[layer], temp.bulk_dens_min[layer] );
              temp.bulk_dens_org[layer] = temp.bulk_dens_min[layer];
            }
          }
        }

        /* read layer soil density */
        for(layer = 0; layer < options.Nlayer; layer++) {
          token = strtok (NULL, delimiters);
          while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
          if( token == NULL ) {
            sprintf(ErrStr,"ERROR: Can't find values for organic SOIL DENSITY for layer %d in soil file\n", layer );
            nrerror(ErrStr);
          }
          sscanf(token, "%lf", &temp.soil_dens_org[layer]);
          if (!options.OUTPUT_FORCE) {
            if(temp.soil_dens_org[layer] <= 0 && temp.organic[layer] > 0) {
              fprintf(stderr,"WARNING: layer %d organic soil density (%f) must be > 0; setting to mineral soil density (%f)\n", layer, temp.soil_dens_org[layer], temp.soil_dens_min[layer] );
              temp.soil_dens_org[layer] = temp.soil_dens_min[layer];
            }
            if(temp.organic[layer] > 0 && temp.bulk_dens_org[layer]>=temp.soil_dens_org[layer]) {
              sprintf(ErrStr,"ERROR: layer %d organic bulk density (%f) must be less than organic soil density (%f)", layer, temp.bulk_dens_org[layer], temp.soil_dens_org[layer] );
              nrerror(ErrStr);
            }
          }
        }

      }
      else {
        for(layer = 0; layer < options.Nlayer; layer++){
          temp.organic[layer] = 0.0;
          temp.bulk_dens_org[layer] = -9999;
          temp.soil_dens_org[layer] = -9999;
        }
      }

      /* read cell gmt offset */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for GMT OFFSET in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &off_gmt);

      /* read layer critical point */
      for(layer=0;layer<options.Nlayer;layer++){
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for CRITICAL POINT for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &(Wcr_FRACT[layer]));
      }

      /* read layer wilting point */
      for(layer=0;layer<options.Nlayer;layer++){
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for WILTING POINT for layer %d in soil file\n", layer );
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &(Wpwp_FRACT[layer]));
      }

      /* read soil roughness */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for SOIL ROUGHNESS in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.rough);
      if (!options.OUTPUT_FORCE) {
        /* Overwrite default bare soil aerodynamic resistance parameters
           with the values taken from the soil parameter file */
        for (j=0; j<12; j++) {
          veg_lib[veg_lib[0].NVegLibTypes].roughness[j] = temp.rough;
          veg_lib[veg_lib[0].NVegLibTypes].displacement[j] = temp.rough*0.667/0.123;
        }
      }

      /* read snow roughness */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for SNOW ROUGHNESS in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.snow_rough);

      /* read cell annual precipitation */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for ANNUAL PRECIPITATION in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%lf", &temp.annual_prec);

      /* read layer residual moisture content */
      for(layer = 0; layer < options.Nlayer; layer++) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for RESIDUAL MOISTURE CONTENT for layer %d in soil file\n", layer);
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &temp.resid_moist[layer]);
      }

      /* read frozen soil active flag */
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for FROZEN SOIL ACTIVE FLAG in soil file\n");
        nrerror(ErrStr);
      }
      sscanf(token, "%d", &tempint);
      temp.FS_ACTIVE = (char)tempint;

      /* read minimum snow depth for full coverage */
      if (options.SPATIAL_SNOW) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for SPATIAL SNOW in soil file\n");
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &tempdbl);
        temp.max_snow_distrib_slope = tempdbl;
      }
      else
        temp.max_snow_distrib_slope = 0;

      /* read slope of frozen soil distribution */
      if (options.SPATIAL_FROST) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for SPATIAL FROST in soil file\n");
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &tempdbl);
        temp.frost_slope = tempdbl;
      }
      else
        temp.frost_slope = 0;

      /* If specified, read cell average July air temperature in the final
         column of the soil parameter file */
      if (options.JULY_TAVG_SUPPLIED) {
        token = strtok (NULL, delimiters);
        while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
        if( token == NULL ) {
          sprintf(ErrStr,"ERROR: Can't find values for average July Tair in soil file\n");
          nrerror(ErrStr);
        }
        sscanf(token, "%lf", &tempdbl);
        temp.avgJulyAirTemp = tempdbl;
      }

      /*******************************************
        End of soil parameters for this grid cell
      *******************************************/

      if (!options.OUTPUT_FORCE) {

        /*******************************************
          Compute Soil Layer Properties
        *******************************************/
        for(layer = 0; layer < options.Nlayer; layer++) {
          temp.bulk_density[layer] = (1-temp.organic[layer])*temp.bulk_dens_min[layer] + temp.organic[layer]*temp.bulk_dens_org[layer];
          temp.soil_density[layer] = (1-temp.organic[layer])*temp.soil_dens_min[layer] + temp.organic[layer]*temp.soil_dens_org[layer];
          if (temp.resid_moist[layer] == MISSING)
              temp.resid_moist[layer] = RESID_MOIST;
          temp.porosity[layer] = 1.0 - temp.bulk_density[layer] / temp.soil_density[layer];
          temp.max_moist[layer] = temp.depth[layer] * temp.porosity[layer] * 1000.;
        }

        /**********************************************
          Validate Soil Layer Thicknesses
        **********************************************/
        for(layer = 0; layer < options.Nlayer; layer++) {
          if(temp.depth[layer] < MINSOILDEPTH) {
            sprintf(ErrStr,"ERROR: Model will not function with layer %d depth %f < %f m.\n",
            layer,temp.depth[layer],MINSOILDEPTH);
            nrerror(ErrStr);
          }
        }
        if(temp.depth[0] > temp.depth[1]) {
          sprintf(ErrStr,"ERROR: Model will not function with layer %d depth (%f m) > layer %d depth (%f m).\n",
            0,temp.depth[0],1,temp.depth[1]);
          nrerror(ErrStr);
        }

        /**********************************************
          Compute Maximum Infiltration for Upper Layers
        **********************************************/
        if(options.Nlayer==2)
          temp.max_infil = (1.0+temp.b_infilt)*temp.max_moist[0];
        else
          temp.max_infil = (1.0+temp.b_infilt)*(temp.max_moist[0]+temp.max_moist[1]);

        /****************************************************************
          Compute Soil Layer Critical and Wilting Point Moisture Contents
        ****************************************************************/
        for(layer=0;layer<options.Nlayer;layer++) {
          temp.Wcr[layer]  = Wcr_FRACT[layer] * temp.max_moist[layer];
          temp.Wpwp[layer] = Wpwp_FRACT[layer] * temp.max_moist[layer];
          if(temp.Wpwp[layer] > temp.Wcr[layer]) {
            sprintf(ErrStr,"Calculated wilting point moisture (%f mm) is greater than calculated critical point moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
            temp.Wpwp[layer], temp.Wcr[layer], layer);
            nrerror(ErrStr);
          }
          if(temp.Wpwp[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) {
            sprintf(ErrStr,"Calculated wilting point moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
            temp.Wpwp[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
            nrerror(ErrStr);
          }
	  temp.Wcr_orig[layer]=temp.Wcr[layer]; //ingjerd jan 2015
	  temp.Wpwp_orig[layer]=temp.Wpwp[layer]; //ingjerd jan 2015
	  temp.Wcr_irrig[layer]=0.557*temp.max_moist[layer]; //ingjerd jan 2015
	  temp.Wpwp_irrig[layer]=0.436*temp.max_moist[layer]; //ingjerd jan 2015
        }

        /**********************************************
          Validate Spatial Snow/Frost Params
        **********************************************/
        if (options.SPATIAL_SNOW) {
          if (temp.max_snow_distrib_slope < 0.0) {
            sprintf(ErrStr,"max_snow_distrib_slope (%f) must be positive.\n", temp.max_snow_distrib_slope);
            nrerror(ErrStr);
          }
        }

        if (options.SPATIAL_FROST) {
          if (temp.frost_slope < 0.0) {
            sprintf(ErrStr,"frost_slope (%f) must be positive.\n", temp.frost_slope);
            nrerror(ErrStr);
          }
        }

        /*************************************************
          If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
          parameters d1, d2, d3, and d4 to ARNO baseflow
          parameters Ds, Dsmax, Ws, and c
        *************************************************/
        if(options.BASEFLOW == NIJSSEN2001) {
          layer = options.Nlayer-1;
          temp.Dsmax = temp.Dsmax *
            pow((double)(1./(temp.max_moist[layer]-temp.Ws)), -temp.c) +
            temp.Ds * temp.max_moist[layer];
          temp.Ds = temp.Ds * temp.Ws / temp.Dsmax;
          temp.Ws = temp.Ws/temp.max_moist[layer];
        }

        /*******************************************************************
          Calculate grid cell area.
        ******************************************************************/

        if (options.EQUAL_AREA) {

          temp.cell_area = global_param.resolution * 1000. * 1000.; /* Grid cell area in m^2. */

        }
        else {

          lat = fabs(temp.lat);
          lng = fabs(temp.lng);

          start_lat = lat - global_param.resolution / 2;
          right_lng = lng + global_param.resolution / 2;
          left_lng  = lng - global_param.resolution / 2;

          delta = get_dist(lat,lng,lat+global_param.resolution/10.,lng);

          dist = 0.;

          for ( i = 0; i < 10; i++ ) {
            dist += get_dist(start_lat,left_lng,start_lat,right_lng) * delta;
            start_lat += global_param.resolution/10;
          }

          temp.cell_area = dist * 1000. * 1000.; /* Grid cell area in m^2. */

        }

        /*************************************************
          Allocate and Initialize Snow Band Parameters
        *************************************************/
        Nbands         = options.SNOW_BAND;
        temp.AreaFract     = (double *)calloc(Nbands,sizeof(double));
        temp.BandElev      = (float *)calloc(Nbands,sizeof(float));
        temp.Tfactor       = (double *)calloc(Nbands,sizeof(double));
        temp.Pfactor       = (double *)calloc(Nbands,sizeof(double));
        temp.AboveTreeLine = (char *)calloc(Nbands,sizeof(char));

        if (temp.Tfactor == NULL || temp.Pfactor == NULL || temp.AreaFract == NULL)
          nrerror("Memory allocation failure in read_snowband");

        if ( Nbands <= 0 ) {
          sprintf(ErrStr,"Number of snow bands must be > 0 (%d)",Nbands);
          nrerror(ErrStr);
        }

        /** Set default values for factors to use unmodified forcing data **/
        for (band = 0; band < Nbands; band++) {
          temp.AreaFract[band] = 0.;
          temp.BandElev[band]  = temp.elevation;
          temp.Tfactor[band]   = 0.;
          temp.Pfactor[band]   = 1.;
        }
        temp.AreaFract[0] = 1.;

        /*************************************************
          Compute soil moistures for various values of water table depth
          Here we use the relationship (e.g., Letts et al., 2000)
            w(z) = { ((zwt-z)/bubble)**(-1/b), z <  zwt-bubble
                   { 1.0,                      z >= zwt-bubble
          where
            z      = depth below surface [cm]
            w(z)   = relative moisture at depth z given by
                     (moist(z) - resid_moist) / (max_moist - resid_moist)
            zwt    = depth of water table below surface [cm]
            bubble = bubbling pressure [cm]
            b      = 0.5*(expt-3)
          Note that zwt-bubble = depth of the free water surface, i.e.
          position below which soil is completely saturated.

          This assumes water in unsaturated zone above water table
          is always in equilibrium between gravitational and matric
          tension (e.g., Frolking et al, 2002).

          So, to find the soil moisture value in a layer corresponding
          to a given water table depth zwt, we integrate w(z) over the
          whole layer:

          w_avg = average w over whole layer = (integral of w*dz) / layer depth

          Then,
            layer moisture = w_avg * (max_moist - resid_moist) + resid_moist

          Instead of the zwt defined above, will actually report free
          water surface elevation zwt' = -(zwt-bubble).  I.e. zwt' < 0
          below the soil surface, and marks the point of saturation
          rather than pressure = 1 atm.

          Do this for each layer individually and also for a) the top N-1 layers
          lumped together, and b) the entire soil column lumped together.

        *************************************************/

        /* Individual layers */
        tmp_depth = 0;
        for (layer=0; layer<options.Nlayer; layer++) {
          b = 0.5*(temp.expt[layer]-3);
          bubble = temp.bubble[layer];
          tmp_resid_moist = temp.resid_moist[layer]*temp.depth[layer]*1000; // in mm
          zwt_prime = 0; // depth of free water surface below top of layer (not yet elevation)
          for (i=0; i<MAX_ZWTVMOIST; i++) {
            temp.zwtvmoist_zwt[layer][i] = -tmp_depth*100-zwt_prime; // elevation (cm) relative to soil surface
            w_avg = ( temp.depth[layer]*100 - zwt_prime
                     - (b/(b-1))*bubble*(1-pow((zwt_prime+bubble)/bubble,(b-1)/b)) )
                    / (temp.depth[layer]*100); // in cm
            if (w_avg < 0) w_avg = 0;
            if (w_avg > 1) w_avg = 1;
            temp.zwtvmoist_moist[layer][i] = w_avg*(temp.max_moist[layer]-tmp_resid_moist)+tmp_resid_moist;
            zwt_prime += temp.depth[layer]*100/(MAX_ZWTVMOIST-1); // in cm
          }
          tmp_depth += temp.depth[layer];
        }

        /* Top N-1 layers lumped together (with average soil properties) */
        tmp_depth = 0;
        b = 0;
        bubble = 0;
        tmp_max_moist = 0;
        tmp_resid_moist = 0;
        for (layer=0; layer<options.Nlayer-1; layer++) {
          b += 0.5*(temp.expt[layer]-3)*temp.depth[layer];
          bubble += temp.bubble[layer]*temp.depth[layer];
          tmp_max_moist += temp.max_moist[layer]; // total max_moist
          tmp_resid_moist += temp.resid_moist[layer]*temp.depth[layer]*1000; // total resid_moist in mm
          tmp_depth += temp.depth[layer];
        }
        b /= tmp_depth; // average b
        bubble /= tmp_depth; // average bubble
        zwt_prime = 0; // depth of free water surface below top of layer (not yet elevation)
        for (i=0; i<MAX_ZWTVMOIST; i++) {
          temp.zwtvmoist_zwt[options.Nlayer][i] = -zwt_prime; // elevation (cm) relative to soil surface
          w_avg = ( tmp_depth*100 - zwt_prime
                     - (b/(b-1))*bubble*(1-pow((zwt_prime+bubble)/bubble,(b-1)/b)) )
                    / (tmp_depth*100); // in cm
          if (w_avg < 0) w_avg = 0;
          if (w_avg > 1) w_avg = 1;
          temp.zwtvmoist_moist[options.Nlayer][i] = w_avg*(tmp_max_moist-tmp_resid_moist)+tmp_resid_moist;
          zwt_prime += tmp_depth*100/(MAX_ZWTVMOIST-1); // in cm
        }

        /* Compute zwt by taking total column soil moisture and filling column from bottom up */
        tmp_depth = 0;
        for (layer=0; layer<options.Nlayer; layer++) {
          tmp_depth += temp.depth[layer];
        }
        zwt_prime = 0; // depth of free water surface below soil surface (not yet elevation)
        for (i=0; i<MAX_ZWTVMOIST; i++) {
          temp.zwtvmoist_zwt[options.Nlayer+1][i] = -zwt_prime; // elevation (cm) relative to soil surface
          // Integrate w_avg in pieces
          if (zwt_prime == 0) {
            tmp_moist = 0;
            for (layer=0; layer<options.Nlayer; layer++)
              tmp_moist += temp.max_moist[layer];
            temp.zwtvmoist_moist[options.Nlayer+1][i] = tmp_moist;
          }
          else {
            tmp_moist = 0;
            layer = options.Nlayer-1;
            tmp_depth2 = tmp_depth-temp.depth[layer];
            while (layer>0 && zwt_prime <= tmp_depth2*100) {
              tmp_moist += temp.max_moist[layer];
              layer--;
              tmp_depth2 -= temp.depth[layer];
            }
            w_avg = (tmp_depth2*100+temp.depth[layer]*100-zwt_prime)/(temp.depth[layer]*100);
            b = 0.5*(temp.expt[layer]-3);
            bubble = temp.bubble[layer];
            tmp_resid_moist = temp.resid_moist[layer]*temp.depth[layer]*1000;
            w_avg += -(b/(b-1))*bubble*( 1 - pow((zwt_prime+bubble-tmp_depth2*100)/bubble,(b-1)/b) ) / (temp.depth[layer]*100);
            tmp_moist += w_avg*(temp.max_moist[layer]-tmp_resid_moist)+tmp_resid_moist;
            b_save = b;
            bub_save = bubble;
            tmp_depth2_save = tmp_depth2;
            while (layer>0) {
              layer--;
              tmp_depth2 -= temp.depth[layer];
              b = 0.5*(temp.expt[layer]-3);
              bubble = temp.bubble[layer];
              tmp_resid_moist = temp.resid_moist[layer]*temp.depth[layer]*1000;
              zwt_prime_eff = tmp_depth2_save*100-bubble+bubble*pow((zwt_prime+bub_save-tmp_depth2_save*100)/bub_save,b/b_save);
              w_avg = -(b/(b-1))*bubble*( 1 - pow((zwt_prime_eff+bubble-tmp_depth2*100)/bubble,(b-1)/b) ) / (temp.depth[layer]*100);
              tmp_moist += w_avg*(temp.max_moist[layer]-tmp_resid_moist)+tmp_resid_moist;
              b_save = b;
              bub_save = bubble;
              tmp_depth2_save = tmp_depth2;
            }
            temp.zwtvmoist_moist[options.Nlayer+1][i] = tmp_moist;
          }
          zwt_prime += tmp_depth*100/(MAX_ZWTVMOIST-1); // in cm
        }

        /* Compute soil albedo in PAR range (400-700nm) following eqn 122 in Knorr 1997 */
        if (options.CARBON) {
          temp.AlbedoPar = 0.92 * BARE_SOIL_ALBEDO - 0.015;
          if (temp.AlbedoPar < AlbSoiParMin) temp.AlbedoPar = AlbSoiParMin;
        }

      } /* !OUTPUT_FORCE */

      /*************************************************
        Miscellaneous terms for MTCLIM disaggregation
      *************************************************/
      /* Central Longitude of Current Time Zone */
      temp.time_zone_lng = off_gmt * 360./24.;
      /* Assume flat grid cell for radiation calculations */
      temp.slope = 0;
      temp.aspect = 0;
      temp.whoriz = 0;
      temp.ehoriz = 0;

    } // end if(!(*MODEL_DONE) && (*RUN_MODEL))

  return temp;

}


