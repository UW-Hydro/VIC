#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

soil_con_struct read_soilparam(FILE *soilparam,
			       int   RUN_MODEL)
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
  double c;                      N/A     exponent
  double expt[MAX_LAYERS];         N/A     exponent in Campbell's eqn for hydraulic conductivity, HBH 5.6
  double Ksat[MAX_LAYERS];         mm/day  saturated hydraulic  conductivity
  double phi_s[MAX_LAYERS];        mm/mm   saturated matrix potential
  double init_moist[MAX_LAYERS];   mm      initial layer moisture level
  float  elevation;	         m       grid cell elevation
  double depth[MAX_LAYERS];        m       thickness of each layer
  double avg_temp;	         C       average soil temperature
  double dp;		         m       soil thermal damping depth
  double bubble;	         cm      bubbling pressure, HBH 5.15
  double quartz;	         fract   quartz content of soil
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
**********************************************************************/
{
  void ttrim( char *string );
  extern option_struct options;
  extern global_param_struct global_param;
  extern veg_lib_struct *veg_lib;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

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
#if EXCESS_ICE
  double          init_ice_fract[MAX_LAYERS];
#endif
  soil_con_struct temp; 

  if( fgets( line, MAXSTRING, soilparam ) == NULL ){
    sprintf(ErrStr,"ERROR: Unexpected EOF while reading soil file\n");
    nrerror(ErrStr);
  }

  if ( RUN_MODEL ) {

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
#if !OUTPUT_FORCE
      if(temp.expt[layer] < 3.0) {
	fprintf(stderr,"ERROR: Exponent in layer %d is %f < 3.0; This must be > 3.0\n", layer, temp.expt[layer]);
	exit(0);
      }
#endif /* !OUTPUT_FORCE */
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
#if !OUTPUT_FORCE
      if(temp.init_moist[layer] < 0.) {
	sprintf(ErrStr,"ERROR: Initial moisture for layer %d cannot be negative (%f)",layer,temp.init_moist[layer]);
	nrerror(ErrStr);
      }
#endif /* !OUTPUT_FORCE */
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

    /* final soil layer thicknesses for !EXCESS_ICE option */
#if !EXCESS_ICE
#if !OUTPUT_FORCE
    for(layer = 0; layer < options.Nlayer; layer++) 
      temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
#endif
#endif /* !EXCESS_ICE */    
    
    /* read average soil temperature */
    token = strtok (NULL, delimiters);
    while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
    if( token == NULL ) {
      sprintf(ErrStr,"ERROR: Can't find values for AVERAGE SOIL TEMPERATURE in soil file\n");
      nrerror(ErrStr);
    }  
    sscanf(token, "%lf", &temp.avg_temp);
#if !OUTPUT_FORCE
    if(options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
      fprintf(stderr,"Need valid average soil temperature in degrees C to run");
      fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	      temp.avg_temp);
      exit(0);
    }
#endif /* !OUTPUT_FORCE */
    
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
#if !OUTPUT_FORCE
      if((options.FULL_ENERGY || options.FROZEN_SOIL) && temp.bubble[layer] < 0) {
	fprintf(stderr,"ERROR: Bubbling pressure in layer %d is %f < 0; This must be positive for FULL_ENERGY = TRUE or FROZEN_SOIL = TRUE\n", layer, temp.bubble[layer]);
	exit(0);
      }
#endif /* !OUTPUT_FORCE */
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
#if !OUTPUT_FORCE
      if(options.FULL_ENERGY 
	 && (temp.quartz[layer] > 1. || temp.quartz[layer] < 0)) {
	fprintf(stderr,"Need valid quartz content as a fraction to run");
	fprintf(stderr," Full Energy model, %f is not acceptable.\n",
		temp.quartz[layer]);
	exit(0);
      }
#endif /* !OUTPUT_FORCE */
    }
    
    /* read layer bulk density */
    for(layer = 0; layer < options.Nlayer; layer++){
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for BULK DENSITY for layer %d in soil file\n", layer );
        nrerror(ErrStr);
      }  
      sscanf(token, "%lf", &temp.bulk_density[layer]);
    }

    /* read layer soil density */
    for(layer = 0; layer < options.Nlayer; layer++) {
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
        sprintf(ErrStr,"ERROR: Can't find values for SOIL DENSITY for layer %d in soil file\n", layer );
        nrerror(ErrStr);
      }  
      sscanf(token, "%lf", &temp.soil_density[layer]);
#if !OUTPUT_FORCE
      if(temp.bulk_density[layer]>=temp.soil_density[layer])
	nrerror("Layer bulk density must be less then soil density");
#endif /* !OUTPUT_FORCE */
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
    /* Overwrite default bare soil aerodynamic resistance parameters
       with the values taken from the soil parameter file */
    for (j=0; j<12; j++) {
      veg_lib[veg_lib[0].NVegLibTypes].roughness[j] = temp.rough;
      veg_lib[veg_lib[0].NVegLibTypes].displacement[j] = temp.rough*0.667/0.123;
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
#if SPATIAL_SNOW
    token = strtok (NULL, delimiters);
    while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
    if( token == NULL ) {
      sprintf(ErrStr,"ERROR: Can't find values for SPATIAL SNOW in soil file\n");
      nrerror(ErrStr);
    }
    sscanf(token, "%lf", &tempdbl);
    temp.depth_full_snow_cover = tempdbl;
#endif // SPATIAL_SNOW
    
    /* read slope of frozen soil distribution */
#if SPATIAL_FROST
    token = strtok (NULL, delimiters);
    while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
    if( token == NULL ) {
      sprintf(ErrStr,"ERROR: Can't find values for SPATIAL FROST in soil file\n");
      nrerror(ErrStr);
    }  
    sscanf(token, "%lf", &tempdbl);
    temp.frost_slope = tempdbl;
#endif // SPATIAL_FROST
    
    /*read volumetric ice fraction for each soil layer */
#if EXCESS_ICE
    for(layer = 0; layer < options.Nlayer; layer++) {
      token = strtok (NULL, delimiters);
      while (token != NULL && (length=strlen(token))==0) token = strtok (NULL, delimiters);
      if( token == NULL ) {
	sprintf(ErrStr,"ERROR: Can't find values for VOLUMETRIC ICE FRACTION (EXCESS_ICE = TRUE) for layer %d in soil file\n", layer);
	nrerror(ErrStr);
      }  
      sscanf(token, "%lf", &init_ice_fract[layer]);
    }
#endif // EXCESS_ICE
 
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

#if !OUTPUT_FORCE
    /*******************************************
      Compute Soil Layer Properties
    *******************************************/
    for(layer = 0; layer < options.Nlayer; layer++) {
      if (temp.resid_moist[layer] == MISSING)
	temp.resid_moist[layer] = RESID_MOIST;
      temp.porosity[layer] = 1.0 - temp.bulk_density[layer] 
	/ temp.soil_density[layer];
#if !EXCESS_ICE      
      temp.max_moist[layer] = temp.depth[layer] * temp.porosity[layer] * 1000.;
#endif      
    }

#if !EXCESS_ICE
    /*******************************************
      Validate Initial Soil Layer Moisture Content for !EXCESS_ICE option.
    *******************************************/
    if (!options.INIT_STATE) { // only do this if we're not getting initial moisture from model state file
      for(layer = 0; layer < options.Nlayer; layer++) {
	if(temp.init_moist[layer] > temp.max_moist[layer]) {
	  fprintf(stderr,"Initial soil moisture (%f mm) is greater than the maximum moisture (%f mm) for layer %d.\n\tResetting soil moisture to maximum.\n",
		temp.init_moist[layer], temp.max_moist[layer], layer);
	  temp.init_moist[layer] = temp.max_moist[layer];
	}
	if(temp.init_moist[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) { 
	  fprintf(stderr,"Initial soil moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %d.\n\tResetting soil moisture to residual moisture.\n",
		temp.init_moist[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
	  temp.init_moist[layer] = temp.resid_moist[layer] * temp.depth[layer] * 1000.;
        }
      }
    }
#endif
    
#if EXCESS_ICE
    /*******************************************
      Compute Soil Layer Properties for EXCESS_ICE option
    *******************************************/
    extra_depth=0;
    for(layer = 0; layer < options.Nlayer; layer++) {
      temp.min_depth[layer]=temp.depth[layer];
      if(init_ice_fract[layer]>MAX_ICE_INIT){ // validate amount based on physical constraints
	fprintf(stderr,"Initial ice fraction (%f) is greater than maximum ice content for layer %d.\n\tResetting to maximum of %f\n",init_ice_fract[layer],layer,MAX_ICE_INIT);
	init_ice_fract[layer]=MAX_ICE_INIT;
      }
      if(init_ice_fract[layer]>=temp.porosity[layer]){//excess ground ice present
	fprintf(stderr,"Excess ground ice present in layer %d:\n",layer+1);
	fprintf(stderr,"\t\tSubsidence will occur when the average soil layer\n\t\t  temperature exceeds %.2f degrees Celsius.\n",powf((1.-ICE_AT_SUBSIDENCE),(3.-temp.expt[layer])/2.)*273.16*9.81*temp.bubble[layer]/(-Lf*100.));	
	temp.effective_porosity[layer]=init_ice_fract[layer];
	fprintf(stderr,"\t\tEffective porosity increased from %.2f to %.2f.\n",temp.porosity[layer],temp.effective_porosity[layer]);
	temp.depth[layer] = temp.min_depth[layer]*(1.0 - temp.porosity[layer])/(1.0 - temp.effective_porosity[layer]); //adjust soil layer depth
	extra_depth += temp.depth[layer]-temp.min_depth[layer]; //net increase in depth due to excess ice
	fprintf(stderr,"\t\tDepth of soil layer adjusted for excess ground ice: from %.2f m to %.2f m.\n",temp.min_depth[layer],temp.depth[layer]);
	fprintf(stderr,"\t\tBulk density adjusted for excess ground ice: from %.2f kg/m^3 to %.2f kg/m^3.\n",temp.bulk_density[layer],(1.0-temp.effective_porosity[layer])*temp.soil_density[layer]);
	temp.bulk_density[layer] = (1.0-temp.effective_porosity[layer])*temp.soil_density[layer]; //adjust bulk density
      }
      else //excess ground ice not present
	temp.effective_porosity[layer]=temp.porosity[layer];
    }
    if(extra_depth>0) {
      fprintf(stderr,"Damping depth adjusted for excess ground ice: from %.2f m to %.2f m.\n",temp.dp,temp.dp+extra_depth);
      temp.dp += extra_depth;  //adjust damping depth
    }
    
    /* final soil layer thicknesses for EXCESS_ICE option */
    for(layer = 0; layer < options.Nlayer; layer++) 
      temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
    
    /* Calculate and Validate Maximum Initial Soil Layer Moisture Content for EXCESS_ICE option */
    for(layer = 0; layer < options.Nlayer; layer++) {
      temp.max_moist[layer] = temp.depth[layer] * temp.effective_porosity[layer] * 1000.;
      if(temp.effective_porosity[layer]>temp.porosity[layer])//excess ground ice present
	temp.init_moist[layer] = temp.max_moist[layer]; 
      else {//excess ground ice not present
	if(temp.depth[layer] * init_ice_fract[layer] * 1000. > temp.init_moist[layer])
	  temp.init_moist[layer] = temp.depth[layer] * init_ice_fract[layer] * 1000.;
      }
    }
    for(layer = 0; layer < options.Nlayer; layer++) {
      if(temp.init_moist[layer] > temp.max_moist[layer]) {
	fprintf(stderr,"Initial soil moisture (%f mm) is greater than the maximum moisture (%f mm) for layer %d.\n\tResetting soil moisture to maximum.\n",
		temp.init_moist[layer], temp.max_moist[layer], layer);
	temp.init_moist[layer] = temp.max_moist[layer];
      }
      if(temp.init_moist[layer] < temp.resid_moist[layer] * temp.depth[layer] * 1000.) { 
	fprintf(stderr,"Initial soil moisture (%f mm) is less than calculated residual moisture (%f mm) for layer %d.\n\tResetting soil moisture to residual moisture.\n",
		temp.init_moist[layer], temp.resid_moist[layer] * temp.depth[layer] * 1000., layer);
	temp.init_moist[layer] = temp.resid_moist[layer] * temp.depth[layer] * 1000.;
      }
    }
#endif // EXCESS_ICE   
    
    /**********************************************
      Validate soil layer depths for top two layers
    **********************************************/
#if !OUTPUT_FORCE
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
#if EXCESS_ICE
    for(layer = 0; layer < options.Nlayer; layer++) {
      if(temp.min_depth[layer] < MINSOILDEPTH) {
	sprintf(ErrStr,"ERROR: Model will not function with layer %d depth %f < %f m.\n",
		layer,temp.min_depth[layer],MINSOILDEPTH);
	nrerror(ErrStr);
      }
    }
    if(temp.min_depth[0] > temp.min_depth[1]) {
      sprintf(ErrStr,"ERROR: Model will not function with layer %d depth (%f m) > layer %d depth (%f m).\n",
	      0,temp.min_depth[0],1,temp.min_depth[1]);
      nrerror(ErrStr);
    }
#endif /* EXCESS_ICE */
#endif /* !OUTPUT_FORCE */

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
#if EXCESS_ICE
      temp.Wcr_FRACT[layer]  = Wcr_FRACT[layer];
      temp.Wpwp_FRACT[layer] = Wpwp_FRACT[layer]; 
#endif
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
    }    
    
    /**********************************************
      Validate Spatial Snow/Frost Params
    **********************************************/
#if SPATIAL_SNOW
    if (temp.depth_full_snow_cover < 0.0) {
      sprintf(ErrStr,"depth_full_snow_cover (%f) must be positive.\n", temp.depth_full_snow_cover);
      nrerror(ErrStr);
    }
#endif // SPATIAL_SNOW
    
#if SPATIAL_FROST
    if (temp.frost_slope < 0.0) {
      sprintf(ErrStr,"frost_slope (%f) must be positive.\n", temp.frost_slope);
      nrerror(ErrStr);
    }
#endif // SPATIAL_FROST
    
    
    /*************************************************
      If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
      parameters d1, d2, d3, and d4 to ARNO baseflow
      parameters Ds, Dsmax, Ws, and c
    *************************************************/
#if EXCESS_ICE
    temp.Dsmax_orig = temp.Dsmax;
    temp.Ds_orig = temp.Ds;
    temp.Ws_orig = temp.Ws;
#endif
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

#endif /* !OUTPUT_FORCE */
    
    
    /*************************************************
      Determine Central Longitude of Current Time Zone 
    *************************************************/
    temp.time_zone_lng = off_gmt * 360./24.;
    
    /* Allocate Layer - Node fraction array */
    temp.layer_node_fract = (float **)malloc((options.Nlayer+1)*sizeof(float *));
    for(layer=0;layer<=options.Nlayer;layer++) 
      temp.layer_node_fract[layer] = (float *)malloc(options.Nnode*sizeof(float));

}
  /* ELSE Grid cell is not active (RUN_MODEL=0), skip soil parameter data */
  
  return temp;
  
} 


