#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

soil_con_struct read_soilparam_arc(FILE *soilparam, 
				   char *soilparamdir, 
				   int  *Ncells,
				   int  *RUN,
				   int   cell)
/**********************************************************************
	read_soilparam_arc     Keith Cherkauer		May 5, 1998

  This routine reads soil parameters for each grid cell from an ASCII
  ARC/INFO output grid.

  Order of ARC/INFO Files:
  CELLNUM
  RUN
  ELEVATION
  B_INFILT
  Ds
  DsMax
  Ws
  c
  AVG_TEMP
  DP
  OFF_GMT
  Wcr_FRACT
  Wpwp_FRACT
  ROUGH
  SNOW_ROUGH
  SAND[]
  CLAY[]
  KSAT[]
  PHI_S[]
  INIT_MOIST[]
  DEPTH[]
  BULK_DENSITY[]
  POROSITY[]

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
  double expt[MAX_LAYERS];         N/A     pore-size distribution, HBH 5.15
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
  07-May-04	Replaced rint(something) with (float)(int)(something + 0.5)
		to handle rounding without resorting to rint().	TJB
  10-May-04	Modified to handle Arno parameters.		TJB
  11-May-04	Removed extraneous tmp variable.		TJB

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  static double *lat;
  static double *lng;
  static int    *cellnum;

  int             layer;
  int             cnt;
  soil_con_struct temp; 
  double          Wcr_FRACT[MAX_LAYERS];
  double          Wpwp_FRACT[MAX_LAYERS];
  double          off_gmt;
  double          clay[MAX_LAYERS];
  double          sand[MAX_LAYERS];
  double          sum_depth;
  char            ErrStr[MAXSTRING];
  char            namestr[MAXSTRING];
  char            tmpstr[MAXSTRING];

  double tmp_bubble;

  tmp_bubble = 0;

  if(cell==0) {
    rewind(soilparam);

    cnt = 0;
    while(!feof(soilparam)) {
      fscanf(soilparam,"%s",tmpstr);
      cnt++;
    }
    if(cnt!=16+10*options.Nlayer) {
      sprintf(ErrStr,"Not the right number of soil parameter files in the ARC/INFO file list.");
      nrerror(ErrStr);
    }
    
    rewind(soilparam);
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    Ncells[0] = read_arcinfo_info(namestr,&lat,&lng,&cellnum);
  }
  else {
    rewind(soilparam);
    fscanf(soilparam,"%s",tmpstr);
  }

  temp.gridcel = cellnum[cell];
  temp.lat = lat[cell];
  temp.lng = lng[cell];

  /** Check if Grid Cell is Run in Model **/
  fscanf(soilparam,"%s",tmpstr);
  strcpy(namestr,soilparamdir);
  strcat(namestr,tmpstr);
  *RUN = (int)read_arcinfo_value(namestr,temp.lat,temp.lng);
  
  if(RUN[0] > 0) {
    /** Get Average Grid Cell Elevation **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.elevation = (float)read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Grid Cell Infiltration Parameter **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.b_infilt = read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Baseflow Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.Ds = read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Baseflow Velocity **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.Dsmax = read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Maximum Soil Moisture Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.Ws = read_arcinfo_value(namestr,temp.lat,temp.lng);
  
    /** Get Exponential **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.c = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Average Soil Temperature **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.avg_temp = read_arcinfo_value(namestr,temp.lat,temp.lng);
    if(options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
      fprintf(stderr,"Need valid average soil temperature in degrees C to run");
      fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	      temp.avg_temp);
      exit(0);
    }

    /** Get Soil Thermal Damping Depth **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.dp = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Data Time Zone Offset from GMT **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    off_gmt = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Critical Soil Moisture Fraction for each layer **/
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      Wcr_FRACT[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Wilting Point Soil Moisture Fraction for each layer **/
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      Wpwp_FRACT[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Bare Soil Roughness **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.rough = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Snow Surface Roughness **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.snow_rough = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Average Annual Precipitation **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.annual_prec = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Layer Percent Sand **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      sand[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Percent Clay **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      clay[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Saturated Hydraulic Conductivity **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.Ksat[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Soil Moisture Diffusion Parameter **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.phi_s[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Initial Layer Moisture **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.init_moist[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Thickness **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.depth[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
      temp.depth[layer] = (float)(int)(temp.depth[layer] * 1000 + 0.5) / 1000;
      if(temp.depth[layer] < MINSOILDEPTH) {
	fprintf(stderr,"Model will not function with layer depth %f < %f m.\n",
		temp.depth[layer],MINSOILDEPTH);
	exit(0);
      }
    }
    if(temp.depth[0] > temp.depth[1]) {
      sprintf(ErrStr,"ERROR: Model will not function with layer %i depth (%f m) < layer %i depth (%f m).\n",
	      0,temp.depth[0],1,temp.depth[1]);
      nrerror(ErrStr);
    }

    /** Get Layer Bulk Density **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.bulk_density[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Get Layer Porosity **/
    for(layer=0;layer<options.Nlayer;layer++) {
      fscanf(soilparam,"%s",tmpstr);
      strcpy(namestr,soilparamdir);
      strcat(namestr,tmpstr);
      temp.porosity[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /** Activate Frozen Soil for Grid Cell **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    temp.FS_ACTIVE = (char)read_arcinfo_value(namestr,temp.lat,temp.lng);

    /*************************************************
    if ARNO_PARAMS == TRUE then convert the baseflow
    parameters d1, d2, d3, d4 to Ds, Dsmax, Ws, and c.  JA
    *************************************************/
    if(options.ARNO_PARAMS) {
      layer = options.Nlayer-1;
      temp.Dsmax = temp.Dsmax *
        pow((double)(1./(temp.max_moist[layer]-temp.Ws)), -temp.c) +
        temp.Ds * temp.max_moist[layer];
      temp.Ds = temp.Ds * temp.Ws / temp.Dsmax;
      temp.Ws = temp.Ws/temp.max_moist[layer];
    }

    /*******************************************
      Compute Soil Layer Properties
      *******************************************/
    sum_depth   = 0.;
    for(layer = 0; layer < options.Nlayer; layer++) {
      sum_depth += temp.depth[layer];
      temp.bulk_density[layer] *= 1000.;
      temp.soil_density[layer] = temp.bulk_density[layer] 
	/ (1.0 - temp.porosity[layer]);
      temp.quartz[layer] = sand[layer] / 100.;
      temp.max_moist[layer] = temp.depth[layer] * temp.porosity[layer] * 1000.;
      if(temp.init_moist[layer] > temp.max_moist[layer]) 
	temp.init_moist[layer] = temp.max_moist[layer];
      temp.bubble[layer] = exp(5.3396738 + 0.1845038*clay[layer] 
			 - 2.48394546*temp.porosity[layer] 
			 - 0.00213853*pow(clay[layer],2.)
			 - 0.04356349*sand[layer]*temp.porosity[layer]
			 - 0.61745089*clay[layer]*temp.porosity[layer]
			 + 0.00143598*pow(sand[layer],2.)
			 * pow(temp.porosity[layer],2.)
			 - 0.00855375*pow(clay[layer],2.)
			 * pow(temp.porosity[layer],2.)
			 - 0.00001282*pow(sand[layer],2.)*clay[layer]
			 + 0.00895359*pow(clay[layer],2.)*temp.porosity[layer]
			 - 0.00072472*pow(sand[layer],2.)*temp.porosity[layer]
			 + 0.00000540*pow(clay[layer],2.)*sand[layer]
			 + 0.50028060*pow(temp.porosity[layer],2.)*clay[layer]);
      temp.expt[layer] = exp(-0.7842831 + 0.0177544*sand[layer] 
			     - 1.062498*temp.porosity[layer] 
			     - 0.00005304*pow(sand[layer],2.)
			     - 0.00273493*pow(clay[layer],2.)
			     + 1.11134946*pow(temp.porosity[layer],2.)
			     - 0.03088295*sand[layer]*temp.porosity[layer]
			     + 0.00026587*pow(sand[layer],2.)
			     * pow(temp.porosity[layer],2.)
			     - 0.00610522*pow(clay[layer],2.)
			     * pow(temp.porosity[layer],2.)
			     - 0.00000235*pow(sand[layer],2.)*clay[layer]
			     + 0.00798746*pow(clay[layer],2.)*temp.porosity[layer]
			     - 0.00674491*pow(temp.porosity[layer],2.)*clay[layer]);
      temp.expt[layer] = 2. / temp.expt[layer] + 3.;
      temp.resid_moist[layer] = - 0.0182482 + 0.00087269 * sand[layer]
	                      + 0.00513488 * clay[layer] 
	                      + 0.02939286 * temp.porosity[layer] 
                              - 0.00015395 * pow(clay[layer],2.) 
                              - 0.00108270 * sand[layer] * temp.porosity[layer] 
                              - 0.00018233 * pow(clay[layer],2.) 
                                           * pow(temp.porosity[layer],2.) 
                              + 0.00030703 * pow(clay[layer],2.0) 
                                           * temp.porosity[layer] 
                              - 0.00235840 * pow(temp.porosity[layer],2.) 
                                           * clay[layer];

      /** Check for valid values of generated parameters **/
      if(temp.bubble[layer]<1.36) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too low (%f), resetting to minimum value (%f).\n",temp.bubble[layer],1.36);
	temp.bubble[layer] = 1.36;
      }
      if(temp.bubble[layer]>187.2) {
	fprintf(stderr,"WARNING: estimated bubbling pressure too high (%f), resetting to maximum value (%f).\n",temp.bubble[layer],187.2);
	temp.bubble[layer] = 187.2;
      }
      if(temp.expt[layer] < 2. / 1.090 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too low (%f), resetting to minimum value (%f).\n", temp.expt[layer], 2. / 1.090 + 3.);
	temp.expt[layer] = 2. / 1.090 + 3.;
      }
      if(temp.expt[layer] > 2. / 0.037 + 3.) {
	fprintf(stderr,"WARNING: estimated exponential (expt) too high (%f), resetting to maximum value (%f).\n",temp.expt[layer], 2. / 0.037 + 3.);
	temp.expt[layer] = 2. / 0.037 + 3.;
      }
      if(temp.resid_moist[layer] < -0.038) {
	fprintf(stderr,"WARNING: estimated residual soil moisture too low (%f), resetting to minimum value (%f).\n",temp.resid_moist[layer],-0.038);
	temp.resid_moist[layer] = -0.038;
      }
      if(temp.resid_moist[layer] > 0.205) {
	fprintf(stderr,"WARNING: estimated residual soil mositure too high (%f), resetting to maximum value (%f).\n",temp.resid_moist[layer],0.205);
	temp.resid_moist[layer] = 0.205;
      }
      tmp_bubble += temp.bubble[layer];
    }
    for(layer=0;layer<options.Nlayer;layer++) temp.bubble[layer] = tmp_bubble/3.;

    /**********************************************
      Compute Maximum Infiltration for Upper Layers
      **********************************************/
    if(options.Nlayer==2)
      temp.max_infil = (1.0+temp.b_infilt)*temp.max_moist[0];
    else
      temp.max_infil = (1.0+temp.b_infilt)*(temp.max_moist[0]
					    + temp.max_moist[1]);

    /****************************************************************
      Compute Soil Layer Critical and Wilting Point Moisture Contents
      ****************************************************************/
    for(layer=0;layer<options.Nlayer;layer++) {
      temp.Wcr[layer]  = Wcr_FRACT[layer] * temp.max_moist[layer];
      temp.Wpwp[layer] = Wpwp_FRACT[layer] * temp.max_moist[layer];
      if(temp.Wpwp[layer] > temp.Wcr[layer])
        nrerror("Wpwp is greater then Wcr");
      if(temp.init_moist[layer] < temp.Wpwp[layer]) { 
	fprintf(stderr,"Initial soil moisture (%f) is less than the wilting point (%f) for layer %i\n\tResetting soil moisture to equal wilting point\n",temp.init_moist[layer],temp.Wpwp[layer],layer);
	temp.init_moist[layer] = temp.Wpwp[layer];
      }
    }

    /*************************************************
      Determine Central Longitude of Current Time Zone 
      *************************************************/
    temp.time_zone_lng = off_gmt * 360./24.;

  }
  else RUN[0] = 0;

  /* Allocate Layer - Node fraction array */
  temp.layer_node_fract = (float **)malloc((options.Nlayer+1)*sizeof(float *));
  for(layer=0;layer<=options.Nlayer;layer++) 
    temp.layer_node_fract[layer] 
      = (float *)malloc(options.Nnode*sizeof(float));

  return temp;
} 
