#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

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
  double expt[MAXlayer];         N/A     pore-size distribution, HBH 5.15
  double Ksat[MAXlayer];         mm/day  saturated hydraulic  conductivity
  double phi_s[MAXlayer];        mm/mm   saturated matrix potential
  double init_moist[MAXlayer];   mm      initial layer moisture level
  float  elevation;	         m       grid cell elevation
  double depth[MAXlayer];        m       thickness of each layer
  double avg_temp;	         C       average soil temperature
  double dp;		         m       soil thermal damping depth
  double bubble;	         cm      bubbling pressure, HBH 5.15
  double quartz;	         fract   quartz content of soil
  double bulk_density[MAXlayer]; kg/m^3  soil bulk density
  double soil_density;		 kg/m^3  soil partical density
  double rough;		         m       soil surface roughness
  double snow_rough;             m       snow surface roughness

  Parameters Computed from those in the File:
  TYPE   NAME                    UNITS   DESCRIPTION
  double max_moist[MAXlayer];    mm      maximum moisture content per layer
  double max_infil;	         N/A     maximum infiltration rate
  double Wcr[MAXlayer];	         mm      critical moisture level for soil
                                         layer, evaporation is no longer
                                         affected moisture stress in the soil
  double Wpwp[MAXlayer];         mm      soil moisture content at permanent
                                         wilting point
  float  time_zone_lng;	         degrees central meridian of the time zone


  Modifications:
  7-19-96	Modified to read through variable layers, and
		read soil depth and average temperature for
		full energy and frozen soil versions of the
		model.						KAC
  4-12-98       Modified to read all parameters from a single
                standard input file.                            KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  static double *lat;
  static double *lng;
  static int    *cellnum;

  int             layer;
  int             cnt;
  soil_con_struct temp; 
  double          Wcr_FRACT;
  double          Wpwp_FRACT;
  double          off_gmt;
  double         *clay;
  double         *sand;
  double         *porosity;
  double          sum_depth;
  char            errstr[MAXSTRING];
  char            namestr[MAXSTRING];
  char            tmpstr[MAXSTRING];

  if(cell==0) {
    rewind(soilparam);
    cnt = 0;
    while(!feof(soilparam)) {
      fscanf(soilparam,"%s",tmpstr);
      cnt++;
    }
    if(cnt!=16+8*options.Nlayer) {
      sprintf(errstr,"Not the right number of soil parameter files in the ARC/INFO file list.");
      vicerror(errstr);
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

    /** Get Critical Soil Moisture Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    Wcr_FRACT = read_arcinfo_value(namestr,temp.lat,temp.lng);

    /** Get Wilting Point Soil Moisture Fraction **/
    fscanf(soilparam,"%s",tmpstr);
    strcpy(namestr,soilparamdir);
    strcat(namestr,tmpstr);
    Wpwp_FRACT = read_arcinfo_value(namestr,temp.lat,temp.lng);

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

    sand     = (double *)calloc(options.Nlayer,sizeof(double));
    clay     = (double *)calloc(options.Nlayer,sizeof(double));
    porosity = (double *)calloc(options.Nlayer,sizeof(double));

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
      if(temp.depth[layer] < MINSOILDEPTH) {
	fprintf(stderr,"Model will not function with layer depth %f < %f m.\n",
		temp.depth[layer],MINSOILDEPTH);
	exit(0);
      }
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
      porosity[layer] = read_arcinfo_value(namestr,temp.lat,temp.lng);
    }

    /*******************************************
      Compute Soil Layer Properties
      *******************************************/
    sum_depth = 0.;
    temp.soil_density = 0.;
    temp.quartz = 0.;
    temp.bubble = 0.;
    for(layer=0;layer<options.Nlayer;layer++) {
      sum_depth += temp.depth[layer];
      temp.bulk_density[layer] *= 1000.;
      temp.soil_density += temp.bulk_density[layer] / (1.0 - porosity[layer]) 
	                 * temp.depth[layer];
      temp.quartz += sand[layer]/100. * temp.depth[layer];
      temp.max_moist[layer] = temp.depth[layer] * porosity[layer] * 1000.;
      if(temp.init_moist[layer] > temp.max_moist[layer]) 
	temp.init_moist[layer] = temp.max_moist[layer];
      temp.bubble += exp(5.3396738 + 0.1845038*clay[layer] 
			 - 2.48394546*porosity[layer] 
			 - 0.00213853*pow(clay[layer],2.)
			 - 0.04356349*sand[layer]*porosity[layer]
			 - 0.61745089*clay[layer]*porosity[layer]
			 + 0.00143598*pow(sand[layer],2.)
			 * pow(porosity[layer],2.)
			 - 0.00855375*pow(clay[layer],2.)
			 * pow(porosity[layer],2.)
			 - 0.00001282*pow(sand[layer],2.)*clay[layer]
			 + 0.00895359*pow(clay[layer],2.)*porosity[layer]
			 - 0.00072472*pow(sand[layer],2.)*porosity[layer]
			 + 0.00000540*pow(clay[layer],2.)*sand[layer]
			 + 0.50028060*pow(porosity[layer],2.)*clay[layer]) 
	             * temp.depth[layer];
      temp.expt[layer] = exp(-0.7842831 + 0.0177544*sand[layer] 
			     - 1.062498*porosity[layer] 
			     - 0.00005304*pow(sand[layer],2.)
			     - 0.00273493*pow(clay[layer],2.)
			     + 1.11134946*pow(porosity[layer],2.)
			     - 0.03088295*sand[layer]*porosity[layer]
			     + 0.00026587*pow(sand[layer],2.)
			     * pow(porosity[layer],2.)
			     - 0.00610522*pow(clay[layer],2.)
			     * pow(porosity[layer],2.)
			     - 0.00000235*pow(sand[layer],2.)*clay[layer]
			     + 0.00798746*pow(clay[layer],2.)*porosity[layer]
			     - 0.00674491*pow(porosity[layer],2.)*clay[layer]);
      temp.expt[layer] = 2. / temp.expt[layer] + 3.;
      temp.resid_moist[layer] = - 0.0182482 + 0.00087269 * sand[layer]
	                      + 0.00513488 * clay[layer] 
	                      + 0.02939286 * porosity[layer] 
                              - 0.00015395 * pow(clay[layer],2.) 
                              - 0.00108270 * sand[layer] * porosity[layer] 
                              - 0.00018233 * pow(clay[layer],2.) 
                                           * pow(porosity[layer],2.) 
                              + 0.00030703 * pow(clay[layer],2.0) 
                                           * porosity[layer] 
                              - 0.00235840 * pow(porosity[layer],2.) 
                                           * clay[layer];
    }
    temp.bubble /= sum_depth;
    temp.soil_density /= sum_depth;
    temp.quartz /= sum_depth;
    /** Need to fix this in the future so that dp can be set to any depth **/
    if(temp.dp > sum_depth) temp.dp = sum_depth;

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
      temp.Wcr[layer]  = Wcr_FRACT * temp.max_moist[layer];
      temp.Wpwp[layer] = Wpwp_FRACT * temp.max_moist[layer];
      if(temp.Wpwp[layer] > temp.Wcr[layer])
        nrerror("Wpwp is greater then Wcr");
    }

    /*************************************************
      Determine Central Longitude of Current Time Zone 
      *************************************************/
    temp.time_zone_lng = off_gmt * 360./24.;

    free((char *)sand);
    free((char *)clay);
    free((char *)porosity);

  }
  else RUN[0] = 0;

  return temp;
} 
