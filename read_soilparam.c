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
  3-13-00       Modified to read more parameters as separate
                layer values                                    KAC
  6-6-2000      Modified to skip individual parameter reads
                if model grid cell is not read.                 KAC

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  char            ErrStr[MAXSTRING];
  char            NotRunStr[2048];
  int             layer, i, tempint;
  double          Wcr_FRACT[MAX_LAYERS];
  double          Wpwp_FRACT[MAX_LAYERS];
  double          off_gmt;
  double          tmp;
  double          tempdbl;
  soil_con_struct temp; 

  if ( RUN_MODEL ) {

    fscanf(soilparam, "%d", &temp.gridcel);
    fscanf(soilparam, "%f", &temp.lat);
    fscanf(soilparam, "%f", &temp.lng);
#if VERBOSE
    /* add print statements for grid cell number -- EDM */
    fprintf(stderr,"\ncell: %d,  lat: %.4f, long: %.4f\n",temp.gridcel,temp.lat,temp.lng);
#endif
    
    /* read infiltration parameter */
    fscanf(soilparam, "%lf", &temp.b_infilt);
    
    /* read fraction of baseflow rate */
    fscanf(soilparam, "%lf", &temp.Ds);
    
    /* read maximum baseflow rate */
    fscanf(soilparam, "%lf", &temp.Dsmax);
    
    /* read fraction of bottom soil layer moisture */
    fscanf(soilparam, "%lf", &temp.Ws);
    
    /* read exponential */
    fscanf(soilparam, "%lf", &temp.c);
    
    /* read expt for each layer */
    for(layer = 0; layer < options.Nlayer; layer++) 
      fscanf(soilparam, "%lf", &temp.expt[layer]);
    
    /* read layer saturated hydraulic conductivity */
    for(layer = 0; layer < options.Nlayer; layer++)
      fscanf(soilparam, "%lf", &temp.Ksat[layer]);
    
    /* read layer phi_s */
    for(layer = 0; layer < options.Nlayer; layer++)
      fscanf(soilparam, "%lf", &temp.phi_s[layer]);
    
    /* read layer initial moisture */
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam, "%lf", &temp.init_moist[layer]);
      if(temp.init_moist[layer] < 0.) {
	sprintf(ErrStr,"ERROR: Initial moisture for layer %i cannot be negative (%f)",layer,temp.init_moist[layer]);
	nrerror(ErrStr);
      }
    }
    
    /* read cell mean elevation */
    fscanf(soilparam, "%f", &temp.elevation);
    
    /* soil layer thicknesses */
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam, "%lf", &temp.depth[layer]);
      temp.depth[layer] = rint(temp.depth[layer] * 1000) / 1000;
      if(temp.depth[layer] < MINSOILDEPTH) {
	sprintf(ErrStr,"ERROR: Model will not function with layer %i depth %f < %f m.\n",
		layer,temp.depth[layer],MINSOILDEPTH);
	nrerror(ErrStr);
      }
    }
    if(temp.depth[0] > temp.depth[1]) {
      sprintf(ErrStr,"ERROR: Model will not function with layer %i depth (%f m) < layer %i depth (%f m).\n",
	      0,temp.depth[0],1,temp.depth[1]);
      nrerror(ErrStr);
    }
    
    /* read average soil temperature */
    fscanf(soilparam, "%lf", &temp.avg_temp);
    if(options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
      fprintf(stderr,"Need valid average soil temperature in degrees C to run");
      fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	      temp.avg_temp);
      exit(0);
    }
    
    /* read soil damping depth */
    fscanf(soilparam, "%lf", &temp.dp);
    
    /* read layer bubbling pressure */
    for(layer = 0; layer < options.Nlayer; layer++) 
      fscanf(soilparam, "%lf", &temp.bubble[layer]);
    
    /* read layer quartz content */
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam, "%lf", &temp.quartz[layer]);
      if(options.FULL_ENERGY 
	 && (temp.quartz[layer] > 1. || temp.quartz[layer] < 0)) {
	fprintf(stderr,"Need valid quartz content as a fraction to run");
	fprintf(stderr," Full Energy model, %f is not acceptable.\n",
		temp.quartz[layer]);
	exit(0);
      }
    }
    
    /* read layer bulk density */
    for(layer = 0; layer < options.Nlayer; layer++)
      fscanf(soilparam, "%lf", &temp.bulk_density[layer]);
    
    /* read layer soil density */
    for(layer = 0; layer < options.Nlayer; layer++) {
      fscanf(soilparam, "%lf", &temp.soil_density[layer]);
      if(temp.bulk_density[layer]>=temp.soil_density[layer])
	nrerror("Layer bulk density must be less then soil density");
    }
    
    /* read cell gmt offset */
    fscanf(soilparam, "%lf", &off_gmt);
    
    /* read layer critical point */
    for(layer=0;layer<options.Nlayer;layer++)
      fscanf(soilparam, "%lf", &(Wcr_FRACT[layer]));
    
    /* read layer wilting point */
    for(layer=0;layer<options.Nlayer;layer++)
      fscanf(soilparam, "%lf", &(Wpwp_FRACT[layer]));
    
    /* read soil roughness */
    fscanf(soilparam, "%lf", &temp.rough);
    
    /* read snow roughness */
    fscanf(soilparam, "%lf", &temp.snow_rough);
    
    /* read cell annual precipitation */
    fscanf(soilparam, "%lf", &temp.annual_prec);
    
    /* read layer residual moisture content */
    for(layer = 0; layer < options.Nlayer; layer++) 
      fscanf(soilparam, "%lf", &temp.resid_moist[layer]);
    
    /* read frozen soil active flag */
    fscanf(soilparam, "%i", &tempint);
    temp.FS_ACTIVE = (char)tempint;
    
    /* read minimum snow depth for full coverage */
    fscanf(soilparam, "%lf", &tempdbl);
#if SPATIAL_SNOW
    temp.depth_full_snow_cover = tempdbl;
#endif // SPATIAL_SNOW
    
    /* read slope of frozen soil distribution */
    fscanf(soilparam, "%lf", &tempdbl);
#if SPATIAL_FROST
    temp.frost_slope = tempdbl;
#endif // SPATIAL_FROST
    
    /*******************************************
      Compute Maximum Soil Layer Moiture Content
    *******************************************/
    for(layer = 0; layer < options.Nlayer; layer++) {
      if (temp.resid_moist[layer] == MISSING)
	temp.resid_moist[layer] = RESID_MOIST;
      temp.porosity[layer] = 1.0 - temp.bulk_density[layer] 
	/ temp.soil_density[layer];
      if(temp.porosity[layer] < temp.resid_moist[layer]) {
	sprintf(ErrStr,"Layer %i porosity (%f mm/mm) must be greater then the defined moisture residue of %f mm/mm", 
		layer, temp.porosity[layer], temp.resid_moist[layer]);
	nrerror(ErrStr);
      }
      temp.max_moist[layer] = temp.depth[layer] * temp.porosity[layer] * 1000.;
      if(temp.init_moist[layer] > temp.max_moist[layer]) 
	temp.init_moist[layer] = temp.max_moist[layer];
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
      if(temp.Wpwp[layer] > temp.Wcr[layer])
	nrerror("Wpwp is greater then Wcr");
      if(temp.init_moist[layer] < temp.Wpwp[layer]) { 
	fprintf(stderr,"Initial soil moisture (%f) is less than the wilting point (%f) for layer %i\n\tResetting soil moisture to wilting point\n",
		temp.init_moist[layer], temp.Wpwp[layer], layer);
	temp.init_moist[layer] = temp.Wpwp[layer];
      }
    }

    /*************************************************
      Determine Central Longitude of Current Time Zone 
    *************************************************/
    temp.time_zone_lng = off_gmt * 360./24.;

    /* Allocate Layer - Node fraction array */
    temp.layer_node_fract = (float **)malloc((options.Nlayer+1)*sizeof(float *));
    for(layer=0;layer<=options.Nlayer;layer++) 
      temp.layer_node_fract[layer] = (float *)malloc(options.Nnode*sizeof(float));

  }
  else {
   
    /* Grid cell is not active, skip soil parameter data */
    fgets(NotRunStr, 2048, soilparam);

  }
    
  return temp;

} 
