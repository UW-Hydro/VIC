#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

soil_con_struct read_soilparam(FILE *soilparam)
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

  int layer;
  soil_con_struct temp; 
  double porosity;
  double max_moist_tot;
  double maxsoilmoist2factor;
  double Wcr_FRACT;
  double Wpwp_FRACT;
  double off_gmt;
  char errstr[MAXSTRING];

  fscanf(soilparam, "%d", &temp.gridcel);
  fscanf(soilparam, "%f", &temp.lat);
  fscanf(soilparam, "%f", &temp.lng);
  fscanf(soilparam, "%lf", &temp.b_infilt);
  fscanf(soilparam, "%lf", &temp.Ds);
  fscanf(soilparam, "%lf", &temp.Dsmax);
  fscanf(soilparam, "%lf", &temp.Ws);
  fscanf(soilparam, "%lf", &temp.c);
  for(layer=0;layer<options.Nlayer;layer++)
    fscanf(soilparam, "%lf", &temp.expt[layer]);
  for(layer=0;layer<options.Nlayer;layer++)
    fscanf(soilparam, "%lf", &temp.Ksat[layer]);
  for(layer=0;layer<options.Nlayer;layer++)
    fscanf(soilparam, "%lf", &temp.phi_s[layer]);
  for(layer=0;layer<options.Nlayer;layer++)
    fscanf(soilparam, "%lf", &temp.init_moist[layer]);
  fscanf(soilparam, "%f", &temp.elevation);
  for(layer=0;layer<options.Nlayer;layer++) {
    fscanf(soilparam, "%lf", &temp.depth[layer]);
    if(temp.depth[layer] < MINSOILDEPTH) {
      fprintf(stderr,"Model will not function with layer depth %f < %f m.\n",
	      temp.depth[layer],MINSOILDEPTH);
      exit(0);
    }
  }
  fscanf(soilparam, "%lf", &temp.avg_temp);
  if(options.FULL_ENERGY && (temp.avg_temp>100. || temp.avg_temp<-50)) {
    fprintf(stderr,"Need valid average soil temperature in degrees C to run");
    fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	    temp.avg_temp);
    exit(0);
  }
  fscanf(soilparam, "%lf", &temp.dp);
  fscanf(soilparam, "%lf", &temp.bubble);
  fscanf(soilparam, "%lf", &temp.quartz);
  if(options.FULL_ENERGY && (temp.quartz>1. || temp.quartz<0)) {
    fprintf(stderr,"Need valid quartz content as a fraction to run");
    fprintf(stderr," Full Energy model, %f is not acceptable.\n",
	    temp.quartz);
    exit(0);
  }
  for(layer=0;layer<options.Nlayer;layer++)
    fscanf(soilparam, "%lf", &temp.bulk_density[layer]);
  fscanf(soilparam, "%lf", &temp.soil_density);
  for(layer=0;layer<options.Nlayer;layer++)
    if(temp.bulk_density[layer]>=temp.soil_density)
      nrerror("Layer bulk density must be less then soil density");
  fscanf(soilparam, "%lf", &off_gmt);
  fscanf(soilparam, "%lf", &Wcr_FRACT);
  fscanf(soilparam, "%lf", &Wpwp_FRACT);
  fscanf(soilparam, "%lf", &temp.rough);
  fscanf(soilparam, "%lf", &temp.snow_rough);

  /*******************************************
    Compute Maximum Soil Layer Moiture Content
    *******************************************/
  for(layer=0;layer<options.Nlayer;layer++) {
    porosity = 1.0 - temp.bulk_density[layer]/temp.soil_density;
    if(porosity < MOIST_RESID) {
      sprintf(errstr,"Layer %i porosity (%lf mm/mm) must be greater then the defined moisture residue of %lf mm/mm",layer,porosity,(double)MOIST_RESID);
      nrerror(errstr);
    }
    temp.max_moist[layer] = temp.depth[layer] * porosity * 1000.;
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
      temp.Wcr[layer]  = Wcr_FRACT * temp.max_moist[layer];
      temp.Wpwp[layer] = Wpwp_FRACT * temp.max_moist[layer];
      if(temp.Wpwp[layer] > temp.Wcr[layer])
        nrerror("Wpwp is greater then Wcr");
  }

  /*************************************************
    Determine Central Longitude of Current Time Zone 
    *************************************************/
  temp.time_zone_lng = off_gmt * 360./24.;

  return temp;
} 
