#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

soil_con_struct read_soilparam(FILE *soilparam)
/**********************************************************************
	read_soilparam		Dag Lohmann		January 1996

  This routine reads soil parameters for each grid cell.

  Modifications:
  7-19-96	Modified to read through variable layers, and
		read soil depth and average temperature for
		full energy and frozen soil versions of the
		model.						KAC

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
  fscanf(soilparam, "%lf", &temp.dp);
  fscanf(soilparam, "%lf", &temp.bubble);
  fscanf(soilparam, "%lf", &temp.quartz);
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

  for(layer=0;layer<options.Nlayer;layer++) {
    porosity = 1.0 - temp.bulk_density[layer]/temp.soil_density;
    if(porosity < MOIST_RESID) {
      sprintf(errstr,"Layer %i porosity (%lf mm/mm) must be greater then the defined moisture residue of %lf mm/mm",layer,porosity,(double)MOIST_RESID);
      nrerror(errstr);
    }
    temp.max_moist[layer] = temp.depth[layer] * porosity * 1000.;
  }

  if(options.Nlayer==2)
    temp.max_infil = (1.0+temp.b_infilt)*temp.max_moist[0];
  else
    temp.max_infil = (1.0+temp.b_infilt)*(temp.max_moist[0]+temp.max_moist[1]);

  for(layer=0;layer<options.Nlayer;layer++) {
      temp.Wcr[layer]  = Wcr_FRACT * temp.max_moist[layer];
      temp.Wpwp[layer] = Wpwp_FRACT * temp.max_moist[layer];
      if(temp.Wpwp[layer] > temp.Wcr[layer])
        nrerror("Wpwp is greater then Wcr");
  }

  temp.time_zone_lng = off_gmt * 360./24.;

  return temp;
} 
