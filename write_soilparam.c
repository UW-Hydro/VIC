#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_soilparam(soil_con_struct soil)
/**********************************************************************
	write_soilparam		Dag Lohmann	January 1996

  This routine writes soil parameters to stdout.  Used to check that
  the correct parameters are in fact being read into the model.

  Modifications:
  5/21/96	Routine rewritten to account for variable number
		of layers					KAC
  4-12-98  Modified to output all standard aoil parameters for the
           VIC-NL model                                         KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  int i;
  int startlayer;

  if(options.FROZEN_SOIL) startlayer=2;
  else startlayer=0;

  printf("Soil Parameters\n");
  printf("\tLat: %f     Lon: %f\n",soil.lat,soil.lng);
  printf("\tbi                   = %lf [Infiltration parameter]\n",soil.b_infilt);
  printf("\tDs                   = %lf [Subsurface flow rate]\n",soil.Ds);
  printf("\tDsmax                = %lf mm/day [Maximum subsurface flow rate]\n",soil.Dsmax);
  printf("\tWs                   = %lf [Soil Water Content]\n",soil.Ws);
  printf("\tc                    = %lf\n",soil.c);
  for(i=0;i<options.Nlayer;i++)
    printf("\tExpt%02d             = %lf [exponential]\n",i+1,soil.expt[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tKsat%02d             = %lf mm/day [Saturated hydraulic conductivity]\n",i+1,soil.Ksat[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tPhi_s%02d            = %lf mm/mm [soil moisture diffusion coefficient]\n",i+1,soil.phi_s[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tinit_moist%02d       = %lf mm [Initial soil layer moisture]\n",i+1,soil.init_moist[i]);
  printf("\televation            = %lf m [Average elevation]\n",soil.elevation);
  for(i=0;i<options.Nlayer;i++)
    printf("\tdepth%02d            = %lf m [Soil layer thickness]\n",i+1,soil.depth[i]);
  printf("\tavg_temp             = %lf C [Average soil temperature]\n",soil.avg_temp);
  printf("\tdp                   = %lf m [Soil thermal damping depth]\n",soil.dp);
  printf("\tbubble               = %lf cm [Bubbling Pressure]\n",soil.bubble);
  printf("\tquartz               = %lf fract [Quartz content]\n",soil.quartz);
  for(i=0;i<options.Nlayer;i++)
    printf("\tbulk_density%02d     = %lf kg/m^3 [Bulk density]\n",i+1,soil.bulk_density[i]);
  printf("\tsoil_density         = %lf kg/m^3 [Soil partical density]\n",soil.soil_density);
  printf("\ttime_zone_lng     = %f degrees [Central longitude of time zone]\n",soil.time_zone_lng);
  for(i=0;i<options.Nlayer;i++)
    printf("\tmax_moist%02d     = %lf mm [Maximum moisture content]\n",soil.max_moist[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tWcr%02d           = %lf mm [Critical moisture content]\n",soil.Wcr[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tWpwp%02d          = %lf mm [Wilting point moisture content]\n",soil.Wpwp[i]);
  printf("\trough             = %lf m [Roughness of bare soil]\n",soil.rough);
  printf("\tsnow_rough     = %lf m [Roughness of snow surface]\n",soil.snow_rough);
}




