#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_soilparam(soil_con_struct soil)
/**********************************************************************
	write_soilparam		Dag Lohmann	January 1996

  This routine writes soil parameters to the screen.

  Modifications:
  5/21/96	Routine rewritten to account for variable number
		of layers					KAC

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
  printf("\tDs                   = %lf [Subsurface flow rate ()]\n",soil.Ds);
  printf("\tDsmax                = %lf [Maximum subsurface flow rate ()]\n",soil.Dsmax);
  printf("\tWs                   = %lf [Soil Water Content (cm^3/cm^3)]\n",soil.Ws);
  printf("\tc                    = %lf\n",soil.c);
  if(options.FROZEN_SOIL || options.FULL_ENERGY)
    printf("\tbubble               = %lf [Bubbling Pressure (cm)]\n",soil.bubble);
  printf("\tSoil Moisture Top    = %lf [mm]\n",soil.init_moist[0]);
  printf("\tSoil Moisture Bottom = %lf [mm]\n",soil.init_moist[1]);
  printf("\televation            = %f [m]\n",soil.elevation);
  for(i=0;i<options.Nlayer;i++)
    printf("\tmax_moist%d           = %lf [mm]\n",i+1,soil.max_moist[i]);
  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    for(i=0;i<options.Nlayer;i++){
      printf("\tdepth%d               = %lf [m]\n",i+1,soil.depth[i]);
      printf("\texponent%d            = %lf\n",i+1, soil.expt[i]);
      printf("\tKsat%d                = %lf [Saturated Hydrologic Conductivity (mm/day)]\n",i+1,soil.Ksat[i]);
      printf("\tphi_s%d               = %lf [Saturated Matrix Potential (mm)]\n",i+1,soil.phi_s[i]);
      printf("\tbulk_density%d        = %lf [Bulk Density (kg/m^3)]\n",i+1,soil.bulk_density[i]);
    }
  }
  if(options.FROZEN_SOIL || options.FULL_ENERGY) {
    printf("\tsoil_density         = %lf [Partical Density (kg/m^3)]\n",soil.soil_density);
    printf("\tTime Zone Central Meridian = %f [degrees]\n",soil.time_zone_lng);
  }
}




