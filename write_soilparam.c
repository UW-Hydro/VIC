#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_soilparam.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void write_soilparam(soil_con_struct *soil)
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

  int i;

  printf("Soil Parameters\n");
  printf("\tLat: %f     Lon: %f\n",soil->lat,soil->lng);
  printf("\tbi                   = %f [Infiltration parameter]\n",
	 soil->b_infilt);
  printf("\tDs                   = %f [Subsurface flow rate]\n",
	 soil->Ds);
  printf("\tDsmax                = %f mm/day [Maximum subsurface flow rate]\n",
	 soil->Dsmax);
  printf("\tWs                   = %f [Soil Water Content]\n",soil->Ws);
  printf("\tc                    = %f\n",soil->c);
  for(i=0;i<options.Nlayer;i++)
    printf("\tExpt%02d             = %f [exponential]\n",i+1,soil->expt[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tKsat%02d             = %f mm/day [Saturated hydraulic conductivity]\n",i+1,soil->Ksat[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tPhi_s%02d            = %f mm/mm [soil moisture diffusion coefficient]\n",i+1,soil->phi_s[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tinit_moist%02d       = %f mm [Initial soil layer moisture]\n",
	   i+1,soil->init_moist[i]);
  printf("\televation            = %f m [Average elevation]\n",
	 soil->elevation);
  for(i=0;i<options.Nlayer;i++)
    printf("\tdepth%02d            = %f m [Soil layer thickness]\n",
	   i+1,soil->depth[i]);
  printf("\tavg_temp             = %f C [Average soil temperature]\n",
	 soil->avg_temp);
  printf("\tdp                   = %f m [Soil thermal damping depth]\n",
	 soil->dp);
  for(i=0;i<options.Nlayer;i++)
    printf("\tbubble%02d           = %f cm [Bubbling Pressure]\n",i,soil->bubble[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tquartz%02d           = %f fract [Quartz content]\n",i,soil->quartz[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tbulk_density%02d     = %f kg/m^3 [Bulk density]\n",
	   i+1,soil->bulk_density[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tsoil_density%02d     = %f kg/m^3 [Soil partical density]\n",
	   i+1,soil->soil_density[i]);
  printf("\ttime_zone_lng     = %f degrees [Central longitude of time zone]\n",soil->time_zone_lng);
  for(i=0;i<options.Nlayer;i++)
    printf("\tmax_moist%02d     = %f mm [Maximum moisture content]\n",
	   i+1,soil->max_moist[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tWcr%02d           = %f mm [Critical moisture content]\n",
	   i+1,soil->Wcr[i]);
  for(i=0;i<options.Nlayer;i++)
    printf("\tWpwp%02d          = %f mm [Wilting point moisture content]\n",
	   i+1,soil->Wpwp[i]);
  printf("\trough             = %f m [Roughness of bare soil]\n",soil->rough);
  printf("\tsnow_rough     = %f m [Roughness of snow surface]\n",
	 soil->snow_rough);
}




