#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void initialize_snow (snow_data_struct *snow, 
		      int veg_num, FILE *fsnow)
/**********************************************************************
	initialize_snow		Keith Cherkauer		January 22, 1997

  This routine initializes the snow variable arrays for each new
  grid cell.

  VARIABLES INITIALIZED:
    snow[i].pack_water 		Liquid water content of snow pack (m) 
    snow[i].surf_water 		Liquid water content of surface layer (m) 
    snow[i].vapor_flux 		Mass flux of water vapor to or from the
					intercepted snow (m) 
    snow[i].pack_temp 		Temperature of snow pack (C) 
    snow[i].melt_energy 	Energy used for melting and heating of snow
					pack 
    snow[i].snow 		1 = snow present, 0 = no snow 
    snow[i].last_snow 		days since last snowfall 
    snow[i].swq 		Snow water equivalent at current pixel (m) 
    snow[i].surf_temp 		Temperature of snow pack surface layer (C) 
    snow[i].density 		Snow pack density (kg/m^3) 
    snow[i].depth		Snow pack depth, computed from density 
				and swq (m).

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char tempstr[512];
  int i;
  int startlayer;

  if(options.FROZEN_SOIL) startlayer=2;
  else startlayer=0;

  for ( i = 0 ; i <= veg_num ; i++) {
    snow[i].pack_water = 0.0;	/* Liquid water content of snow pack (m) */
    snow[i].surf_water = 0.0;	/* Liquid water content of surface layer (m) */
    snow[i].vapor_flux = 0.0;	/* Mass flux of water vapor to or from the
                               intercepted snow (m) */
    snow[i].pack_temp = 0.0;	/* Temperature of snow pack (C) */
    snow[i].melt_energy = 0.0;	/* Energy used for melting and heating of snow
                                   pack */
    snow[i].snow_canopy = 0.0;	/* Snow on canopy */
    snow[i].tmp_int_storage = 0.0;
    if(options.INIT_SNOW) {
      rewind(fsnow);
      fgets(tempstr,512,fsnow);
      fscanf(fsnow,"%*s %c\n",&snow[i].snow);
      fscanf(fsnow,"%*s %i\n",&snow[i].last_snow);
      fscanf(fsnow,"%*s %lf\n",&snow[i].swq);
      fscanf(fsnow,"%*s %lf\n",&snow[i].surf_temp);
      fscanf(fsnow,"%*s %lf\n",&snow[i].density);
    }
    else {	/* Assume no snow present */
      snow[i].snow = 0;		/* 1 = snow present, 0 = no snow */
      snow[i].last_snow = 0;	/* days since last snowfall */
      snow[i].swq = 0.0;	/* Snow water equivalent at current pixel (m) */
      snow[i].surf_temp = 0.0;	/* Temperature of snow pack surface layer (C) */
      snow[i].density = 0.0;	/* Snow pack density (kg/m^3) */
    }
    if(snow[i].density>0.) snow[i].depth = 1000. * snow[i].swq
                                         / snow[i].density;
    else snow[i].depth = 0.;
  }

  if(debug.DEBUG || debug.PRT_SNOW) {
    fprintf(debug.fg_snow,"Date\tSWE TOT\tSWE SRF\tSWE PCK\tGRND T\tlyr1 T\tSURF T\tPACK T\tMELT\tVPR FLX\tAIR T\tSNOW\tRAIN\tGRNDFLX\tDEPTH\tKAPPA\tCANOPY\tCNPYFLUX\n");
  }

}
