#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id: calc_root_fraction.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void calc_root_fractions(veg_con_struct  *veg_con,
			 soil_con_struct  *soil_con)
/**********************************************************************
  calc_root_fraction.c    Keith Cherkauer      September 24, 1998

  This routine computes the fraction of roots in each soil layer based
  on the root zone distribution defined in the vegetation parameter
  file.  Roots are assumed to be linearly distributed within each
  root zone.

**********************************************************************/
{
  extern option_struct options;

  char   ErrStr[MAXSTRING];
  int    Nveg;
  int    veg;
  int    layer;
  int    zone;
  int    i;
  float  sum_depth;
  float  sum_fract;
  float  dum;
  double Zstep;
  double Zsum;
  double Lstep;
  double Lsum;
  double Zmin_fract;
  double Zmin_depth;
  double Zmax;

  Nveg      = veg_con[0].vegetat_type_num;

  for(veg=0;veg<Nveg;veg++) {
    sum_depth  = 0;
    sum_fract  = 0;
    layer      = 0;
    Lstep      = soil_con->depth[layer];
    Lsum       = Lstep;
    Zsum       = 0;
    zone       = 0;
    
    while(zone<options.ROOT_ZONES) {
      Zstep = (double)veg_con[veg].zone_depth[zone];
      if((Zsum + Zstep) <= Lsum && Zsum >= Lsum - Lstep) {
	/** CASE 1: Root Zone Completely in Soil Layer **/
	sum_fract += veg_con[veg].zone_fract[zone];
      }
      else {
	/** CASE 2: Root Zone Partially in Soil Layer **/
	if(Zsum < Lsum - Lstep) {
	  /** Root zone starts in previous soil layer **/
	  Zmin_depth = Lsum - Lstep;
	  Zmin_fract = linear_interp(Zmin_depth,Zsum,Zsum+Zstep,0,
				     veg_con[veg].zone_fract[zone]);
	}
	else {
	  /** Root zone starts in current soil layer **/
	  Zmin_depth = Zsum;
	  Zmin_fract = 0.;
	}
	if(Zsum + Zstep <= Lsum) {
	  /** Root zone ends in current layer **/
	  Zmax = Zsum + Zstep;
	}
	else {
	  /** Root zone extends beyond bottom of current layer **/
	  Zmax = Lsum;
	}
	sum_fract += linear_interp(Zmax,Zsum,Zsum+Zstep,0,
				   veg_con[veg].zone_fract[zone]) - Zmin_fract;
      }

      /** Update Root Zone and Soil Layer **/
      if(Zsum + Zstep < Lsum) {
	Zsum += Zstep;
	zone ++;
      }
      else if(Zsum + Zstep == Lsum) {
	Zsum += Zstep;
	zone ++;
	if(layer<options.Nlayer) {
	  veg_con[veg].root[layer] = sum_fract;
	  sum_fract = 0.;
	}
	layer++;
	if(layer<options.Nlayer) {
	  Lstep  = soil_con->depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==options.Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<options.ROOT_ZONES-1) {
	    for(i=zone+1;i<options.ROOT_ZONES;i++) {
	      Lstep += veg_con[veg].zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
      else if(Zsum + Zstep > Lsum) {
	if(layer<options.Nlayer) {
	  veg_con[veg].root[layer] = sum_fract;
	  sum_fract = 0.;
	}
	layer++;
	if(layer<options.Nlayer) {
	  Lstep  = soil_con->depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==options.Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<options.ROOT_ZONES-1) {
	    for(i=zone+1;i<options.ROOT_ZONES;i++) {
	      Lstep += veg_con[veg].zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
	
    }

    if(sum_fract > 0 && layer >= options.Nlayer) {
      veg_con[veg].root[options.Nlayer-1] += sum_fract;
    }
    else if(sum_fract > 0) {
      veg_con[veg].root[layer] += sum_fract;
    }

    dum=0.;
    for (layer=0;layer<options.Nlayer;layer++) {
      if(veg_con[veg].root[layer] < 1.e-4) veg_con[veg].root[layer] = 0.;
      dum += veg_con[veg].root[layer];
    }
    if(dum == 0.0){
      sprintf(ErrStr,"Root fractions sum equals zero: %f , Vege Class: %d\n",
	      dum, veg_con[veg].veg_class);
      nrerror(ErrStr);
    }
    for (layer=0;layer<options.Nlayer;layer++) {
      veg_con[veg].root[layer] /= dum;
    }

  }

}

