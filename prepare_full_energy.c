#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void prepare_full_energy(int               iveg,
			 int               Nveg,
			 int               Nnodes,
			 dist_prcp_struct *prcp,
			 soil_con_struct  *soil_con,
			 double           *moist,
			 double           *ice0) {
/*******************************************************************
  prepare_full_energy.c      Keith Cherkauer       January 20, 2000

  This subroutine returns the soil thermal properties, moisture 
  and ice contents for the top two layers for use with the QUICK_FLUX
  ground heat flux solution.

  Modifications:
  01-20-00 split into separate file, formerly at the end of 
           full_energy.c                                      KAC

*******************************************************************/

  extern option_struct options;

  int                i, band;
  double            *null_ptr;
  layer_data_struct *layer;

  layer = (layer_data_struct *)calloc(options.Nlayer,
				      sizeof(layer_data_struct));

  for(band=0;band<options.SNOW_BAND;band++) {

    /* Compute average soil moisture values for distributed precipitation */

    for(i=0;i<options.Nlayer;i++) 
      layer[i] = find_average_layer(&(prcp->cell[WET][iveg][band].layer[i]),
				    &(prcp->cell[DRY][iveg][band].layer[i]),
				    soil_con->depth[i], prcp->mu[iveg]);
    
    /* Compute top soil layer moisture content (mm/mm) */

    (*moist) = layer[0].moist / ( soil_con->depth[0] * 1000. );

    /* Compute top soil layer ice content (mm/mm) */

    if(options.FROZEN_SOIL){
      if((prcp->energy[iveg][band].T[0] 
	  + prcp->energy[iveg][band].T[1])/2.<0.) {
	(*ice0) = (*moist) 
	  - maximum_unfrozen_water((prcp->energy[iveg][band].T[0]
				    + prcp->energy[iveg][band].T[1]) / 2.,
				   soil_con->max_moist[0]
				   / (soil_con->depth[0] * 1000.),
				   soil_con->bubble[0], soil_con->expt[0]);
	if((*ice0)<0.) (*ice0)=0.;
      }
      else (*ice0)=0.;
    }
    else {
      (*ice0) = 0.;
    }

    /** Compute Soil Thermal Properties **/
    compute_soil_layer_thermal_properties(layer,soil_con->depth,
					  soil_con->bulk_density,
					  soil_con->soil_density,
					  soil_con->quartz,
#if SPATIAL_FROST
					  soil_con->frost_fract,
#endif
					  options.Nlayer);
    
    /** Save Thermal Conductivities for Energy Balance **/
    prcp->energy[iveg][band].kappa[0] = layer[0].kappa; 
    prcp->energy[iveg][band].Cs[0]    = layer[0].Cs; 
    prcp->energy[iveg][band].kappa[1] = layer[1].kappa; 
    prcp->energy[iveg][band].Cs[1]    = layer[1].Cs; 
  }

  free((char *)layer);

}
