#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void prepare_full_energy(int               iveg,
			 int               Nveg,
			 int               Nnodes,
			 all_vars_struct  *all_vars,
			 soil_con_struct  *soil_con,
			 double           *moist0,
			 double           *ice0) {
/*******************************************************************
  prepare_full_energy.c      Keith Cherkauer       January 20, 2000

  This subroutine returns the soil thermal properties, moisture 
  and ice contents for the top two layers for use with the QUICK_FLUX
  ground heat flux solution.

  Modifications:
  01-20-00 split into separate file, formerly at the end of 
           full_energy.c                                      KAC
  03-12-03 modified so that ice content is set to zero unless
           the frozen soil algorithm is implemented and active
           in the current grid cell.                          KAC
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2008-Jan-23 Changed ice0 from a scalar to an array.  Previously,
	      when options.SNOW_BAND > 1, the value of ice0 computed
	      for earlier bands was always overwritten by the value
	      of ice0 computed for the final band (even if the final
	      band had 0 area).						JS via TJB
  2008-May-05 Changed moist from a scalar to an array (moist0).  Previously,
	      when options.SNOW_BAND > 1, the value of moist computed
	      for earlier bands was always overwritten by the value
	      of moist computed for the final band (even if the final
	      band had 0 area).						KAC via TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
*******************************************************************/

  extern option_struct options;

  int                i, band;
  double            *null_ptr;
  layer_data_struct *layer;

  layer = (layer_data_struct *)calloc(options.Nlayer,
				      sizeof(layer_data_struct));

  for(band=0;band<options.SNOW_BAND;band++) {

    if (soil_con->AreaFract[band] > 0.0) {

      for(i=0;i<options.Nlayer;i++) 
        layer[i] = all_vars->cell[iveg][band].layer[i];
    
      /* Compute top soil layer moisture content (mm/mm) */

      moist0[band] = layer[0].moist / ( soil_con->depth[0] * 1000. );

      /* Compute top soil layer ice content (mm/mm) */

      if(options.FROZEN_SOIL && soil_con->FS_ACTIVE){
        if((all_vars->energy[iveg][band].T[0] 
	    + all_vars->energy[iveg][band].T[1])/2.<0.) {
	  ice0[band] = moist0[band] 
	    - maximum_unfrozen_water((all_vars->energy[iveg][band].T[0]
				      + all_vars->energy[iveg][band].T[1]) / 2.,
				     soil_con->max_moist[0]
				     / (soil_con->depth[0] * 1000.),
				     soil_con->bubble[0], soil_con->expt[0]);
	  if(ice0[band]<0.) ice0[band]=0.;
        }
        else ice0[band]=0.;
      }
      else {
        ice0[band] = 0.;
      }

      /** Compute Soil Thermal Properties **/
      compute_soil_layer_thermal_properties(layer,soil_con->depth,
					    soil_con->bulk_dens_min,
					    soil_con->soil_dens_min,
					    soil_con->quartz,
					    soil_con->bulk_density,
					    soil_con->soil_density,
					    soil_con->organic,
					    soil_con->frost_fract,
					    options.Nlayer);
    
      /** Save Thermal Conductivities for Energy Balance **/
      all_vars->energy[iveg][band].kappa[0] = layer[0].kappa; 
      all_vars->energy[iveg][band].Cs[0]    = layer[0].Cs; 
      all_vars->energy[iveg][band].kappa[1] = layer[1].kappa; 
      all_vars->energy[iveg][band].Cs[1]    = layer[1].Cs; 

    }

    else {

      ice0[band] = 0.;

    }

  }

  free((char *)layer);

}
