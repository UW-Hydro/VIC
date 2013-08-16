#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id: initialize_new_storm.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void initialize_new_storm(cell_data_struct ***cell,
			  veg_var_struct   ***veg_var,
			  int                 veg,
			  int                 Nveg,
			  int                 rec,
			  double              old_mu,
			  double              new_mu) {
/**********************************************************************
  initialize_new_storm.c	Keith Cherkauer		June 24, 1998

  This subroutine averages all soil moisture components before the
  start of a new storm, so that the model loses memory of the previous
  storm's location (which was never really known).

  Modifications:
  07-13-98 modified to redistribute Wdew within all defined
           elevation bands                                         KAC

**********************************************************************/
 
  extern option_struct options;

  unsigned char error;
  char          ErrorString[MAXSTRING];
  int           layer;
  int           band;
  double        temp_wet;
  double        temp_dry;

  /** Redistribute Soil Moisture **/
  for(layer = 0; layer < options.Nlayer; layer++) {

    for(band = 0; band < options.SNOW_BAND; band++) {

      temp_wet = cell[WET][veg][band].layer[layer].moist;
      temp_dry = cell[DRY][veg][band].layer[layer].moist;
      error = average_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"moist does not balance before new storm: %f -> %f record %i\n",
		cell[WET][veg][band].layer[layer].moist * new_mu
		+ cell[DRY][veg][band].layer[layer].moist * (1. - new_mu),
		temp_wet + temp_dry, rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].moist = temp_wet;
      cell[DRY][veg][band].layer[layer].moist = temp_dry;
      
      temp_wet = cell[WET][veg][band].layer[layer].ice;
      temp_dry = cell[DRY][veg][band].layer[layer].ice;
      error = average_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"ice does not balance before new storm: %f -> %f record %i\n",
		cell[WET][veg][band].layer[layer].ice * new_mu
		+ cell[DRY][veg][band].layer[layer].ice * (1. - new_mu),
		temp_wet + temp_dry, rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].ice = temp_wet;
      cell[DRY][veg][band].layer[layer].ice = temp_dry;
      
    }
  }

  /****************************************
    Redistribute Stored Water in Vegetation
  ****************************************/
  if(veg<Nveg) {
    for(band=0;band<options.SNOW_BAND;band++) {
      temp_wet = veg_var[WET][veg][band].Wdew;
      temp_dry = veg_var[DRY][veg][band].Wdew;
      error = average_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"Wdew does not balance before new storm: %f -> %f record %i\n",
		veg_var[WET][veg][band].Wdew * new_mu
		+ veg_var[DRY][veg][band].Wdew * (1. - new_mu),
		temp_wet + temp_dry, rec);
	vicerror(ErrorString);
      }
      veg_var[WET][veg][band].Wdew = temp_wet;
      veg_var[DRY][veg][band].Wdew = temp_dry;
    }
  }
}

unsigned char average_moisture_for_storm(double *wet_value,
					 double *dry_value,
					 double  old_mu,
					 double  new_mu) {
/**********************************************************************
  This subroutine averages total soil moisture between the wet and dry
  fractions of the soil column
**********************************************************************/

  unsigned char error;
  double temp_wet;
  double temp_dry;
  double diff;

  temp_wet = *wet_value * old_mu;
  temp_dry = *dry_value * (1. - old_mu);
  *wet_value = temp_wet + temp_dry;
  *dry_value = temp_wet + temp_dry;
  diff = (temp_wet+temp_dry) - (*wet_value*new_mu + *dry_value*(1.-new_mu));
  if(fabs(diff) > 1.e-10) error = TRUE;
  else error = FALSE;

  return (error);
}
