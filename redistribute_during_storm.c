#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void redistribute_during_storm(cell_data_struct ***cell,
			       veg_var_struct   ***veg_var,
			       int                 veg,
			       int                 Nveg,
			       int                 rec,
			       double              old_mu,
			       double              new_mu) {
/**********************************************************************
  redistribute_during_storm.c     Keith Cherkauer     January 13, 1998

  This subroutine redistributes soil moisture when the current storm 
  changes intensity.

  Modified:
  06-24-98 Changed to run on only one vegetation type at a time    KAC
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
  for(layer=0;layer<options.Nlayer;layer++) {

    for(band=0;band<options.SNOW_BAND;band++) {
      temp_wet = cell[WET][veg][band].layer[layer].moist_thaw;
      temp_dry = cell[DRY][veg][band].layer[layer].moist_thaw;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, old_mu, 
					      new_mu);
      if(error) {
	sprintf(ErrorString,"Moist_thaw does not balance after change in storm: %lf -> %lf record %i\n",
		cell[WET][veg][band].layer[layer].moist_thaw*new_mu
		+ cell[DRY][veg][band].layer[layer].moist_thaw*(1.-new_mu),
		temp_wet+temp_dry,rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].moist_thaw = temp_wet;
      cell[DRY][veg][band].layer[layer].moist_thaw = temp_dry;
      
      temp_wet = cell[WET][veg][band].layer[layer].moist_froz;
      temp_dry = cell[DRY][veg][band].layer[layer].moist_froz;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"moist_froz does not balance after change in storm: %lf -> %lf record %i\n",
		cell[WET][veg][band].layer[layer].moist_froz*new_mu
		+ cell[DRY][veg][band].layer[layer].moist_froz*(1.-new_mu),
		temp_wet+temp_dry,rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].moist_froz = temp_wet;
      cell[DRY][veg][band].layer[layer].moist_froz = temp_dry;
      
      temp_wet = cell[WET][veg][band].layer[layer].moist;
      temp_dry = cell[DRY][veg][band].layer[layer].moist;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"moist does not balance after change in storm: %lf -> %lf record %i\n",
		cell[WET][veg][band].layer[layer].moist*new_mu
		+ cell[DRY][veg][band].layer[layer].moist*(1.-new_mu),
		temp_wet+temp_dry,rec);
	vicerror(ErrorString);
      }
      cell[WET][veg][band].layer[layer].moist = temp_wet;
      cell[DRY][veg][band].layer[layer].moist = temp_dry;
      
      temp_wet = cell[WET][veg][band].layer[layer].ice;
      temp_dry = cell[DRY][veg][band].layer[layer].ice;
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"ice does not balance after change in storm: %lf -> %lf record %i\n",
		cell[WET][veg][band].layer[layer].ice*new_mu
		+ cell[DRY][veg][band].layer[layer].ice*(1.-new_mu),
		temp_wet+temp_dry,rec);
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
      error = redistribute_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
	sprintf(ErrorString,"Wdew does not balance after change in storm: %lf -> %lf record %i\n",
		veg_var[WET][veg][band].Wdew*new_mu
		+ veg_var[DRY][veg][band].Wdew*(1.-new_mu),
		temp_wet+temp_dry,rec);
	vicerror(ErrorString);
      }
      veg_var[WET][veg][band].Wdew = temp_wet;
      veg_var[DRY][veg][band].Wdew = temp_dry;
    }
  }
}

unsigned char redistribute_moisture_for_storm(double *wet_value,
					      double *dry_value,
					      double old_mu,
					      double new_mu) {
/**********************************************************************
  This subroutine redistributes the given parameter between wet and
  dry cell fractions when the precipitation changes in intensity.
**********************************************************************/

  unsigned char error;
  double temp_wet;
  double temp_dry;
  double diff;

  temp_wet = *wet_value * old_mu;
  temp_dry = *dry_value * (1. - old_mu);
  if(old_mu>new_mu && (new_mu!=1. && new_mu!=0.)) {
    *wet_value
        = (((1. - (old_mu-new_mu)/old_mu)*temp_wet) / new_mu);
    *dry_value
        = (((old_mu-new_mu)/old_mu*temp_wet + temp_dry) / (1.-new_mu));
  }
  else if(new_mu!=1. && new_mu!=0.) {
    *wet_value
        = (((new_mu-old_mu)/(1.-old_mu)*temp_dry + temp_wet) / new_mu);
    *dry_value
        = (((1. - (new_mu-old_mu)/(1.-old_mu))*temp_dry) / (1.-new_mu));
  }
  else {
    *wet_value
        = (temp_dry + temp_wet);
    *dry_value
        = (temp_dry + temp_wet);
  }
  diff = (temp_wet+temp_dry) - (*wet_value*new_mu + *dry_value*(1.-new_mu));
  if(fabs(diff) > 1.e-10) error = TRUE;
  else error = FALSE;

  return (error);
}
