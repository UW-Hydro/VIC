#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_data(out_data_struct *out_data,
		outfiles_struct outfiles,
		dmy_struct      *dmy)

/**********************************************************************
	write_data	Dag Lohmann		Janurary 1996

  This subroutine writes out the final results for soil moisture, 
  evaporation, runoff and baseflow.

  OUTPUT:
	evaporation in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in mm

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  int j;

  if(options.FROZEN_SOIL) {
    fprintf(outfiles.fdepth  ,"%04i\t%02i\t%02i\t%02i\t%.5lf\t%.5lf",
	    dmy->year, dmy->month, dmy->day, dmy->hour, out_data->fdepth[0], 
            out_data->fdepth[1]);
    for(j=0;j<options.Nlayer;j++) {
      fprintf(outfiles.fdepth,"\t%lf", out_data->ice[j] + out_data->moist[j]);
    }
    fprintf(outfiles.fdepth,"\n"); 
  }

  if(options.FULL_ENERGY && options.SNOW_MODEL) {
    fprintf(outfiles.snow  ,"%04i\t%02i\t%02i\t%02i\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n",
	    dmy->year, dmy->month, dmy->day, dmy->hour,
	    out_data->swq, out_data->snow_depth, out_data->snow_canopy,
            out_data->advection, out_data->coldcontent, out_data->melt_energy);
  }
  else if(!options.FULL_ENERGY && options.SNOW_MODEL) {
    fprintf(outfiles.snow  ,"%04i\t%02i\t%02i\t%02i\t%.5lf\t%.5lf\t%.5lf\n",
	    dmy->year, dmy->month, dmy->day, dmy->hour,
	    out_data->swq, out_data->snow_depth, out_data->snow_canopy);
  }


  if(options.FULL_ENERGY) {
    fprintf(outfiles.fluxes,"%04i\t%02i\t%02i\t%02i\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t\n",
	    dmy->year, dmy->month, dmy->day, dmy->hour,
	    out_data->prec, out_data->evap, out_data->runoff,
	    out_data->baseflow, out_data->Wdew, out_data->moist[0],
	    out_data->moist[1], out_data->moist[2], out_data->rad_temp,
	    out_data->net_short, out_data->r_net, out_data->latent,
	    out_data->latent_canop, out_data->latent_trans,
	    out_data->latent_bare, out_data->latent_pet, 
	    out_data->latent_pet_mm, out_data->sensible,
	    out_data->grnd_flux, out_data->aero_resist, out_data->surf_cond, 
	    out_data->albedo);
  }
  else {
    fprintf(outfiles.fluxes,"%04i\t%02i\t%02i\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t\n",
	    dmy->year, dmy->month, dmy->day,
	    out_data->prec, out_data->evap, out_data->runoff, 
	    out_data->baseflow, out_data->Wdew, out_data->moist[0], 
	    out_data->moist[1], out_data->moist[2],
	    out_data->net_short, out_data->r_net, out_data->latent_canop, 
	    out_data->latent_trans, out_data->latent_bare, 
	    out_data->latent_pet, out_data->latent_pet_mm, 
	    out_data->aero_resist, out_data->surf_cond, out_data->albedo);
  }
}

void calc_water_balance_error(int    rec,
			      double inflow,
			      double outflow,
			      double storage) {
/***************************************************************
  This subroutine computes the overall model water balance
***************************************************************/

  static double last_storage;
  static double cum_error;

  double error;

  if(rec<0) {
    last_storage = storage;
    cum_error = 0.;
  }
  else {
    error = inflow - outflow - (storage - last_storage);
    cum_error += error;
    last_storage = storage;
    if(fabs(error)>1.e-5)
      fprintf(stderr,"Moist Error:\t%i\t%.5lf\t%.5lf\n",rec,error,cum_error);
  }

}

void calc_energy_balance_error(int    rec,
			       double net_rad,
			       double latent,
			       double sensible,
			       double grnd_flux) {
/***************************************************************
  This subroutine computes the overall model energy balance
***************************************************************/

  static double cum_error;

  double error;

  if(rec<0) cum_error=0;
  else {
    error = net_rad + latent + sensible + grnd_flux;
    cum_error += error;
    if(fabs(error)>1.e-5)
      fprintf(stderr,"Energy Error:\t%i\t%.5lf\t%.5lf\n",rec,error,cum_error);
  }
}
