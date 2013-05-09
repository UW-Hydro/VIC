#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_data.c,v 2.8 1998/10/02 01:34:50 vicadmin Exp $";

void write_data(out_data_struct *out_data,
		outfiles_struct  outfiles,
		dmy_struct      *dmy,
		int              dt)
/**********************************************************************
	write_data	Dag Lohmann		Janurary 1996

  This subroutine writes all energy and moisture balance parameters to
  output files.

  OUTPUT:
	evaporation and vapor fluxes in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in cm
	snow depth in cm
	snow water equivlence in mm
	all energy fluxes are in W/m^2

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC
  3/98          Routine modified to output fluxes in PILPS2c 
                ASCII column format                             Dag
  4/30/98       Routine modified to add binary output options for
                improved file speed, and less disk usage for large
		model basins                                    KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  int   band, j;
  int   *tmp_iptr;
  float *tmp_fptr;

/*   tmp_iptr = (int *)calloc(1,sizeof(int)); */
/*   tmp_fptr = (float *)calloc(1,sizeof(float)); */

  /************************************
    Output Frozen Soil Variables
  ************************************/
/*   if(options.FROZEN_SOIL && options.BINARY_OUTPUT) { */
    /***** Write Binary Frozen Soil Output File *****/
/*     tmp_iptr[0] = dmy->year; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fdepth); */
/*     tmp_iptr[0] = dmy->month; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fdepth); */
/*     tmp_iptr[0] = dmy->day; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fdepth); */
/*     tmp_iptr[0] = dmy->hour; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fdepth); */
/*     tmp_fptr[0] = (float)out_data->fdepth[0]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fdepth); */
/*     tmp_fptr[0] = (float)out_data->fdepth[1]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fdepth); */
/*     for(j=0;j<options.Nlayer;j++) { */
/*       tmp_fptr[0] = (float)(out_data->ice[j] + out_data->moist[j]); */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fdepth); */
/*     } */
/*   } */
/*   else if(options.FROZEN_SOIL) { */
    /***** Write ASCII Frozen Soil Output File *****/
/*     fprintf(outfiles.fdepth  ,"%04i\t%02i\t%02i\t%02i\t%.4lf\t%.4lf", */
/* 	    dmy->year, dmy->month, dmy->day, dmy->hour, out_data->fdepth[0],  */
/*             out_data->fdepth[1]); */
/*     for(j=0;j<options.Nlayer;j++) { */
/*       fprintf(outfiles.fdepth,"\t%lf", out_data->ice[j] + out_data->moist[j]); */
/*     } */
/*     fprintf(outfiles.fdepth,"\n");  */
/*   } */

  /**************************************
    Ouput Snow Variables
  **************************************/
/*   if(options.SNOW_MODEL && options.BINARY_OUTPUT) { */
    /***** Write Binary Snow Output File *****/
/*     tmp_iptr[0] = dmy->year; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snow); */
/*     tmp_iptr[0] = dmy->month; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snow); */
/*     tmp_iptr[0] = dmy->day; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snow); */
/*     if(options.FULL_ENERGY) { */
/*       tmp_iptr[0] = dmy->hour; */
/*       fwrite(tmp_iptr,1,sizeof(int),outfiles.snow); */
/*     } */
/*     tmp_fptr[0] = (float)out_data->swq[0]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*     tmp_fptr[0] = (float)out_data->snow_depth[0]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*     tmp_fptr[0] = (float)out_data->snow_canopy[0]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*     tmp_fptr[0] = (float)out_data->coverage[0]; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*     if(options.FULL_ENERGY) { */
/*       tmp_fptr[0] = (float)out_data->advection[0]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*       tmp_fptr[0] = (float)out_data->deltaCC[0]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*       tmp_fptr[0] = (float)out_data->snow_flux[0]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*       tmp_fptr[0] = (float)out_data->refreeze_energy[0]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snow); */
/*     } */
/*   } */
/*   else if(options.FULL_ENERGY && options.SNOW_MODEL) { */
    /***** Write ASCII full energy snow output file *****/
/*     fprintf(outfiles.snow  ,"%04i\t%02i\t%02i\t%02i\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n", */
/* 	    dmy->year, dmy->month, dmy->day, dmy->hour, */
/* 	    out_data->swq[0], out_data->snow_depth[0],  */
/* 	    out_data->snow_canopy[0], */ /* out_data->coverage[0], */ 
/* 	    out_data->advection[0], out_data->deltaCC[0],  */
/* 	    out_data->snow_flux[0], out_data->refreeze_energy[0]); */
/*   } */
/*   else if(!options.FULL_ENERGY && options.SNOW_MODEL) { */
    /***** Write ASCII water balance snow output file *****/
/*     fprintf(outfiles.snow  ,"%04i\t%02i\t%02i\t%02i\t%.4lf\t%.4lf\t%.4lf\n", */
/* 	    dmy->year, dmy->month, dmy->day, dmy->hour, */
/* 	    out_data->swq[0], out_data->snow_depth[0],  */
/* 	    out_data->snow_canopy[0]); */
/*   } */

  /**************************************
    Ouput Snow Band Variables
  **************************************/
/*   if(options.SNOW_MODEL && options.PRT_SNOW_BAND && options.BINARY_OUTPUT) { */
    /***** Write Binary Snow Band Output File *****/
/*     tmp_iptr[0] = dmy->year; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snowband); */
/*     tmp_iptr[0] = dmy->month; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snowband); */
/*     tmp_iptr[0] = dmy->day; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.snowband); */
/*     if(options.FULL_ENERGY) { */
/*       tmp_iptr[0] = dmy->hour; */
/*       fwrite(tmp_iptr,1,sizeof(int),outfiles.snowband); */
/*     } */
/*     for(band=0;band<options.SNOW_BAND;band++) { */
/*       tmp_fptr[0] = (float)out_data->swq[band+1]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/*       tmp_fptr[0] = (float)out_data->snow_depth[band+1]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/*       tmp_fptr[0] = (float)out_data->snow_canopy[band+1]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
      /*     tmp_fptr[0] = (float)out_data->coverage; */
      /*     fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/*       if(options.FULL_ENERGY) { */
/* 	tmp_fptr[0] = (float)out_data->advection[band+1]; */
/* 	fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/* 	tmp_fptr[0] = (float)out_data->deltaCC[band+1]; */
/* 	fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/* 	tmp_fptr[0] = (float)out_data->snow_flux[band+1]; */
/* 	fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/* 	tmp_fptr[0] = (float)out_data->refreeze_energy[band+1]; */
/* 	fwrite(tmp_fptr,1,sizeof(float),outfiles.snowband); */
/*       } */
/*     } */
/*   } */
/*   else if(options.FULL_ENERGY && options.SNOW_MODEL && options.PRT_SNOW_BAND) { */
    /***** Write ASCII full energy snow band output file *****/
/*     fprintf(outfiles.snowband  ,"%04i\t%02i\t%02i\t%02i", */
/* 	    dmy->year, dmy->month, dmy->day, dmy->hour); */
/*     for(band=0;band<options.SNOW_BAND;band++)  */
/*       fprintf(outfiles.snowband  ,"\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf", */
/* 	      out_data->swq[band+1], out_data->snow_depth[band+1],  */
/* 	      out_data->snow_canopy[band+1], */ /* out_data->coverage[band+1], */ 
/* 	      out_data->advection[band+1], out_data->deltaCC[band+1],  */
/* 	      out_data->snow_flux[band+1], out_data->refreeze_energy[band+1]); */
/*     fprintf(outfiles.snowband,"\n"); */
/*   } */
/*   else if(!options.FULL_ENERGY && options.SNOW_MODEL  */
/* 	  && options.PRT_SNOW_BAND) { */
    /***** Write ASCII water balance snow band output file *****/
/*     fprintf(outfiles.snowband  ,"%04i\t%02i\t%02i\t%02i", */
/* 	    dmy->year, dmy->month, dmy->day, dmy->hour); */
/*     for(band=0;band<options.SNOW_BAND;band++) */
/*       fprintf(outfiles.snowband  ,"\t%.4lf\t%.4lf\t%.4lf", */
/* 	      out_data->swq[band+1], out_data->snow_depth[band+1],  */
/* 	      out_data->snow_canopy[band+1]); */
/*     fprintf(outfiles.snowband,"\n"); */
/*   } */

  /************************************
    Output Standard Energy and Moisture Flux Variables
  ************************************/
/*   if(options.BINARY_OUTPUT) { */
    /***** Write Binary Fluxes File *****/
/*     tmp_iptr[0] = dmy->year; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fluxes); */
/*     tmp_iptr[0] = dmy->month; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fluxes); */
/*     tmp_iptr[0] = dmy->day; */
/*     fwrite(tmp_iptr,1,sizeof(int),outfiles.fluxes); */
/*     if(dt<24) { */
/*       tmp_iptr[0] = dmy->hour; */
/*       fwrite(tmp_iptr,1,sizeof(int),outfiles.fluxes); */
/*     } */
/*     tmp_fptr[0] = (float)out_data->prec; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->evap; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->runoff; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->baseflow; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->Wdew; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     for(j=0;j<options.Nlayer;j++) {     */
/*       tmp_fptr[0] = (float)out_data->moist[j]; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     } */
/*     if(options.FULL_ENERGY) { */
/*       tmp_fptr[0] = (float)out_data->rad_temp; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     } */
/*     tmp_fptr[0] = (float)out_data->net_short; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->r_net; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     if(options.FULL_ENERGY) { */
/*       tmp_fptr[0] = (float)out_data->latent; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     } */
/*     tmp_fptr[0] = (float)out_data->evap_canop; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->evap_veg; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->evap_bare; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->sub_canop; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->sub_snow; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     if(options.FULL_ENERGY) { */
/*       tmp_fptr[0] = (float)out_data->sensible; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*       tmp_fptr[0] = (float)out_data->grnd_flux; */
/*       fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     } */
/*     tmp_fptr[0] = (float)out_data->aero_resist; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->surf_temp; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*     tmp_fptr[0] = (float)out_data->albedo; */
/*     fwrite(tmp_fptr,1,sizeof(float),outfiles.fluxes); */
/*   } */
/*   else if(options.FULL_ENERGY) { */
    /***** Write ASCII energy balance fluxes file *****/
    fprintf(outfiles.fluxes,"%04i\t%02i\t%02i\t%.4le\t%.4le\n",
	    dmy->year, dmy->month, dmy->day, out_data->runoff,
	    out_data->baseflow);
/*     for(j=0;j<options.Nlayer;j++)  */
/*       fprintf(outfiles.fluxes,"\t%.4lf", out_data->moist[j]); */
/*     fprintf(outfiles.fluxes,"\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t\n", */
/* 	    out_data->rad_temp, */
/* 	    out_data->net_short, out_data->r_net, out_data->latent, */
/* 	    out_data->evap_canop, out_data->evap_veg, */
/* 	    out_data->evap_bare, out_data->sub_canop,  */
/* 	    out_data->sub_snow, out_data->sensible, */
/* 	    out_data->grnd_flux, out_data->aero_resist, out_data->surf_temp,  */
/* 	    out_data->albedo); */
/*   } */
/*   else { */
    /***** Write ASCII Water Balance Fluxes File *****/
/*     fprintf(outfiles.fluxes,"%04i\t%02i\t%02i", */
/* 	    dmy->year, dmy->month, dmy->day); */
/*     if(dt<24) */
/*       fprintf(outfiles.fluxes,"\t%02i", dmy->hour); */
/*     fprintf(outfiles.fluxes,"\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf", */
/* 	    out_data->prec, out_data->evap, out_data->runoff,  */
/* 	    out_data->baseflow, out_data->Wdew); */
/*     for(j=0;j<options.Nlayer;j++)  */
/*       fprintf(outfiles.fluxes,"\t%.4lf", out_data->moist[j]); */
/*     fprintf(outfiles.fluxes,"\t%.4lf\t%.4lf\t%.4lg\t%.4lg\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t\n", */
/* 	    out_data->net_short, out_data->r_net, out_data->evap_canop,  */
/* 	    out_data->evap_veg, out_data->evap_bare,  */
/* 	    out_data->sub_canop, out_data->sub_snow,  */
/* 	    out_data->aero_resist, out_data->surf_temp, out_data->albedo); */
/*   } */

/*   free((char *)tmp_iptr); */
/*   free((char *)tmp_fptr); */

}

void calc_water_balance_error(int    rec,
			      double inflow,
			      double outflow,
			      double storage) {
/***************************************************************
  calc_water_balance_error  Keith Cherkauer        April 1998

  This subroutine computes the overall model water balance, and 
  warns the model user if large errors are found.
***************************************************************/

  static double last_storage;
  static double cum_error;
  static double max_error;
  static int    error_cnt;
  static int    Nrecs;

  double error;

  if(rec<0) {
    last_storage = storage;
    cum_error    = 0.;
    max_error    = 0.;
    error_cnt    = 0;
    Nrecs        = -rec;
  }
  else {
    error = inflow - outflow - (storage - last_storage);
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>1e-5) {
      max_error = error;
      fprintf(stderr,"Maximum Moist Error:\t%i\t%.5lf\t%.5lf\n",
	      rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Water Error for Grid Cell = %.4lf\n",
	      cum_error);
    }
    last_storage = storage;
  }

}

void calc_energy_balance_error(int    rec,
			       double net_rad,
			       double latent,
			       double sensible,
			       double grnd_flux,
			       double snow_fluxes) {
/***************************************************************
  calc_energy_balance_error   Keith Cherkauer     April 1998

  This subroutine computes the overall model energy balance, and
  reports the maximum time step error above a thresehold to the
  user.  The total cumulative error for the grid cell is also 
  computed and reported at the end of the model run.
***************************************************************/

  static double cum_error;
  static double max_error;
  static int    Nrecs;

  double error;

  if(rec<0) {
    cum_error = 0;
    Nrecs     = -rec;
    max_error = 0;
  }
  else {
    error = net_rad - latent - sensible - grnd_flux + snow_fluxes;
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>0.0001) {
      max_error = error;
      fprintf(stderr,"Maximum Energy Error:\t%i\t%.4lf\t%.4lf\n",
		rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Energy Error for Grid Cell = %.4lf\n",
	      cum_error);
    }
  }
}
