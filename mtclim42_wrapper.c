/*
 * Purpose: Initialize and call mtclim42 routines to estimate meteorological
 *          variables 
 * Usage  :
 * Author : Bart Nijssen
 * E-mail : nijssen@u.washington.edu
 * Created:
 * Last Changed: Fri Apr 25 11:51:24 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * Notes  : Acts as interface between mtclim42 and vic
 * Disclaimer: Feel free to use or adapt any part of this program for your own
 *             convenience.  However, this only applies with the understanding
 *             that YOU ARE RESPONSIBLE TO ENSURE THAT THE PROGRAM DOES WHAT YOU
 *             THINK IT SHOULD DO.  The author of this program does not in any
 *             way, shape or form accept responsibility for problems caused by
 *             the use of this code.  At any time, please feel free to discard
 *             this code and WRITE YOUR OWN, it's what I would do.
 */

  // Modified 04-25-03 to checnge calls to vicerror into calls to nrerror.  
  //          Vicerror may cause errors as all variables are not yet 
  //          defined.                                                KAC
  // Modified 04-25-03 to compute pressures in Pa, rather than kPa.   KAC

/******************************************************************************/
/*			    PREPROCESSOR DIRECTIVES                           */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>
#include <mtclim42_vic.h>

/******************************************************************************/
/*			TYPE DEFINITIONS, GLOBALS, ETC.                       */
/******************************************************************************/

/******************************************************************************/
/*			      FUNCTION PROTOTYPES                             */
/******************************************************************************/
void mtclim42_init(int have_dewpt, double elevation, double annual_prcp, 
		   double lat, global_param_struct *vic_global, dmy_struct *dmy, 
		   double *prec, double *tmax, double *tmin, double *hourlyrad, 
		   double *tiny_radfract, control_struct *ctrl, 
		   parameter_struct *p, data_struct *mtclim42_data); 

void mtclim42_to_vic(int have_dewpt, int have_shortwave, double hour_offset, 
		     global_param_struct *vic_global, dmy_struct *dmy, 
		     double *tiny_radfract, control_struct *ctrl, 
		     data_struct *mtclim42_data, double *tskc, double *vp, 
		     double *hourlyrad);

void mtclim42_wrapper(int have_dewpt, int have_shortwave, double hour_offset,
		      double elevation, double annual_prcp, double lat, 
		      global_param_struct *vic_global, dmy_struct *dmy, 
		      double *prec, double *tmax, double *tmin, double *tskc,
		      double *vp, double *hourlyrad) 
{
  control_struct ctrl;
  parameter_struct p;
  data_struct mtclim42_data;
  long ntinys;
  double *tiny_radfract;

  /* allocate space for the tiny_radfract array */
  ntinys = (long) 86400L/(long)SRADDT;
  ntinys *= 366L;
  tiny_radfract = (double *) calloc(ntinys, sizeof(double));
  if (tiny_radfract == NULL) {
    nrerror("Memory allocation error in mtclim42_init() ...\n");
  }

  /* initialize the mtclim 4.2 data structures */ 
  mtclim42_init(have_dewpt, elevation, annual_prcp, lat, vic_global, dmy, prec,
		tmax, tmin, hourlyrad, tiny_radfract, &ctrl, &p,
		&mtclim42_data);  

  /* calculate daily air temperatures */
  if (calc_tair(&ctrl, &p, &mtclim42_data)) {
    nrerror("Error in calc_tair()... exiting\n");
  }
  
  /* calculate daily precipitation */
  if (calc_prcp(&ctrl, &p, &mtclim42_data)) {
    nrerror("Error in calc_prcp()... exiting\n");
  }
  
  /* test for the presence of Tdew observations, and branch to the
     appropriate srad and humidity algorithms */
  if (ctrl.indewpt) {
    /* calculate srad and humidity using real Tdew data */
    if (calc_srad_humidity(&ctrl, &p, &mtclim42_data, tiny_radfract)) {
      nrerror("Error in calc_srad_humidity()... exiting\n");
    }
  }
  else { /* no dewpoint temperature data */
    /* calculate srad and humidity with iterative algorithm */
    if (calc_srad_humidity_iterative(&ctrl, &p, &mtclim42_data,
				     tiny_radfract)) { 
      nrerror("Error in calc_srad_humidity_iterative()... exiting\n");
    }
  }

  /* translate the mtclim 4.2 structures back to the VIC data structures */
  mtclim42_to_vic(have_dewpt, have_shortwave, hour_offset, vic_global,
		  dmy, tiny_radfract, &ctrl,&mtclim42_data, tskc, vp,
		  hourlyrad);

  /* clean up */
  if (data_free(&ctrl, &mtclim42_data)) {
    nrerror("Error in data_free()... exiting\n");
  }  
  free(tiny_radfract);
}
  
void mtclim42_init(int have_dewpt, double elevation, double annual_prcp, 
		   double lat, global_param_struct *vic_global, dmy_struct *dmy, 
		   double *prec, double *tmax, double *tmin, double *hourlyrad, 
		   double *tiny_radfract, control_struct *ctrl, 
		   parameter_struct *p, data_struct *mtclim42_data)
{
  int i;
  int stepspday;
  long tinystepspyear;

  /* initialize the control structure */

  ctrl->ndays = (vic_global->nrecs*vic_global->dt)/24;
  
  if (have_dewpt)
    nrerror("have_dewpt not yet implemented...\n");
  else
    ctrl->indewpt = 0;
  ctrl->outhum = 1;		/* output vapor pressure */
  ctrl->inyear = 0;
  
  /* initialize the parameter structure.  Meteorological variables are only
     calculated for the mean grid cell elevation.  The temperatures are lapsed
     outside of the mtclim code.  Therefore p->base_elev and p->site_elev are
     set to the same value.  The same is true for p->base_isoh and
     p->site_isoh. */
  p->base_elev   = elevation;
  p->base_isoh   = annual_prcp/10.; /* MTCLIM prcp in cm */
  p->site_lat    = lat;
  p->site_elev   = elevation;
  p->site_slp    = 0.;
  p->site_asp    = 0.;
  p->site_isoh   = annual_prcp/10.; /* MTCLIM prcp in cm */
  p->site_ehoriz = 0.;
  p->site_whoriz = 0.;
  p->tmax_lr     = T_lapse;	    /* not used since site_elev == base_elev */
  p->tmin_lr     = T_lapse;	    /* not used since site_elev == base_elev */

  /* allocate space in the data arrays for input and output data */
  if (data_alloc(ctrl, mtclim42_data)) {
    nrerror("Error in data_alloc()... exiting\n");
  }

  /* First populate the solar day array with the tmin, tmax, prcp and possibly
     the Tdew values for the corresponding local days.  At
     this point we will not worry about the first and last incomplete solar days 
     (if there are any).  It could be argued that the mtclim method does not
     make all that much sense if you only have data for one or two days
     anyway. */

  /* in this first version we just take care of the daily data, subdaily data
     will be implemented later */
     
  /* initialize the data arrays with the vic input data */
  stepspday = 24/vic_global->dt;
  for (i = 0; i < ctrl->ndays; i++) {
    mtclim42_data->yday[i] = dmy[i*stepspday].day_in_year;
    mtclim42_data->tmax[i] = tmax[i];
    mtclim42_data->tmin[i] = tmin[i];
    /* MTCLIM prcp in cm */
    mtclim42_data->prcp[i] = prec[i]/10.; 
    if (have_dewpt)
      nrerror("have_dewpt not yet implemented ...\n");
  }
  tinystepspyear = 366L*(86400L/(long)SRADDT);
  for (i = 0; i < tinystepspyear; i++)
    tiny_radfract[i] = 0;
}

void mtclim42_to_vic(int have_dewpt, int have_shortwave, double hour_offset, 
		     global_param_struct *vic_global, dmy_struct *dmy, 
		     double *tiny_radfract, control_struct *ctrl, 
		     data_struct *mtclim42_data, double *tskc, double *vp, 
		     double *hourlyrad)
{
  int i;
  int j;
  long tinystep;
  long tinystepsphour;
  long tinystepspyear;
  long nhours;
  double dailyrad;
  
  /* First, fill the hourlyrad structure with the amount of actual radiation
     during that hour.  The actual amount of radiation during the hour is found
     by taking the sum of the actual radiation for all SRADDT intervals during
     that hour.  The actual radiation during an SRADDT interval is the product
     of the daily radiation and the fraction of the potential radiation occurring
     during the SRADDT interval.  (SRADDT is defined in mtclim42_vic.h). 

     The time in the hourlyrad structure coincides with the timesteps on which
     vic runs, i.e. the first hour of the hourlyrad array coincides with the
     first hour of the VIC run and the last hour coincides with the last hour of
     the vic run.
     
     Note: If you run on a timestep dt greater than 1 hour, the first dt
     elements of hourlyrad correspond to the first model timestep, the
     second set of dt elements to the second model timestep, etc. 
     
     Note: The tiny_radfract array has values for one complete year, so we can
     "wrap around" if tinystep < 0 or tinystep > tinystepsperyear 

     Note: If this seems confusing and/or overkill, ponder the issue while
     running in vic time with the radiation in solar time and input in GMT */
  
  tinystepsphour = 3600L/(long)SRADDT;
  tinystepspyear = 366L*24L*tinystepsphour;
  nhours = vic_global->nrecs * vic_global->dt;

  tinystep = ((dmy[0].day_in_year-1)*24 - hour_offset)*tinystepsphour;
  if (tinystep < 0)
    tinystep += tinystepspyear;
  for (i = 0; i < nhours; i++) {
    dailyrad = mtclim42_data->s_srad[i/24] * mtclim42_data->s_dayl[i/24]/3600.;
    hourlyrad[i] = 0;
    for (j = 0; j < tinystepsphour; j++, tinystep++) {
      if (tinystep >= tinystepspyear)
	tinystep -= tinystepspyear;
      hourlyrad[i] += tiny_radfract[tinystep];
    }
    hourlyrad[i] *= dailyrad;
  }
  
  for (i = 0; i < ctrl->ndays; i++) {
    tskc[i] = mtclim42_data->s_tskc[i];
    vp[i] = mtclim42_data->s_hum[i]; /* convert to kPa */
  }
}
