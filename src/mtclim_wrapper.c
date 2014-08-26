/*******************************************************************************
 * Purpose: Initialize and call mtclim routines to estimate meteorological
 *          variables 
 * Usage  :
 * Author : Bart Nijssen
 * E-mail : nijssen@u.washington.edu
 * Created:
 * Last Changed: Fri Apr 25 11:51:24 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * Notes  : Acts as interface between mtclim and vic
 * Disclaimer: Feel free to use or adapt any part of this program for your own
 *             convenience.  However, this only applies with the understanding
 *             that YOU ARE RESPONSIBLE TO ENSURE THAT THE PROGRAM DOES WHAT YOU
 *             THINK IT SHOULD DO.  The author of this program does not in any
 *             way, shape or form accept responsibility for problems caused by
 *             the use of this code.  At any time, please feel free to discard
 *             this code and WRITE YOUR OWN, it's what I would do.
  Modifications:
  2003-Apr-25 Change calls to vicerror into calls to nrerror.  
              Vicerror may cause errors as all variables are not yet 
	      defined.							KAC
  2003-Apr-25 Compute pressures in Pa, rather than kPa.			KAC
  2011-Nov-04 Updated to mtclim 4.3. 					TJB
  2011-Nov-04 Modified to handle user-supplied observed SW and/or VP.	TJB
  2011-Nov-04 Fixed bug in tinyradfract array indexing and transfer to
	      hourlyrad array.  Previously the code looped over all 366
	      days' worth of the tinyradfract array, causing the diurnal
	      cycle of radiation to be shifted 1 day later in the year
	      every non-leap year.  At high latitudes this resulted in
	      substantial errors in the diurnal cycle after 20-30 years
	      of simulation.  This has been fixed.			TJB
  2013-Jul-19 Fixed bug in shortwave computation for case when daily shortwave
	      is supplied by the user.					HFC via TJB
  2013-Jul-25 Added data->s_fdir.					TJB

*******************************************************************************/
/******************************************************************************/
/*			    PREPROCESSOR DIRECTIVES                           */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>
#include <mtclim_constants_vic.h>
#include <mtclim_parameters_vic.h>

static char vcid[] = "$Id$";

/******************************************************************************/
/*			TYPE DEFINITIONS, GLOBALS, ETC.                       */
/******************************************************************************/

/******************************************************************************/
/*			      FUNCTION PROTOTYPES                             */
/******************************************************************************/
void mtclim_init(int have_dewpt, int have_shortwave, double elevation, double slope, double aspect,
                   double ehoriz, double whoriz, double annual_prcp, 
		   double lat, int Ndays, dmy_struct *dmy, 
		   double *prec, double *tmax, double *tmin, double *vp, double *hourlyrad, 
		   double **tiny_radfract, control_struct *ctrl, 
		   parameter_struct *p, data_struct *mtclim_data); 

void mtclim_to_vic(double hour_offset, 
		     int Ndays, dmy_struct *dmy, 
		     double **tiny_radfract, control_struct *ctrl, 
		     data_struct *mtclim_data, double *tskc, double *vp, 
		     double *hourlyrad, double *fdir);

void mtclim_wrapper(int have_dewpt, int have_shortwave, double hour_offset,
		      double elevation, double slope, double aspect,
                      double ehoriz, double whoriz,
                      double annual_prcp, double lat, 
		      int Ndays, dmy_struct *dmy, 
		      double *prec, double *tmax, double *tmin, double *tskc,
		      double *vp, double *hourlyrad, double *fdir) 
/******************************************************************************
  mtclim_wrapper: interface between VIC and MTCLIM.

  Modifications:
  2012-Feb-16 Cleaned up commented code.					TJB
******************************************************************************/
{
  control_struct ctrl;
  parameter_struct p;
  data_struct mtclim_data;
  double **tiny_radfract;
  int i;

  /* allocate space for the tiny_radfract array */
  tiny_radfract = (double **) calloc(366, sizeof(double*));
  if (tiny_radfract == NULL) {
    nrerror("Memory allocation error in mtclim_init() ...\n");
  }
  for (i=0; i<366; i++) {
    tiny_radfract[i] = (double *) calloc(86400, sizeof(double));
    if (tiny_radfract[i] == NULL) {
      nrerror("Memory allocation error in mtclim_init() ...\n");
    }
  }

  /* initialize the mtclim data structures */ 
  mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz, whoriz,
                annual_prcp, lat, Ndays, dmy, prec,
		tmax, tmin, vp, hourlyrad, tiny_radfract, &ctrl, &p,
		&mtclim_data);  

  /* calculate daily air temperatures */
  if (calc_tair(&ctrl, &p, &mtclim_data)) {
    nrerror("Error in calc_tair()... exiting\n");
  }
  
  /* calculate daily precipitation */
  if (calc_prcp(&ctrl, &p, &mtclim_data)) {
    nrerror("Error in calc_prcp()... exiting\n");
  }
  
  /* calculate daily snowpack using simple model (this is only for radiation correction, *not* the same as the VIC snowpack estimate) */
  if (snowpack(&ctrl, &p, &mtclim_data)) {
    nrerror("Error in snowpack()... exiting\n");
  }
  
  /* calculate srad and humidity with iterative algorithm */
  if (calc_srad_humidity_iterative(&ctrl, &p, &mtclim_data, tiny_radfract)) { 
    nrerror("Error in calc_srad_humidity_iterative()... exiting\n");
  }

  /* translate the mtclim structures back to the VIC data structures */
  mtclim_to_vic(hour_offset, Ndays,
		  dmy, tiny_radfract, &ctrl,&mtclim_data, tskc, vp,
		  hourlyrad, fdir);

  /* clean up */
  if (data_free(&ctrl, &mtclim_data)) {
    nrerror("Error in data_free()... exiting\n");
  }
  for (i=0; i<366; i++) {
    free(tiny_radfract[i]);
  }
  free(tiny_radfract);
}
  
void mtclim_init(int have_dewpt, int have_shortwave, double elevation, double slope, double aspect,
                   double ehoriz, double whoriz, double annual_prcp, 
		   double lat, int Ndays, dmy_struct *dmy, 
		   double *prec, double *tmax, double *tmin, double *vp, double *hourlyrad, 
		   double **tiny_radfract, control_struct *ctrl, 
		   parameter_struct *p, data_struct *mtclim_data)
{
  int i,j;
  int tinystepspday;

  /* initialize the control structure */

  ctrl->ndays = Ndays;
  
  ctrl->indewpt = 0;
  ctrl->invp = 0;
  if (have_dewpt) {
    if (have_dewpt == 1) {
      nrerror("have_dewpt not yet implemented for tdew; however you can supply observed vapor pressure and set have_dewpt to 2\n");
    }
    else if (have_dewpt == 2) {
      ctrl->invp = 1;
    }
  }
  if (have_shortwave) ctrl->insw = 1;
  else ctrl->insw = 0;
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
  p->site_slp    = slope;
  p->site_asp    = aspect;
  p->site_isoh   = annual_prcp/10.; /* MTCLIM prcp in cm */
  p->site_ehoriz = ehoriz;
  p->site_whoriz = whoriz;
  p->tmax_lr     = -1*T_LAPSE*METERS_PER_KM;	/* not used since site_elev == base_elev */
  p->tmin_lr     = -1*T_LAPSE*METERS_PER_KM;	/* not used since site_elev == base_elev */

  /* allocate space in the data arrays for input and output data */
  if (data_alloc(ctrl, mtclim_data)) {
    nrerror("Error in data_alloc()... exiting\n");
  }

  /* initialize the data arrays with the vic input data */
  for (i = 0; i < ctrl->ndays; i++) {
    mtclim_data->yday[i] = dmy[i*24].day_in_year;
    mtclim_data->tmax[i] = tmax[i];
    mtclim_data->tmin[i] = tmin[i];
    if (ctrl->insw) {
      mtclim_data->s_srad[i] = 0;
      for (j=0; j<24; j++) {
        mtclim_data->s_srad[i] += hourlyrad[i*24+j];
      }
      mtclim_data->s_srad[i] /= 24;
    }
    if (ctrl->invp) mtclim_data->s_hum[i] = vp[i];
    /* MTCLIM prcp in cm */
    mtclim_data->prcp[i] = prec[i]/10.; 
    if (have_dewpt==1)
      nrerror("have_dewpt not yet implemented ...\n");
  }
  tinystepspday = 86400/SRADDT;
  for (i = 0; i < 366; i++) {
    for (j = 0; j < tinystepspday; j++) {
      tiny_radfract[i][j] = 0;
    }
  }
}

void mtclim_to_vic(double hour_offset, 
		     int Ndays, dmy_struct *dmy, 
		     double **tiny_radfract, control_struct *ctrl, 
		     data_struct *mtclim_data, double *tskc, double *vp, 
		     double *hourlyrad, double *fdir)
/******************************************************************************
  mtclim_to_vic: Store MTCLIM variables in VIC arrays.

  Modifications:
  2012-Feb-16 Removed check on mtclim_data->insw for storing tinyradfract data
	      in hourlyrad array.						TJB
******************************************************************************/
{
  int i,j,k;
  int tinystepsphour;
  int tinystep;
  int tiny_offset;
  double tmp_rad;
  
  tinystepsphour = 3600/SRADDT;

  tiny_offset = (int)((float)tinystepsphour * hour_offset);
  for (i = 0; i < ctrl->ndays; i++) {
    // s_srad = avg SW flux (W/m2) over daylight hours
    // s_dayl = number of seconds of daylight in current day
    // total_daily_sw = s_srad*s_dayl  (J/m2)
    // tiny_radfrac = fraction of total daily sw falling in each SRADDT interval
    // hourlyrad = SW flux (W/m2) over each hour = total_daily_sw * sum_over_hour(tiny_radfract) / 3600
    //                                           = tmp_rad * sum_over_hour(tiny_radfract)

    //if radiation read from input file, assume it's a 24 hours average, 
    //else (i.e., MTCLIM calculated), assume it's a daylight period average
    if (ctrl->insw) {
      tmp_rad = mtclim_data->s_srad[i] * 24.;
    }
    else {
      tmp_rad = mtclim_data->s_srad[i] * mtclim_data->s_dayl[i] / 3600.;
    }    
    for (j = 0; j < 24; j++) {
      hourlyrad[i*24+j] = 0;
      for (k = 0; k < tinystepsphour; k++) {
        tinystep = j*tinystepsphour + k - tiny_offset;
        if (tinystep < 0) {
          tinystep += 24*tinystepsphour; 
        }
        if (tinystep > 24*tinystepsphour-1) {
          tinystep -= 24*tinystepsphour; 
        }
        hourlyrad[i*24+j] += tiny_radfract[dmy[i*24+j].day_in_year-1][tinystep];
      }
      hourlyrad[i*24+j] *= tmp_rad;
    }
  }
  
  for (i = 0; i < ctrl->ndays; i++) {
    fdir[i] = mtclim_data->s_fdir[i];
    tskc[i] = mtclim_data->s_tskc[i];
    vp[i] = mtclim_data->s_hum[i];
  }
}
