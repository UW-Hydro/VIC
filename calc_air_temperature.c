#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/****************************************************************************
  Subroutines developed by Bart Nijssen to estimate the daily temperature
  cycle from maximum and minimum daily temperature measurements.  Modified
  June 23, 1998 by Keith Cherkauer to be run within the VIC-NL model.
  ***************************************************************************/

/****************************************************************************/
/*				    hermite                                 */
/****************************************************************************/
/* calculate the coefficients for the Hermite polynomials */
void hermite(int n, 
	     double *x, 
	     double *yc1, 
	     double *yc2, 
	     double *yc3, 
	     double *yc4)
{
  int i;
  double dx;
  double divdf1;
  double divdf3;
  
  for (i = 0; i < n-1; i++) {
    dx = x[i+1] - x[i];
    divdf1 = (yc1[i+1] - yc1[i])/dx;
    divdf3 = yc2[i] + yc2[i+1] - 2 * divdf1;
    yc3[i] = (divdf1 - yc2[i] - divdf3)/dx;
    yc4[i] = divdf3/(dx*dx);
  }
}

/**************************************************************************/
/*				    hermint                               */
/**************************************************************************/
/* use the Hermite polynomials, to find the interpolation function value at 
   xbar */
double hermint(double xbar, int n, double *x, double *yc1, double *yc2, 
	       double *yc3, double *yc4)
{
  int klo,khi,k;
  double dx;
  double result;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (x[k] > xbar) khi=k;
    else klo=k;
  }

  dx = xbar - x[klo];
  result = yc1[klo] + dx * (yc2[klo] + dx * (yc3[klo] + dx * yc4[klo]));
  return result;
}

/****************************************************************************/
/*				    HourlyT                                 */
/****************************************************************************/
double HourlyT(int Dt, int *TmaxHour, double *Tmax, 
	       int *TminHour, double *Tmin, double *Tair)
{
  double *x;
  double *Tyc1;
  double *yc2;
  double *yc3;
  double *yc4;
  double new_Tmin;
  int i;
  int j;
  int n;
  int hour;
  int nHours;

  nHours = HOURSPERDAY;
  n     = 6;
  x     = (double *) calloc(n, sizeof(double));
  Tyc1  = (double *) calloc(n, sizeof(double));
  yc2   = (double *) calloc(n, sizeof(double));
  yc3   = (double *) calloc(n, sizeof(double));
  yc4   = (double *) calloc(n, sizeof(double));

  /* First fill the x vector with the times for Tmin and Tmax, and fill the 
     Tyc1 with the corresponding temperature and humidity values */
  for (i = 0, j = 0, hour = 0.5 * Dt; i < 3; i++, hour += HOURSPERDAY) {
    if (TminHour[i] < TmaxHour[i]) {
      x[j]       = TminHour[i] + hour;
      Tyc1[j++]  = Tmin[i];
      x[j]       = TmaxHour[i] + hour;
      Tyc1[j++]  = Tmax[i];
    }
    else {
      x[j]       = TmaxHour[i] + hour;
      Tyc1[j++]  = Tmax[i];
      x[j]       = TminHour[i] + hour;
      Tyc1[j++]  = Tmin[i];
    }

  }

  /* we want to preserve maxima and minima, so we require that the first 
     derivative at these points is zero */
  for (i = 0; i < n; i++)
    yc2[i] = 0.;

  /* calculate the coefficients for the splines for the temperature */
  hermite(n, x, Tyc1, yc2, yc3, yc4);

  /* interpolate the temperatures */
  new_Tmin=100;
  for (i = 0, hour = 0.5*Dt+HOURSPERDAY; i < nHours; i++, hour += Dt) {
    Tair[i] = hermint(hour, n, x, Tyc1, yc2, yc3, yc4);
    if(Tair[i]<new_Tmin) new_Tmin=Tair[i];
  }

  free((char*)x);
  free((char*)Tyc1);
  free((char*)yc2);
  free((char*)yc3);
  free((char*)yc4);

  return(new_Tmin);

}

void store_max_min_temp(atmos_data_struct *atmos,
                        double            *tmax,
			int               *tmax_hour,
                        double            *tmin,
			int               *tmin_hour,
                        int                rec,
                        int                Nrecs,
			int                skip_recs) {
/**********************************************************************
  This subroutine sets prepares arrays of maximum and minimum 
  daily air temperature, and the hour of day at which the maximum
  and minimum temperatures are achieved.  These arrays are used 
  when estimating the daily temperature cycle from tmax and tmin.

  atmos_data_struct *atmos,        atmospheric forcing data structure
  double            *tmax,         maximum temperatures for three days
  int               *tmax_hour,    hour of maximum temperature
  double            *tmin,         minimum temperatures for three days
  int               *tmin_hour,    hour of minimum temperature
  int                rec,          current record number
  int                Nrecs,        total number of records
  int                skip_recs     number of subdaily records (24/dt)

**********************************************************************/

  static int last_rec;

  if(rec==0) {
    tmax[0] = tmax[1] = atmos[0].tmax;
    tmin[0] = tmin[1] = atmos[0].tmin;
    tmax[2] = atmos[skip_recs].tmax;
    tmin[2] = atmos[skip_recs].tmin;
    tmax_hour[0] = (int)(2. / 3. * (atmos[0].set_hour 
				    - atmos[0].rise_hour)) 
      + atmos[0].rise_hour;
    tmin_hour[0] = atmos[0].rise_hour - 1;
    tmax_hour[1] = tmax_hour[0];
    tmin_hour[1] = tmin_hour[0];
    tmax_hour[2] = (int)(2. / 3. * (atmos[skip_recs].set_hour 
				    - atmos[skip_recs].rise_hour)) 
      + atmos[0].rise_hour;
    tmin_hour[2] = atmos[skip_recs].rise_hour - 1;
    
    last_rec = rec;
  }
  else if(rec>=last_rec+skip_recs) {
    last_rec = rec;
    if(rec>=Nrecs-skip_recs) {
      tmax[0] = tmax[1];
      tmax[1] = tmax[2];
      tmin[0] = tmin[1];
      tmin[1] = tmin[2];
      tmax_hour[0] = tmax_hour[1];
      tmax_hour[1] = tmax_hour[2];
      tmin_hour[0] = tmin_hour[1];
      tmin_hour[1] = tmin_hour[2];
    }
    else {
      tmax[0] = tmax[1];
      tmax[1] = tmax[2];
      tmax[2] = atmos[skip_recs].tmax;
      tmin[0] = tmin[1];
      tmin[1] = tmin[2];
      tmin[2] = atmos[skip_recs].tmin;
      tmax_hour[0] = tmax_hour[1];
      tmax_hour[1] = tmax_hour[2];
      tmax_hour[2] = (int)(2. / 3. * (atmos[skip_recs].set_hour 
				      - atmos[skip_recs].rise_hour)) 
	+ atmos[skip_recs].rise_hour;
      tmin_hour[0] = tmin_hour[1];
      tmin_hour[1] = tmin_hour[2];
      tmin_hour[2] = atmos[skip_recs].rise_hour - 1;
    }
  }
}

