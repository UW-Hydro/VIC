/*
  Note:  The following code is largely taken from mtclim4.3, a weather
         preprocessor developed by the NTSG group in the School of Forestry at
	 the University of Montana.  The code has been left intact as much as
	 possible, but small adaptations have been made to allow integration
	 with VIC.  Adaptations have been made so the functions, etc. are local
	 to this routine, which should allow easier updating.
	 The original comments have largely been left intact
	 typedefs and function definitions have been moved to mtclim42_vic.h
	 functions are called from mtclim42_wrapper.c
  Changes:  The only changes to the original MTCLIM code are 
            - cloudiness is returned for each day as well through the variable
	    tskc (calculated using equation 2.29 in Bras, R. L., "Hydrology, an
	    introduction to hydrologic science", Addison-Wesley, Reading,
	    Massachusetts, 1990).
	    - For each day of the year the fraction of the total daily potential
	    radiation that occurs during each radiation time step is also
	    returned through the variable tiny_radfract
	    - SRADDT has been changed to be 300 sec
	    Changes are preceded by the comment * start vic_change *  and
	    followed by the comment * end vic_change *
  Author: Most of the code was written by Peter E. Thornton at the Univeristy of 
          Montana (see below).
	  Adapted by: Bart Nijssen
          Adaptation started on Sat Aug 21 using the mtclim4.2 code from 5/7/1999 
  Modifications:
  2004-Jun-04 If data->s_srad becomes negative, we set it to 0.0.		TJB
  2005-Jun-15 If p->base_isoh and p->site_isoh are both small, set ratio = 1.
	      This handles the case in which annual precip for the grid cell
	      is 0, resulting in both p->base_isoh and p->site_isoh being 0,
	      and their ratio being undefined.					TJB
  2007-Apr-21 Initialize ratio to -1.						TJB
  2011-Nov-05 MTCLIM code updated to MTCLIM v. 4.3.				TJB
  2011-Nov-05 Modified tiny_radfract array to fix a bug in VIC's looping over
	      the array outside of MTCLIM.					TJB
  2011-Nov-05 Added SW_PREC_THRESH, the minimum daily precip necessary to
	      trigger dimming of shortwave.  The purpose of this is to handle
	      situations where gridded forcings are used for input, and trace
	      amounts of precip from neighboring cells are "leaking" into the
	      current cell due to the regridding scheme.			TJB
  2012-Feb-16 Removed function calc_srad_humidity(), which is no longer called.
	      Cleaned up unneeded comment statements.				TJB
  2012-Mar-28 Moved computation of data->s_tskc to after conditional overwriting of
	      data->s_tfmax, so that data->s_tskc reflects observed SW if observed
	      SW is provided.							TJB
  2012-Apr-03 Added check to make sure t_tmax is always positive.		TJB
  2012-Apr-13 Simplified the relationship between tskc and tfmax for LW_CLOUD
	      == LW_CLOUD_DEARDORFF.						TJB
  2013-Jul-19 Fixed bug in shortwave computation for case when daily shortwave
	      is supplied by the user.						HFC via TJB
  2013-Jul-25 Added data->s_fdir.						TJB
*/

/*
mtclim43.c

MTCLIM
VERSION 4.3  

Peter Thornton
NTSG, School of Forestry
University of Montana
1/20/00

***************************************
** Questions or comments? Contact... **
** Dr. Peter E. Thornton             **
** NTSG, School of Forestry          **
** University of Montana             **
** Missoula, MT 59812                **
** email: peter@ntsg.umt.edu         **
** phone: 406-243-4326               **
***************************************

--------------------------------------------------------------------------
Code History
------------
Original code written by R.R. Nemani
Updated   4/1/1989 by J.C. Coughlan
Updated 12/23/1989 by Joe Glassy
Updated   1/4/1993 by Raymond Hunt (version 2.1)
Updated  3/26/1997 by Peter Thornton (version 3.0)
Updated  7/28/1997 by Peter Thornton (version 3.1)
Updated   5/7/1998 by Peter Thornton (version 4.0)
Updated   8/1/1998 by Peter Thornton (version 4.1)
Updated  4/20/1999 by Peter Thornton (version 4.2)
Updated  1/20/2000 by Peter Thornton (version 4.3)

CHANGES FROM VERSION 4.2 TO VERSION 4.3
These changes based on the results of:
Thornton, P.E., H. Hasenauer, and M.A. White, 2000. Simultaneous
estimation of daily solar radiation and humidity from observed temperature
and precipitation: an application over complex terrain in Austria.
For submission to Ag For Met.
1) ADDED THE SNOW ACCUMULATION/MELT MODEL FOR CORRECTING RADIATION
FOR INFLUENCE OF SNOWPACK.
2) ADDED ALBEDO CORRECTION TO HORIZON OBSTRUCTION CALCULATIONS.
3) HUMIDITY ALGORITHM NOW DEPENDS ON THE ANNUAL RATION PET/PRCP:
FOR VALUES < 2.5, USING THE TMIN=TDEW ALGORITHM, FOR VALUES > 2.5,
USING THE ARID CORRECTION FROM:
Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for
estimating surface humidity from daily minimum temperature. Ag For Met
85:87-98.
4) MOVED ALL MODEL PARAMETERS INTO THE PARAMETERS.H FILE

CHANGES FROM VERSION 4.1 TO VERSION 4.2
1) PUT THE SLOPE AND ASPECT CORRECTIONS TO RADIATION BACK IN. THEY
HAD BEEN REMOVED DURING THE DEVELOPMENT OF THE NEW RADIATION CODE.
This includes the estimation for diffuse radiation on sloping surfaces
during periods when the sun is above a flat horizon, but not providing
direct illumination to the slope.  Site east and west horizon corrections
have also been reintroduced.

CHANGES FROM VERSION 4.0 TO VERSION 4.1
1) ADDITIONAL REVISION OF RADIATION CODE, FOLLOWING SUBMISSION OF AG FOR MET
MANUSCRIPT DESCRIBING ANALYSIS OF SAMSON DATA.  The current code
follows exactly the algorithm detailed in
Thornton, P.E. and S.W. Running, 1999. An improved algorithm for estimating 
incident daily solar radiation from measurements of temperature, humidity,
and precipitation. Ag For Met 93:211-228.

changes from version 3.1 to version 4.0:
1) radiation code completely revamped, following analysis of samson
observations for the vemap2 project.
2) includes an iterative algorithm for estimating vapor pressure and
radiation
3) par output no longer an option, since the old par algorithm is suspect
4) 2-day tmin smoothing no longer an option
5) solar constant now calculated as an analytical function of yearday,
instead of the monthly array of values in earlier versions. removes the
discontinuities between months. 
6) boxcar function now uses the previous n days of data, making it a
"pulled boxcar" instead of a "centered boxcar"

Changes from version 2.0 to version 3.1
Modified to include the following improvements to the original MTCLIM logic:
1) Improved vapor pressure calculation
2) Improved transmissivity calculation
3) Improved daylength calculation

Some other differences between version 3.1 and previous versions of MTCLIM:
1) No english units option
2) No total or average radiation option
3) No threshold radiation option
5) No pre-rainy days correction for transmissivity
6) No radiation correction to air temperatures
7) No LAI corrections to air temperatures
8) Only one precipitation station allowed
9) Some parameters formerly in initialization file are now in mtclim_const.h

--------------------------------------------------------------------------
UNITS
-----
Temperatures           degrees C
Temp. lapse rates      degrees C / 1000 m
Precipitation          cm / day
Vapor pressure         Pa
Vapor pressure deficit Pa
Radiation              W/m2, average over daylight period
Daylength              s (sunrise to sunset, flat horizons)
Elevation              m
Latitude               decimal degrees
Aspect                 decimal degrees
Slope                  decimal degrees
E/W horizons           decimal degrees

--------------------------------------------------------------------------
Files
-----

Parameter initialization file
Input meteorological data file
Output meteorological data file  (*.mtcout)

Example initialization file:

---top of init file-------------------------------------------------------
sample               (this entire line written to outfiles)
other comments       (this entire line discarded)

IOFILES                    Keyword, don't edit this line 
test.mtcin                 Name of input meteorological data file
test                       Prefix for output file
                      
CONTROL                    Keyword, don't edit this line
3                          (int) Number of header lines in input file
365                        (int) Number of days of data in input file
0                          (int) Dewpoint temperature input? (0=NO, 1=YES)
1                          (int) Humidity output flag        (0=VPD, 1=VP)
1                          (int) Year field in input file?   (0=NO, 1=YES)

PARAMETERS                 Keyword, don't edit this line
500.0                      (double) Base station elevation, meters
50.0                       (double) Base station annual precip isohyet, cm
48.0                       (double) Site latitude, degrees (- for south)
1500.0                     (double) Site elevation, meters
15.0                       (double) Site slope, degrees
180.0                      (double) Site aspect, degrees (0=N,90=E,180=S,270=W)
75.0                       (double) Site annual precip isohyet, cm
2.0                        (double) Site east horizon, degrees
5.0                        (double) Site west horizon, degrees
-6.0                       (double) Lapse rate for max temperature, deg C/1000m
-3.0                       (double) Lapse rate for min temperature, deg C/1000m

END                        Keyword, don't edit this line
---end of init file-------------------------------------------------------

*.init FILE INFO
For all lines, except the header and comment lines, the  parameter value on the
left can be followed by comments on the right, as long as there is whitespace
separating the value on the left from the following comment. Blank lines can be
inserted or deleted, but all keyword lines and parameter lines must be
included.  The order and number of non-blank lines in this file is crucial. The
keywords are there as a rudimentary quality check on the format of the
initialization file, and they ensure that the proper number of lines are read.
They DON'T ensure that the parameters are in the proper order.

NOTE: The dashed lines at the top and bottom of the example file shown above
are NOT part of the file.


INPUT FILE INFO
All input temperatures are in degrees Celcius, and precipitation is in 
centimeters.
Input data file format:
<some number of header lines (can be zero)>
<year (optional)> <yearday> <Tmax> <Tmin> <Tdew (optional)> <precipitation>
.
.
.

Example input data file... without year field, and without dewpoint temperature
---start of input data file --------------
This is a header line, which is discarded
101 16.0 2.0 1.2
102 17.0 3.0  0.1
103 16.5 5.2  0.0
104 20.1 6.4  0.0
---end of input data file ----------------



OUTPUT FILE INFO
The output file is created by appending ".mtcout" to the output filename
prefix specified in the initialization file. Existing files are overwritten,
so be careful to rename files between runs if you want to save the results, 
and for safety, don't use ".mtcout" as the extension for the input data file.

******************************
Input and Output CONTROL FLAGS
******************************
There is a flag in the initialization file that controls the input of dewpoint
temperature information. If your input file does not contain dewpoint data,
set this flag to 0. Otherwise, if you have dewpoint temperature information 
in your input file, set this flag to 1, and be sure that the dewpoint
temperature field is between the tmin and prcp fields, as specified under
the heading "INPUT FILE INFO", above.

There is another flag in the initialization file that controls the output of
humidity information. When set to 0, humidity output is the vapor pressure
deficit from the saturated vapor pressure at the daylight average temperature
to the ambient vapor pressure, in Pa. When the flag is set to 1, the output is
simply the  ambient vapor pressure in Pa.  By using the ambient vapor pressure
(flag  set to 1), you can avoid possible errors related to a difference between
the temperature chosen for saturation calculations (in this case tday) and the
actual temperature at any given time of day.

Another flag controls the treatment of a year-field in the input and output
files. When the flag is set (1), it is assumed that the first field in the
input file contains the year, which can be any integer between -32000 and
32000 (approx). This field is simply copied line-for-line into the output
file, where it preceeds the standard yearday output field.

Example output data file... (using VPD output, no PAR, no year field)
---start of output data file -----------------------------------
MTCLIM v3.1 OUTPUT FILE : Mon Mar 24 14:50:51 1997
MTCLIM version 3.1 testing
  yday    Tmax    Tmin    Tday    prcp      VPD     srad  daylen
       (deg C) (deg C) (deg C)    (cm)     (Pa)  (W m-2)     (s)
   101   10.00   -0.50    7.39    1.80   439.51   379.95   47672
   102   11.00    1.10    7.67    0.15   387.27   364.95   47880
   103   10.50    2.80    6.84    0.00   243.78   302.15   48087
   104   14.10    3.40    9.29    0.00   390.67   390.22   48293
---end of output data file -------------------------------------


*/

/*********************
**                  **
**  START OF CODE   **
**                  **
*********************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vicNl.h>

#include <mtclim_constants_vic.h>   /* physical constants */
#include <mtclim_parameters_vic.h>  /* model parameters */

static char vcid[] = "$Id$";

/****************************
 **                         **
 **    START OF FUNCTION    **
 **       DEFINITIONS       **
 **                         **
 ****************************/

/* data_alloc() allocates space for input and output data arrays */
int data_alloc(const control_struct *ctrl, data_struct *data)
{
  int ok=1;
  int ndays;
  
  ndays = ctrl->ndays;
  
  if (ok && ctrl->inyear && !(data->year = (int*) malloc(ndays * sizeof(int)))) {
    printf("Error allocating for year array\n");
    ok=0;
  } 
  if (ok && !(data->yday = (int*) malloc(ndays * sizeof(int)))) {
    printf("Error allocating for yearday array\n");
    ok=0;
  } 
  if (ok && !(data->tmax = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for tmax array\n");
    ok=0;
  }
  if (ok && !(data->tmin = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for tmin array\n");
    ok=0;
  }
  if (ok && !(data->prcp = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for prcp array\n");
    ok=0;
  }
  if (ok && ctrl->indewpt && 
      !(data->tdew = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for input humidity array\n");
    ok=0;
  }
  if (ok && !(data->s_tmax = (double*) malloc(ndays * sizeof(double))))	{
    printf("Error allocating for site Tmax array\n");
    ok=0;
  }
  if (ok && !(data->s_tmin = (double*) malloc(ndays * sizeof(double))))	{
    printf("Error allocating for site Tmin array\n");
    ok=0;
  }
  if (ok && !(data->s_tday = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site Tday array\n");
    ok=0;
  }
  if (ok && !(data->s_prcp = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site prcp array\n");
    ok=0;
  }
  if (ok && !(data->s_hum = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site VPD array\n");
    ok=0;
  }
  if (ok && !(data->s_srad = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site radiation array\n");
    ok=0;
  }
  if (ok && !(data->s_dayl = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site daylength array\n");
    ok=0;
  }
  if (ok && !(data->s_swe = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site snowpack array\n");
    ok=0;
  }
  /* start vic_change */
  if (ok && !(data->s_fdir = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for direct fraction array\n");
    ok=0;
  }
  if (ok && !(data->s_tskc = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site cloudiness array\n");
    ok=0;
  }
  if (ok && !(data->s_ppratio = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site pet/prcp ratio array\n");
    ok=0;
  }
  if (ok && !(data->s_tfmax = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site pet/prcp ratio array\n");
    ok=0;
  }
  if (ok && !(data->s_ttmax = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for site pet/prcp ratio array\n");
    ok=0;
  }
  /* end vic_change */
  return (!ok);
} /* end of data_alloc() */

/* calc_tair() calculates daily air temperatures */
int calc_tair(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data)
{
  int ok=1;
  int i,ndays;
  double dz;
  double tmean, tmax, tmin;
  
  ndays = ctrl->ndays;
  /* calculate elevation difference in kilometers */
  dz = (p->site_elev - p->base_elev)/1000.0;
  
  /* apply lapse rate corrections to tmax and tmin */
  /* Since tmax lapse rate usually has a larger absolute value than tmin
     lapse rate, it is possible at high elevation sites for these corrections
     to result in tmin > tmax. Check for that occurrence and force
     tmin = corrected tmax - 0.5 deg C. */
  for (i=0 ; i<ndays ; i++) {
    /* lapse rate corrections */
    data->s_tmax[i] = tmax = data->tmax[i] + (dz * p->tmax_lr);
    data->s_tmin[i] = tmin = data->tmin[i] + (dz * p->tmin_lr);
    
    /* derived temperatures */
    tmean = (tmax + tmin)/2.0;
    data->s_tday[i] = ((tmax - tmean)*TDAYCOEF) + tmean;
  }
  
  return (!ok);
}
/* end of calc_tair() */

/* calc_prcp() calculates daily total precipitation */
int calc_prcp(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data)
{
  int ok=1;
  int i,ndays;
  double ratio;
  
  ndays = ctrl->ndays;
  
  /* start vic_change */
  ratio = -1.;
  if ( p->site_isoh < 1e-10 && p->base_isoh < 1e-10 ) {
    /* If base_isoh and site_isoh are both small, set the ratio to 1.
       This handles the case in which annual precip is 0, resulting in
       base_isoh and site_isoh being 0 and their ratio being undefined. */
    ratio = 1.;
  }
  else if (p->base_isoh == 0) {
    vicerror("Error in calc_prcp(): base_isoh == 0 and site_isoh/base_isoh == NaN.");
  }
  else {
    ratio = p->site_isoh / p->base_isoh;
  }
  /* end vic_change */
  
  if (ok) {
    for (i=0 ; i<ndays ; i++) {
      data->s_prcp[i] = data->prcp[i] * ratio;
    }
  }
  
  return (!ok);
}
/* end of calc_prcp() */

/* snowpack() estimates the accumulation and melt of snow for radiation
algorithm corrections */
int snowpack(const control_struct *ctrl, const parameter_struct *p,
	     data_struct *data)
{
  int ok=1;
  int i,ndays,count;
  int start_yday,prev_yday;
  double snowpack,newsnow,snowmelt,sum;
  
  ndays = ctrl->ndays;

  /* first pass to initialize SWE array */
  snowpack = 0.0;
  for (i=0 ; i<ndays ; i++)
  {
    newsnow = 0.0;
    snowmelt = 0.0;
    if (data->s_tmin[i] <= SNOW_TCRIT) newsnow = data->s_prcp[i];
    else snowmelt = SNOW_TRATE * (data->s_tmin[i] - SNOW_TCRIT);
    snowpack += newsnow - snowmelt;
    if (snowpack < 0.0) snowpack = 0.0;
    data->s_swe[i] = snowpack;
  }

  /* use the first pass to set the initial snowpack conditions for the
  first day of data */
  start_yday = data->yday[0];
  if (start_yday == 1) prev_yday = 365;
  else prev_yday = start_yday-1;
  count = 0;
  sum = 0.0;
  for (i=1 ; i<ndays ; i++)
  {
    if (data->yday[i] == start_yday || data->yday[i] == prev_yday)
    {
      count ++;
      sum += data->s_swe[i];
    }
  }
  /* Proceed with correction if there are valid days to reinitialize
  the snowpack estiamtes. Otherwise use the first-pass estimate. */
  if (count)
  {
    snowpack = sum/(double)count;
    for (i=0 ; i<ndays ; i++)
    {
      newsnow = 0.0;
      snowmelt = 0.0;
      if (data->s_tmin[i] <= SNOW_TCRIT) newsnow = data->s_prcp[i];
      else snowmelt = SNOW_TRATE * (data->s_tmin[i] - SNOW_TCRIT);
      snowpack += newsnow - snowmelt;
      if (snowpack < 0.0) snowpack = 0.0;
      data->s_swe[i] = snowpack;
    }
  }

  return (!ok);
}
/* end of snowpack() */

/* iterative estimation of shortwave radiation and humidity */
/* Note: too many changes to maintain the start/end vic change comments */
int calc_srad_humidity_iterative(const control_struct *ctrl,
				 const parameter_struct *p, data_struct *data,
				 double **tiny_radfract)
{
  int ok=1;
  int i,j,ndays;
  int start_yday,end_yday,isloop;
  int ami,yday;
  double ttmax0[366];
  double flat_potrad[366];
  double slope_potrad[366];
  double daylength[366];
  double *dtr, *sm_dtr;
  double *parray, *window, *t_fmax, *tdew;
  double *pet;
  double sum_prcp,ann_prcp,effann_prcp;
  double sum_pet,ann_pet;
  double tmax,tmin;
  double t1,t2;
  double pratio;
  double lat,coslat,sinlat,dt,h,dh;
  double cosslp,sinslp,cosasp,sinasp;
  double bsg1,bsg2,bsg3;
  double decl,cosdecl,sindecl,cosegeom,sinegeom,coshss,hss;
  double sc,dir_beam_topa;
  double sum_flat_potrad,sum_slope_potrad,sum_trans;
  double cosh,sinh;
  double cza,cbsa,coszeh,coszwh;
  double dir_flat_topa,am;
  double t_tmax,b;
  double tmink,ratio,ratio2,ratio3,tdewk;
  double pvs,vpd;
  double trans1,trans2;
  double t_final,pdif,pdir,srad1,srad2; 
  double pa;
  double sky_prop;
  double avg_horizon, slope_excess;
  double horizon_scalar, slope_scalar;
  int update_pva;

  /* optical airmass by degrees */
  double optam[21] = {2.90,3.05,3.21,3.39,3.69,3.82,4.07,4.37,4.72,5.12,5.60,
		      6.18,6.88,7.77,8.90,10.39,12.44,15.36,19.79,26.96,30.00};
  
  /* start vic_change */
  int tinystep;
  int tinystepspday;
  extern option_struct options;
  double tfmax_tmp;
  /* end vic_change */
  
  int iter;
  int max_iter;
  double rmse_tdew;
  double tol;
  double *pva;
  double *tdew_save;
  double *pva_save;

  /* number of simulation days */
  ndays = ctrl->ndays;
  
  /* local array memory allocation */
  /* allocate space for DTR and smoothed DTR arrays */
  if (!(dtr = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for DTR array\n");
    ok=0;
  }
  if (!(sm_dtr = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for smoothed DTR array\n");
    ok=0;
  }
  /* allocate space for effective annual precip array */
  if (!(parray = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for effective annual precip array\n");
    ok=0;
  }
  /* allocate space for the prcp totaling array */
  if (!(window = (double*) malloc((ndays+90)*sizeof(double)))) {
    printf("Error allocating for prcp totaling array\n");
    ok = 0;
  }
  /* allocate space for t_fmax */
  if (!(t_fmax = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for p_tt_max array\n");
      ok=0;
  }
  /* allocate space for Tdew array */
  if (!(tdew = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for Tdew array\n");
    ok=0;
  }
  /* allocate space for pet array */
  if (!(pet = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for pet array\n");
    ok=0;
  }
  /* allocate space for pva array */
  if (!(pva = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for pva array\n");
    ok=0;
  }
  /* allocate space for tdew_save array */
  if (!(tdew_save = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for tdew_save array\n");
    ok=0;
  }
  /* allocate space for pva_save array */
  if (!(pva_save = (double*) malloc(ndays * sizeof(double)))) {
    printf("Error allocating for pva_save array\n");
    ok=0;
  }

  /* calculate diurnal temperature range for transmittance calculations */
  for (i=0 ; i<ndays ; i++) {
    tmax = data->tmax[i];
    tmin = data->tmin[i];
    if (tmax < tmin) 
      tmax = tmin;
    dtr[i] = tmax-tmin;
  }
  
  /* smooth dtr array: After Bristow and Campbell, 1984 */
  if (ndays >= 30) { /* use 30-day antecedent smoothing window */
    if (pulled_boxcar(dtr, sm_dtr, ndays, 30, 0)) {
      printf("Error in boxcar smoothing, calc_srad_humidity()\n");
      ok=0;
    }
  }
  else { /* smoothing window width = ndays */
    if (pulled_boxcar(dtr, sm_dtr, ndays, ndays, 0)) {
      printf("Error in boxcar smoothing, calc_srad_humidity()\n");
      ok=0;
    }
  }

  /* calculate the annual total precip */
  sum_prcp = 0.0;
  for (i=0 ; i<ndays ; i++) {
    sum_prcp += data->s_prcp[i];
  }
  ann_prcp = (sum_prcp/(double)ndays) * 365.25;
  if (ann_prcp == 0.0) ann_prcp = 1.0;
	
  /* Generate the effective annual precip, based on a 3-month
     moving-window. Requires some special case handling for the
     beginning of the record and for short records. */
  /* check if there are at least 90 days in this input file, if not,
     use a simple total scaled to effective annual precip */
  if (ndays < 90) {
    sum_prcp = 0.0;
    for (i=0 ; i<ndays ; i++) {
      sum_prcp += data->s_prcp[i];
    }
    effann_prcp = (sum_prcp/(double)ndays) * 365.25;
    /* if the effective annual precip for this period
       is less than 8 cm, set the effective annual precip to 8 cm
       to reflect an arid condition, while avoiding possible
       division-by-zero errors and very large ratios (PET/Pann) */
    if (effann_prcp < 8.0) 
      effann_prcp = 8.0;
    for (i=0 ; i<ndays ; i++) {
      parray[i] = effann_prcp;
    }
  }
  else {
    /* Check if the yeardays at beginning and the end of this input file
       match up. If so, use parts of the three months at the end
       of the input file to generate effective annual precip for
       the first 3-months. Otherwise, duplicate the first 90 days
       of the record. */
    start_yday = data->yday[0];
    end_yday = data->yday[ndays-1];
    if (start_yday != 1) {
      isloop = (end_yday == start_yday-1) ? 1 : 0;
    }
    else {
      isloop = (end_yday == 365 || end_yday == 366) ? 1 : 0;
    }
    
    /* fill the first 90 days of window */
    for (i=0 ; i<90 ; i++) {
      if (isloop) 
	window[i] = data->s_prcp[ndays-90+i];
      else 
	window[i] = data->s_prcp[i];
    }
    /* fill the rest of the window array */
    for (i=0 ; i<ndays ; i++) {
      window[i+90] = data->s_prcp[i];
    }
    
    /* for each day, calculate the effective annual precip from 
       scaled 90-day total */
    for (i=0 ; i<ndays ; i++)	{
      sum_prcp = 0.0;
      for (j=0 ; j<90 ; j++) {
	sum_prcp += window[i+j];
      }
      sum_prcp = (sum_prcp/90.0) * 365.25;
      /* if the effective annual precip for this 90-day period
	 is less than 8 cm, set the effective annual precip to 8 cm
	 to reflect an arid condition, while avoiding possible
	 division-by-zero errors and very large ratios (PET/Pann) */
      parray[i] = (sum_prcp < 8.0) ? 8.0 : sum_prcp;
    }
  } /* end if ndays >= 90 */	
  
  /*****************************************
   *                                       *
   * start of the main radiation algorithm *
   *                                       *
   *****************************************/
  
  /* before starting the iterative algorithm between humidity and 
     radiation, calculate all the variables that don't depend on 
     humidity so they only get done once. */
  
  /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
  t1 = 1.0 - (LR_STD * p->site_elev)/T_STD;
  t2 = G_STD / (LR_STD * (R/MA));
  pratio = pow(t1,t2);
  
  /* STEP (2) correct initial transmittance for elevation */ 
  trans1 = pow(TBASE,pratio);
  
  /* STEP (3) build 366-day array of ttmax0, potential rad, and daylength */
  
  /* precalculate the transcendentals */
  lat = p->site_lat;
  /* check for (+/-) 90 degrees latitude, throws off daylength calc */
  lat *= RADPERDEG;
  if (lat > 1.5707) 
    lat = 1.5707;
  if (lat < -1.5707) 
    lat = -1.5707;
  coslat = cos(lat);
  sinlat = sin(lat);
  cosslp = cos(p->site_slp * RADPERDEG);
  sinslp = sin(p->site_slp * RADPERDEG);
  cosasp = cos(p->site_asp * RADPERDEG);
  sinasp = sin(p->site_asp * RADPERDEG);
  /* cosine of zenith angle for east and west horizons */
  coszeh = cos(1.570796 - (p->site_ehoriz * RADPERDEG));
  coszwh = cos(1.570796 - (p->site_whoriz * RADPERDEG));
  
  /* sub-daily time and angular increment information */
  dt = SRADDT;                /* set timestep */ 
  dh = dt / SECPERRAD;        /* calculate hour-angle step */
  /* start vic_change */
  tinystepspday = 86400/SRADDT;
  /* end vic_change */

  /* begin loop through yeardays */
  for (i=0 ; i<365 ; i++) {
    /* calculate cos and sin of declination */
    decl = MINDECL * cos(((double)i + DAYSOFF) * RADPERDAY);
    cosdecl = cos(decl);
    sindecl = sin(decl);
    
    /* do some precalculations for beam-slope geometry (bsg) */
    bsg1 = -sinslp * sinasp * cosdecl;
    bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
    bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;
    
    /* calculate daylength as a function of lat and decl */
    cosegeom = coslat * cosdecl;
    sinegeom = sinlat * sindecl;
    coshss = -(sinegeom) / cosegeom;
    if (coshss < -1.0) 
      coshss = -1.0;  /* 24-hr daylight */
    if (coshss > 1.0) 
      coshss = 1.0;    /* 0-hr daylight */
    hss = acos(coshss);                /* hour angle at sunset (radians) */
    /* daylength (seconds) */
    daylength[i] = 2.0 * hss * SECPERRAD;
    
    /* start vic_change */
    if (daylength[i] > 86400)
      daylength[i] = 86400;
    /* end vic_change */
    
    /* solar constant as a function of yearday (W/m^2) */
    sc = 1368.0 + 45.5*sin((2.0*PI*(double)i/365.25) + 1.7);
    /* extraterrestrial radiation perpendicular to beam, total over
       the timestep (J) */
    dir_beam_topa = sc * dt;
    
    sum_trans = 0.0;
    sum_flat_potrad = 0.0;
    sum_slope_potrad = 0.0;
    
    /* begin sub-daily hour-angle loop, from -hss to hss */
    for (h=-hss ; h<hss ; h+=dh) {
      /* precalculate cos and sin of hour angle */
      cosh = cos(h);
      sinh = sin(h);
      
      /* calculate cosine of solar zenith angle */
      cza = cosegeom * cosh + sinegeom;
      
      /* calculate cosine of beam-slope angle */
      cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;
      
      /* check if sun is above a flat horizon */
      if (cza > 0.0) {
	/* when sun is above the ideal (flat) horizon, do all the
	   flat-surface calculations to determine daily total
	   transmittance, and save flat-surface potential radiation
	   for later calculations of diffuse radiation */
	
	/* potential radiation for this time period, flat surface,
	   top of atmosphere */
	dir_flat_topa = dir_beam_topa * cza;
	
	/* determine optical air mass */
	am = 1.0/(cza + 0.0000001);
	if (am > 2.9) {
	  ami = (int)(acos(cza)/RADPERDEG) - 69;
	  if (ami < 0) 
	    ami = 0;
	  if (ami > 20) 
	    ami = 20;
	  am = optam[ami];
	}
	
	/* correct instantaneous transmittance for this optical
	   air mass */
	trans2 = pow(trans1,am);
	
	/* instantaneous transmittance is weighted by potential
	   radiation for flat surface at top of atmosphere to get
	   daily total transmittance */
	sum_trans += trans2 * dir_flat_topa;
	
	/* keep track of total potential radiation on a flat
	   surface for ideal horizons */
	sum_flat_potrad += dir_flat_topa;
	
	/* keep track of whether this time step contributes to
	   component 1 (direct on slope) */
	if ((h<0.0 && cza>coszeh && cbsa>0.0) ||
	    (h>=0.0 && cza>coszwh && cbsa>0.0)) {
	  
	  /* sun between east and west horizons, and direct on
	     slope. this period contributes to component 1 */
	  sum_slope_potrad += dir_beam_topa * cbsa;
	}
	
      } /* end if sun above ideal horizon */
      else dir_flat_topa = -1;

      /* start vic_change */
      tinystep = (12L * 3600L + h * SECPERRAD)/SRADDT;
      if (tinystep < 0)
	tinystep = 0;
      if (tinystep > tinystepspday-1)
	tinystep = tinystepspday-1;
      if (dir_flat_topa > 0)
	tiny_radfract[i][tinystep] = dir_flat_topa;
      else
	tiny_radfract[i][tinystep] = 0;
      /* end vic_change */
      
    } /* end of sub-daily hour-angle loop */
    
    /* start vic_change */
    if (daylength[i] && sum_flat_potrad > 0) {
      for (j = 0; j < tinystepspday; j++)
	tiny_radfract[i][j] /= sum_flat_potrad;
    }
    /* end vic_change */

    /* calculate maximum daily total transmittance and daylight average
       flux density for a flat surface and the slope */
    if (daylength[i]) {
      ttmax0[i] = sum_trans / sum_flat_potrad;
      flat_potrad[i] = sum_flat_potrad / daylength[i];
      slope_potrad[i] = sum_slope_potrad / daylength[i];
    }
    else {
      ttmax0[i] = 0.0;
      flat_potrad[i] = 0.0;
      slope_potrad[i] = 0.0;
    }

  } /* end of i=365 days loop */
  
  /* force yearday 366 = yearday 365 */
  ttmax0[365] = ttmax0[364];
  flat_potrad[365] = flat_potrad[364];
  slope_potrad[365] = slope_potrad[364];
  daylength[365] = daylength[364];
  
  /* start vic_change */
  for (j = 0 ; j < tinystepspday; j++)
    tiny_radfract[365][j] = tiny_radfract[364][j];
  /* end vic_change */

  /* STEP (4)  calculate the sky proportion for diffuse radiation */
  /* uses the product of spherical cap defined by average horizon angle
     and the great-circle truncation of a hemisphere. this factor does not
     vary by yearday. */
  avg_horizon = (p->site_ehoriz + p->site_whoriz)/2.0;
  horizon_scalar = 1.0 - sin(avg_horizon * RADPERDEG);
  if (p->site_slp > avg_horizon) 
    slope_excess = p->site_slp - avg_horizon;
  else 
    slope_excess = 0.0;
  if (2.0*avg_horizon > 180.0) 
    slope_scalar = 0.0;
  else {
    slope_scalar = 1.0 - (slope_excess/(180.0 - 2.0*avg_horizon));
    if (slope_scalar < 0.0) 
      slope_scalar = 0.0;
  }
  sky_prop = horizon_scalar * slope_scalar;
  
  /* b parameter, and t_fmax not varying with Tdew, so these can be
     calculated once, outside the iteration between radiation and humidity
     estimates. Requires storing t_fmax in an array. */
  for (i=0 ; i<ndays ; i++) {	
    /* b parameter from 30-day average of DTR */
    b = B0 + B1 * exp(-B2 * sm_dtr[i]);
    
    /* proportion of daily maximum transmittance */
    t_fmax[i] = 1.0 - 0.9 * exp(-b * pow(dtr[i],C));
    
    /* correct for precipitation if this is a rain day */
    if (data->prcp[i] > options.SW_PREC_THRESH) t_fmax[i] *= RAIN_SCALAR;
    data->s_tfmax[i] = t_fmax[i];

  }
 
  /* Initial values of vapor pressure, etc */
  if (ctrl->indewpt) {
    /* Observed Tdew supplied */
    for (i=0 ; i<ndays ; i++) {
      tdew[i] = data->tdew[i];
    }
  }
  else {
    /* Estimate Tdew */
    for (i=0 ; i<ndays ; i++) {
      tdew[i] = data->s_tmin[i];
    }
  }
  if (ctrl->invp) {
    /* Observed vapor pressure supplied */
    for (i=0 ; i<ndays ; i++) {
      pva[i] = data->s_hum[i];
    }
  }
  else {
    /* convert dewpoint to vapor pressure */
    for (i=0 ; i<ndays ; i++) {
      /* start vic_change */
      /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
      pva[i] = svp(tdew[i]);
      /* end vic_change */
    }
  }

  /* Other values needed for srad_humidity calculation */
  pa = atm_pres(p->site_elev);
  for (i=0 ; i<ndays ; i++) {
    yday = data->yday[i]-1;
    data->s_dayl[i] = daylength[yday];
    tdew_save[i] = tdew[i];
    pva_save[i] = pva[i];
  }
  
  /* Initial estimates of solar radiation, cloud fraction, etc. */
  compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength, pet, parray, pa, dtr);

  /* estimate annual PET */
  sum_pet = 0.0;
  for (i=0 ; i<ndays ; i++) {
    sum_pet += pet[i];
  }
  ann_pet = (sum_pet/(double)ndays) * 365.25;

  /* Reset humidity terms if no iteration desired */
  if (ctrl->indewpt || ctrl->invp || (options.VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) ) {
    for (i=0 ; i<ndays ; i++) {
      tdew[i] = tdew_save[i];
      pva[i] = pva_save[i];
    }
  }

  /* Set up srad-humidity iterations */
  if (options.VP_ITER == VP_ITER_ALWAYS || (options.VP_ITER == VP_ITER_ANNUAL && ann_pet/ann_prcp >= 2.5) || options.VP_ITER == VP_ITER_CONVERGE) {
    //printf("Using arid-climate humidity algorithm\n");
    if (options.VP_ITER == VP_ITER_CONVERGE) {
      max_iter = 100;
    }
    else {
      max_iter = 2;
    }
  }
  else {
    //printf("Using Tdew=Tmin humidity algorithm\n");
    max_iter = 1;
  }

  /* srad-humidity iterations */
  tol = 0.01;
  iter = 1;
  rmse_tdew = tol+1;
  while (rmse_tdew > tol && iter < max_iter) {
    update_pva = 1;
    for (i=0 ; i<ndays ; i++) {
      tdew_save[i] = tdew[i];
    }
    compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength, pet, parray, pa, dtr);
    rmse_tdew = 0;
    for (i=0 ; i<ndays ; i++) {
      rmse_tdew += (tdew[i]-tdew_save[i])*(tdew[i]-tdew_save[i]);
    }
    rmse_tdew /= ndays;
    rmse_tdew = pow(rmse_tdew,0.5);
    iter++;
  }

  /* save humidity in output data structure */
  if (ctrl->outhum) {
    if (!ctrl->invp) {
      for (i=0 ; i<ndays ; i++) {
        data->s_hum[i] = pva[i];
      }
    }
  }
  else {
    /* output humidity as vapor pressure deficit (Pa) */
    for (i=0 ; i<ndays ; i++) {
      /* calculate saturated VP at tday */
      /* start vic_change */
      /* pvs = 610.7 * exp(17.38 * data->s_tday[i]/(239.0+data->s_tday[i])); */
      pvs = svp(data->s_tday[i]);
      /* end vic_change */
      vpd = pvs - pva[i];
      if (vpd < 0.0) vpd = 0.0;
      data->s_hum[i] = vpd;
    }
  }

  /* free local array memory */
  free(dtr);
  free(sm_dtr);
  free(parray);
  free(window);
  free(t_fmax);
  free(tdew);
  free(pet);
  free(pva);
  free(tdew_save);
  free(pva_save);

  return (!ok);

} /* end of calc_srad_humidity_iterative() */

/* New function, not originally part of MTCLIM code */
void compute_srad_humidity_onetime(int ndays, const control_struct *ctrl, data_struct *data, double *tdew, double *pva, double *ttmax0, double *flat_potrad, double *slope_potrad, double sky_prop, double *daylength, double *pet, double *parray, double pa, double *dtr) {

  extern option_struct options;
  int i;
  int yday;
  double t_tmax;
  double t_final;
  double pdif;
  double pdir;
  double srad1;
  double srad2;
  double sc;
  double potrad;
  double tmink;
  double ratio;
  double ratio2;
  double ratio3;
  double tdewk;

  for (i=0 ; i<ndays ; i++) {

    yday = data->yday[i]-1;

    /*** Compute SW radiation ***/

    t_tmax = ttmax0[yday] + ABASE * pva[i];
    if (t_tmax < 0.0001) t_tmax = 0.0001; // this is mainly for the case of observed VP supplied, for which t_tmax sometimes ends up being negative (when potential radiation is low and VP is high)
    data->s_ttmax[i] = t_tmax;

    /* final daily total transmittance */
    t_final = t_tmax * data->s_tfmax[i];

    /* estimate fraction of radiation that is diffuse, on an
    instantaneous basis, from relationship with daily total
    transmittance in Jones (Plants and Microclimate, 1992)
    Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
    Fig 6.14, p. 122. */
    pdif = -1.25*t_final + 1.25;
    if (pdif > 1.0) pdif = 1.0;
    if (pdif < 0.0) pdif = 0.0;

    /* estimate fraction of radiation that is direct, on an
    instantaneous basis */
    pdir = 1.0 - pdif;

    /* the daily total radiation is estimated as the sum of the
    following two components:
    1. The direct radiation arriving during the part of
       the day when there is direct beam on the slope.
    2. The diffuse radiation arriving over the entire daylength
       (when sun is above ideal horizon).
    */

    /* component 1 */
    srad1 = slope_potrad[yday] * t_final * pdir;

    /* component 2 (diffuse) */
    /* includes the effect of surface albedo in raising the diffuse
    radiation for obstructed horizons */
    srad2 = flat_potrad[yday] * t_final * pdif * (sky_prop + DIF_ALB*(1.0-sky_prop)); 

    /* snow pack influence on radiation */
    if (options.MTCLIM_SWE_CORR && data->s_swe[i] > 0.0) {
      /* snow correction in J/m2/day */
      sc = (1.32 + 0.096 * data->s_swe[i]) * 1e6;
      /* convert to W/m2 and check for zero daylength */
      if (daylength[yday] > 0.0) sc /= daylength[yday];
      else sc = 0.0;
      /* set a maximum correction of 100 W/m2 */
      if (sc > 100.0) sc = 100.0;
    }
    else sc = 0.0;

    /* save daily radiation */
    /* save cloud transmittance when rad is an input */
    if (ctrl->insw) {
      potrad = (srad1+srad2+sc)*daylength[yday]/t_final/86400;
      if (potrad>0 && data->s_srad[i]>0 && daylength[yday]>0) {
	data->s_tfmax[i] = (data->s_srad[i])/(potrad*t_tmax);//both of these are 24hr mean rad. here
	if (data->s_tfmax[i] > 1.0) data->s_tfmax[i] = 1.0;
      }
      else {
	data->s_tfmax[i] = 1.0;
      }
    }
    else {
      data->s_srad[i] = srad1 + srad2 +sc;
    }

    /* start vic_change */
    if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
      data->s_tskc[i] = (1.-data->s_tfmax[i]);
    }
    else {
      data->s_tskc[i] = sqrt((1.-data->s_tfmax[i])/0.65);
    }
    data->s_fdir[i] = pdir;
    /* end vic_change */

  }

  /*** Compute PET using SW radiation estimate, and update Tdew, pva ***/
  for (i=0 ; i<ndays ; i++) {

    tmink = data->s_tmin[i] + KELVIN;
    pet[i] = calc_pet(data->s_srad[i],data->s_tday[i],pa,data->s_dayl[i]);

    /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
    ratio = pet[i]/parray[i];
    data->s_ppratio[i] = ratio*365.25;
    ratio2 = ratio*ratio;
    ratio3 = ratio2*ratio;
    tdewk = tmink*(-0.127 + 1.121*(1.003 - 1.444*ratio + 12.312*ratio2 
                - 32.766*ratio3) + 0.0006*(dtr[i]));
    tdew[i] = tdewk - KELVIN;

    /* start vic_change */
    /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
    pva[i] = svp(tdew[i]);
    /* end vic_change */

  }

  return;

}

/* data_free frees the memory previously allocated by data_alloc() */
int data_free(const control_struct *ctrl, data_struct *data)
{
  int ok=1;
  if (ctrl->inyear) free(data->year);
  free(data->yday);
  free(data->tmax);
  free(data->tmin);
  free(data->prcp);
  if (ctrl->indewpt) free(data->tdew);
  free(data->s_tmax);
  free(data->s_tmin);
  free(data->s_tday);
  free(data->s_prcp);
  free(data->s_hum);
  free(data->s_srad);
  free(data->s_dayl);
  free(data->s_swe);
  /* start vic_change */
  free(data->s_fdir);
  free(data->s_tskc);
  free(data->s_ppratio);
  free(data->s_tfmax);
  free(data->s_ttmax);
  /* end vic_change */
  return (!ok);
}

/* calc_pet() calculates the potential evapotranspiration for aridity 
   corrections in calc_vpd(), according to Kimball et al., 1997 */
double calc_pet(double rad, double ta, double pa, double dayl)
{
  /* input parameters and units :
     double rad      (W/m2)  daylight average incident shortwave radiation
     double ta       (deg C) daylight average air temperature
     double pa       (Pa)    air pressure
     double dayl     (s)     daylength 
  */
  
  double rnet;       /* (W m-2) absorbed shortwave radiation avail. for ET */
  double lhvap;      /* (J kg-1) latent heat of vaporization of water */ 
  double gamma;      /* (Pa K-1) psychrometer parameter */
  double dt = 0.2;   /* offset for saturation vapor pressure calculation */
  double t1, t2;     /* (deg C) air temperatures */
  double pvs1, pvs2; /* (Pa)   saturated vapor pressures */
  double pet;        /* (kg m-2 day-1) potential evapotranspiration */
  double s;          /* (Pa K-1) slope of saturated vapor pressure curve */
  
  /* calculate absorbed radiation, assuming albedo = 0.2  and ground
     heat flux = 10% of absorbed radiation during daylight */
  rnet = rad * 0.72;
  
  /* calculate latent heat of vaporization as a function of ta */
  lhvap = 2.5023e6 - 2430.54 * ta;
  
  /* calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
     where:
     cp       (J/kg K)   specific heat of air
     epsilon  (unitless) ratio of molecular weights of water and air
  */
  gamma = CP * pa / (lhvap * EPS);
  
  /* estimate the slope of the saturation vapor pressure curve at ta */
  /* temperature offsets for slope estimate */
  t1 = ta+dt;
  t2 = ta-dt;
  
  /* calculate saturation vapor pressures at t1 and t2, using formula from 
     Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
     Meteorol. Mag., 114:49-56.
  */
  /* start vic_change */
  /* pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1)); */
  pvs1 = svp(t1);
  /* pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2)); */
  pvs2 = svp(t2);
  /* end vic_change */
  
  /* calculate slope of pvs vs. T curve near ta */
  s = (pvs1-pvs2) / (t1-t2);
  
  /* calculate PET using Priestly-Taylor approximation, with coefficient
     set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day */
  pet = (1.26 * (s/(s+gamma)) * rnet * dayl)/lhvap;
  
  /* return a value in centimeters/day, because this value is used in a ratio
     to annual total precip, and precip units are centimeters */
  return (pet/10.0);
}    

/* atm_pres() calculates the atmospheric pressure as a function of elevation */
double atm_pres(double elev)
{
  /* daily atmospheric pressure (Pa) as a function of elevation (m) */
  /* From the discussion on atmospheric statics in:
     Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
     Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
     (p. 168)
  */
  
  int ok=1;
  double t1,t2;
  double pa;
  
  t1 = 1.0 - (LR_STD * elev)/T_STD;
  t2 = G_STD / (LR_STD * (R / MA));
  pa = P_STD * pow(t1,t2);
  
  return(pa);
}

/* pulled_boxcar() calculates a moving average of antecedent values in an
   array, using either a ramped (w_flag=1) or a flat (w_flag=0) weighting */	
int pulled_boxcar(double *input,double *output,int n,int w,int w_flag)
{
  int ok=1;
  int tail,i,j;
  double *wt;
  double total,sum_wt;
  
  if (w > n) {
    printf("Boxcar longer than array...\n");
    printf("Resize boxcar and try again\n");
    ok=0;
  }
  
  if (ok && !(wt = (double*) malloc(w * sizeof(double)))) {
    printf("Allocation error in boxcar()\n");
    ok=0;
  }
  
  if (ok) {
    /* when w_flag != 0, use linear ramp to weight tails,
       otherwise use constant weight */
    sum_wt = 0.0;
    if (w_flag) {
      for (i=0 ; i<w ; i++) {
	wt[i] = (double)(i+1);
	sum_wt += wt[i];
      }
    }
    else {
      for (i=0 ; i<w ; i++) { 	
	wt[i] = 1.0;
	sum_wt += wt[i];
      }
    }
    
    /* fill the output array, starting with the point where a full
       boxcar can be calculated */
    for (i=w-1 ; i<n ; i++) {
      total = 0.0;
      for (j=0 ; j<w ; j++) {
	total += input[i-w+j+1] * wt[j];
      }
      output[i] = total/sum_wt;
    }
    
    /* fill the first w elements of the output array with the value from
       the first full boxcar */
    for (i=0 ; i<w-1 ; i++) {
      output[i] = output[w-1];
    }
    
    free(wt);
    
  } /* end if ok */
  
  return (!ok);
}
/* end of pulled_boxcar() */  

