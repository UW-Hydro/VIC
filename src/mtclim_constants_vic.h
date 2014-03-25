/*
  Changes for the VIC implementation are preceded by the comment * start
  vic_change * and followed by the comment * end vic_change * */

/* RCS Id String
 * $Id$
 */

/* 
mtclim_constants_vic.h
constants typedefs, and function prototypes for MTCLIM 4.3

Peter Thornton
NTSG, School of Forestry
University of Montana
5/10/98 

(dim) stands for dimensionless values

Adapted for inclusion in VIC-code:
Bart Nijssen
Sat Aug 21 16:58:43 1999

  Modified:
  2011-Nov-04 Updated to MTCLIM 4.3				TJB
  2012-Feb-16 Removed calc_srad_humidity().			TJB
  2013-Jul-25 Added data->s_fdir.				TJB

*/



#define TDAYCOEF 0.45     /* daylight air temperature coefficient (dim) */

#define SECPERRAD 13750.9871     /* seconds per radian of hour angle */
#define RADPERDAY 0.017214       /* radians of Earth orbit per julian day */
#define RADPERDEG 0.01745329     /* radians per degree */
#define MINDECL -0.4092797       /* minimum declination (radians) */
#define DAYSOFF 11.25            /* julian day offset of winter solstice */
/* start vic_change */
#define SRADDT 30.0             /* timestep for radiation routine (seconds) */
				/* Note:  Make sure that 3600 % SRADDT == 0 */
/* end vic_change */

#define MA       28.9644e-3      /* (kg mol-1) molecular weight of air */
#define MW       18.0148e-3      /* (kg mol-1) molecular weight of water */
#define R        8.3143          /* (m3 Pa mol-1 K-1) gas law constant */
#define G_STD    9.80665         /* (m s-2) standard gravitational accel. */ 
#define P_STD    101325.0        /* (Pa) standard pressure at 0.0 m elevation */
#define T_STD    288.15          /* (K) standard temp at 0.0 m elevation  */ 
#define CP       1010.0          /* (J kg-1 K-1) specific heat of air */
#define LR_STD   0.0065          /* (-K m-1) standard temperature lapse rate */
/* start vic_change */
#ifndef PI
#define PI       3.14159265
#endif
/* end vic_change */

/****************************
 **                         ** 
 **  STRUCTURE DEFINITIONS  **
 **                         **
 ****************************/
typedef struct
{
  int ndays;             /* number of days of data in input file */
  int insw;              /* input shortwave radiation flag (0=NO, 1=YES) */
  int indewpt;           /* input dewpoint temperature flag (0=NO, 1=YES) */
  int invp;              /* input vapor pressure flag (0=NO, 1=YES) */
  int outhum;            /* output humidity flag            (0=VPD, 1=VP) */
  int inyear;            /* input year flag                 (0=NO, 1=YES) */
} control_struct;

typedef struct
{
  double base_elev;      /* base elevation, meters */
  double base_isoh;      /* base annual precip isohyet, cm */
  double site_lat;       /* site latitude, dec. degrees (- for south) */
  double site_elev;      /* site elevation, meters */
  double site_slp;       /* site slope, degrees */
  double site_asp;       /* site aspect, degrees */
  double site_isoh;      /* site annual precip isohyet, cm */
  double site_ehoriz;    /* site east horizon, degrees */
  double site_whoriz;    /* site west horizon, degrees */
  double tmax_lr;        /* maximum temperature lapse rate, deg C/1000m */
  double tmin_lr;        /* minimum temperature lapse rate, deg C/1000m */
} parameter_struct;

typedef struct
{
  int *year;             /* array of year values */
  int *yday;             /* array of yearday values */
  double *tmax;          /* array of base maximum temperature values */
  double *tmin;          /* array of base minimum temperature values */
  double *prcp;          /* array of base daily precipitation values */
  double *tdew;          /* array of base dewpoint temperature values */
  double *s_tmax;        /* array of site tmax values */
  double *s_tmin;        /* array of site tmin values */
  double *s_tday;        /* array of site daylight temperature values */
  double *s_prcp;        /* array of site prcp values */
  double *s_hum;         /* array of site humidity values (VPD or VP, Pa) */
  double *s_srad;        /* array of site shortwave radiation values */
  double *s_dayl;        /* array of site daylength values */
  double *s_swe;         /* array of site snowpack values */
  /* start vic_change */
  double *s_fdir;	 /* array of site values of direct fraction of shortwave radiation */
  double *s_tskc;	 /* array of cloudiness values */
  double *s_ppratio; /* array of pet/prcp ratio values */
  double *s_ttmax; /* array of clear sky transmittance values */
  double *s_tfmax; /* array of cloud transmittance factor values */
  /* end vic_change */
} data_struct;

/********************************
 **                             **
 **    FUNCTION PROTOTYPES      **
 **                             **
 ********************************/
int calc_tair(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data);
int calc_prcp(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data);
/* start vic_change */
int calc_srad_humidity_iterative(const control_struct *ctrl,
				 const parameter_struct *p, data_struct *data,
				 double **tiny_radfract);
int snowpack(const control_struct *ctrl, const parameter_struct *p, 
	      data_struct *data);
void compute_srad_humidity_onetime(int ndays, const control_struct *ctrl, data_struct *data, double *tdew, double *pva, double *ttmax0, double *flat_potrad, double *slope_potrad, double sky_prop, double *daylength, double *pet, double *parray, double pa, double *dtr);
/* end vic_change */
int data_alloc(const control_struct *ctrl, data_struct *data);
int data_free(const control_struct *ctrl, data_struct *data);
double calc_pet(double rad, double ta, double pa, double dayl);
double atm_pres(double elev);
int pulled_boxcar(double *input,double *output,int n,int w,int w_flag);

