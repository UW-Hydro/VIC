#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double estimate_dew_point(dmy_struct *dmy,
			  double     *yearly_prec,
			  double      tmin,
			  double      tmax,
			  double      priest,
			  double      day_len,
			  int         day) {
/*******************************************************************
  estimate_dew_point        Keith Cherkauer      November 11, 1998

  This routine was written to replace the pre-publication version
  of dew point temperature estimation routine from Kimball et. al., 
  with the version that was published.  The set-up code for this
  routine which was orginally in rad_and_vpd.c has been relocated 
  into the initialize_atmos.c routines which now handle initializing
  all of the atmospheric variables.  

  As described in Kimball, J. S., S. W. Running, and R. Nemani, "An
  improved method for estimating surface humidity from daily minimum
  temperature", Agricultural and Forest Engineering, 85(1997) 87-98.

  Variables:
  dmy_struct *dmy          N/A  Structure containing date information
  double     *yearly_prec  mm/year  Precipitation per year
  double      tmin         C        Minimum daily temperature
  double      tmax         C        Maximum daily temperature 
  double      priest       kg-m^2/s Priestly potential evaporation for the day
  double      day_len      s        Length of the current day in seconds
  int         day          N/A      Day in simulation

*******************************************************************/

  double Tdew;
  double EF;

  EF = ((priest) * day_len) / (yearly_prec[dmy[day].year-dmy[0].year]);
  Tdew = (tmin + KELVIN) * ( -0.127 + 1.121 * (1.003 - 1.444 * EF
					       + 12.312 * EF * EF
					       -32.766 * EF * EF * EF) 
			     + 0.0006 * (tmax - tmin) );

  return(Tdew-KELVIN);
}
