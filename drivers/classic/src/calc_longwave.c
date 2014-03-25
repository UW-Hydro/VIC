#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void calc_longwave(double *longwave, 
		   double  tskc, 
		   double  air_temp, 
		   double  vp)
/***************************************************************************
  This routine estimates long wave radiation based on the fractional cloud
  cover (tskc), the current air temperature (C), and the atmospheric vapor
  pressure (Pa).

  Modifications:
  2000-Oct-27 Modified to use vapor pressure in Pa instead of kPa.	KAC
  2011-Nov-04 Added new option for handling the effect of clouds.	TJB
  2011-Nov-04 Added alternative clear-sky longwave algorithm.		TJB
  2012-Apr-13 Changed relationship between cloud_fraction and tskc for
	      LW_CLOUD==LW_CLOUD_DEARDORFF, to account for new, simplified
	      meaning of tskc.						TJB
  2013-Dec-26 Removed LWAVE_COR.					TJB
***************************************************************************/
{
  extern option_struct options;
  double emissivity;
  double emissivity_clear;
  double cloudfactor;
  double cloudfrac;
  double x;
  
  air_temp += KELVIN; // convert to Kelvin
  vp /= 100; // convert to mbar

  /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
     Addison-Wesley, 1990, p. 42-45 */

  if (options.LW_TYPE == LW_TVA) {
    emissivity_clear = 0.740 + 0.0049 * vp; // TVA (1972) - see Bras 2.35
  }
  else if (options.LW_TYPE == LW_ANDERSON) {
    emissivity_clear = 0.68 + 0.036*pow(vp,0.5); // Anderson (1964)
  }
  else if (options.LW_TYPE == LW_BRUTSAERT) {
    x = vp/air_temp;
    emissivity_clear = 1.24*pow(x,0.14285714); // Brutsaert (1975)
  }
  else if (options.LW_TYPE == LW_SATTERLUND) {
    emissivity_clear = 1.08*(1-exp(-1*pow(vp,(air_temp/2016)))); // Satterlund (1979)
  }
  else if (options.LW_TYPE == LW_IDSO) {
    emissivity_clear = 0.7 + 5.95e-5*vp*exp(1500/air_temp); // Idso (1981)
  }
  else if (options.LW_TYPE == LW_PRATA) {
    x = 46.5*vp/air_temp;
    emissivity_clear = 1 - (1+x)*exp(-1*pow((1.2+3*x),0.5)); // Prata (1996)
  }

  if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
    /* Assume emissivity of clouds is 1.0, and that total emissivity is weighted
       average of cloud emission plus clear-sky emission, weighted by fraction of
       sky occupied by each (method of Deardorff, 1978) */
    cloudfrac = tskc;
    emissivity = cloudfrac*1.0 + (1-cloudfrac)*emissivity_clear; // Deardorff (1978)
  }
  else {
    cloudfactor = 1.0 + 0.17 * tskc * tskc; // see Bras 2.43
    emissivity = cloudfactor*emissivity_clear;
  }

  *longwave = emissivity * STEFAN_B * air_temp * air_temp * air_temp * air_temp; 

}

