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
  10-27-00 Modified to use vapor pressure in Pa instead of kPa   KAC
***************************************************************************/
{
  /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
     Addison-Wesley, 1990, p. 42-45 */
  double emissivity;
  double cloudfactor;
  
  emissivity = 0.740 + 0.0049 * vp / 100.; /* Bras 2.35 */
  cloudfactor = 1.0 + 0.17 * tskc * tskc; /* Bras 2.43 */
  air_temp += KELVIN;
  
  *longwave = emissivity * cloudfactor * STEFAN_B * air_temp * air_temp *
    air_temp * air_temp / LWAVE_COR; 

}

void calc_netlongwave(double *longwave, 
		      double  tskc, 
		      double  air_temp, 
		      double  vp)
/****************************************************************************
  Function to calculate the daily net incoming longwave radiation in W/m2 
  (Bras 1990).  Uses the fractional cloud cover (tskc), air temperature (C)
  and vapor pressure (Pa).

  Modifcations:
  10-25-00 Modifeid to use vp in Pa instead of kPa               KAC
****************************************************************************/
{
  double emissivity;            /* emissivity of the atmosphere */
  double cloudfactor;           /* cloudiness correction factor */
  
  emissivity = 0.740 + 0.0049 * vp / 100.; /* Bras 2.35 */
  cloudfactor = 1.0 + 0.17 * tskc * tskc; /* Bras 2.43 */
  air_temp += KELVIN;
 
  *longwave = (emissivity * cloudfactor / LWAVE_COR - 1 ) * STEFAN_B 
    * air_temp * air_temp * air_temp * air_temp ; 
}
