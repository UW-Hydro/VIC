#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>
 
static char vcid[] = "";

void calc_longwave(double *longwave, double tskc, double air_temp, double vp)
{
  /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
     Addison-Wesley, 1990, p. 42-45 */
  double emissivity;
  double cloudfactor;
  
  emissivity = 0.740 + 0.0049 * vp * 10.0; /* Bras 2.35 */
  cloudfactor = 1.0 + 0.17 * tskc * tskc; /* Bras 2.43 */
  air_temp += KELVIN;
  
  *longwave = emissivity * cloudfactor * STEFAN_B * air_temp * air_temp *
    air_temp * air_temp / LWAVE_COR; 

}

/****************************************************************************
  Function to calculate the daily net incoming longwave radiation in W/m2 
  (Bras 1990)
****************************************************************************/
void calc_netlongwave(double *longwave, double tskc, double air_temp, double vp)
{
  double emissivity;            /* emissivity of the atmosphere */
  double cloudfactor;           /* cloudiness correction factor */
  
  emissivity = 0.740 + 0.0049 * vp * 10.0; /* Bras 2.35 */
  cloudfactor = 1.0 + 0.17 * tskc * tskc; /* Bras 2.43 */
  air_temp += KELVIN;
 
  *longwave = (emissivity * cloudfactor / LWAVE_COR - 1 ) * STEFAN_B 
    * air_temp * air_temp * air_temp * air_temp ; 
}
