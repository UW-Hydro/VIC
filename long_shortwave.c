#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

/*******************************************************************************
*
  Function to calculate the daily incoming shortwave in W/m2 (Shuttleworth 1993)
*******************************************************************************/
double in_shortwave(float lat, int day, double trans)
{
  double omega_s;               /* sunset hour angle in radians */
  double declination;           /* solar declination in radians */
  double extra_solar;           /* extraterrestrial solar radiation (W/m^2) */
  double dist;                  /* relative distance bewteen earth and sun */
  double lat_rad;               /* latitude in radians */

  declination = 0.4093 * sin(2.0 * PI/DAYS_PER_YEAR * day - 1.405);
  lat_rad = lat * DtoR;
 
  omega_s = acos(-tan(lat_rad)*tan(declination));
 
  dist = 1.0 + 0.033 * cos(2.0 * PI/DAYS_PER_YEAR * day);
 
  extra_solar = SOLAR_CONSTANT / PI * dist * 
    (omega_s * sin(declination) * sin(lat_rad) + 
     cos(declination) * cos(lat_rad) * sin(omega_s));
  
  return (trans*extra_solar);
}
 
/*******************************************************************************
*
  Function to calculate the daily net outgoing longwave radiation in W/m2 
  (Bras 1990)
*******************************************************************************/
double net_out_longwave(double trans,
                        double trans_clear,
                        double tair,
                        double vapor,
                        double *cloudiness)
{
  double emissivity;            /* emissivity of the atmosphere */
  double cloudfactor;           /* cloudiness correction factor */
  double longwave;              /* net ourgoing longwave */
  double t_kelvin;              /* temperature in Kelvin */
  double v_mbar;                /* vapor pressure in mbar */
 
  t_kelvin = tair + 273.15;
  v_mbar = vapor/100.0;
 
  *cloudiness = 1.0/0.65 * sqrt(1.0 - trans/trans_clear);
  *cloudiness = (*cloudiness > 1.0) ? 1.0 : *cloudiness;
  emissivity = 0.70 + 0.0000595 * v_mbar * exp(1500.0/t_kelvin);
  cloudfactor = 1.0 + 0.17 * pow(*cloudiness,2.0);
  
  longwave = (1.0 - emissivity * cloudfactor) * STEFAN * pow(t_kelvin,4.0);
  
/* printf("%lf  %lf  %lf  %lf  %lf  %lf\n",cloudiness,emissivity,
	 longwave,tair,trans/trans_clear,vapor);*/
  /* return longwave in W/m^2 */
  
  return longwave;
}


