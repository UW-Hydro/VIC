/*
  filename  : rad_and_vpd.h
  purpose   : include file for rad_and_vpd.h, part of the water balance version
              of the VIC-2L model.
  programmer: Bart Nijssen
  date      : August 3, 1995
  changes   :
*/

/* define statements */

/* miscellaneous constants */
#ifndef MAXSTRING
#define MAXSTRING 512
#endif

#ifndef MINSTRING
#define MINSTRING 20
#endif

#define DAYS_PER_YEAR 365.
#define MAX_YEARS 100
#define DtoR 0.017453293	/* degrees to radians */
#define PI 3.1415927
#define STEFAN 5.6696e-8	/* Stefan boltzmann constant */
#define SOLAR_CONSTANT 1400.0	/* Solar constant in W/m^2 */
#define SEC_PER_DAY 86400.	/* seconds per day */

/* constants for daily transmittance, values from Kimball et al., 1995 */
#define A1_TRANS 0.60		/* maximum clear sky transmittance at MSL */
#define A2_TRANS 0.0000295	/* lapse rate per meter of maximum clear sky
				   transmittance */
#define B_TRANS	-0.0030		/* empirical constant */
#define C_TRANS	2.4		/* empirical constant */

/* constants for Priestley-Taylor potential evaporation, values from Kimball 
   et al., 1995 */
#define ALPHA_PT 1.26		/* Priestley-Taylor parameter */
#define GAMMA_PT 6.6		/* psychrometric constant (kPa/K) */
#define LV_PT    2.5E6		/* latent heat of vaporization (J/kg) */

/* define constants for saturated vapor pressure curve (kPa) */
#define A_SVP 0.61078
#define B_SVP 17.269
#define C_SVP 237.3

/* define constants for penman evaporation */
#define CP_PM 1013		/* specific heat of moist air J/kg/C 
				   (Handbook of Hydrology) */
#define PS_PM 101300		/* sea level air pressure in Pa */
#define LAPSE_PM -0.006		/* environmental lapse rate in C/m */

/* define constants for regression equation to correct vapor pressure estimate 
   based on minimum temperature */

#define HUMID_RATIO 2.25	/* ratio of total yearly potential evaporation 
				   to total yearly precipitation*/




