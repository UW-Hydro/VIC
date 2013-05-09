/********************************************************************************
  filename  : penman.c
  purpose   : Calculate daily evapotranspiration using the combination equation
  interface : - input : - net radiation (W/m2)
                        - vapor pressure deficit (Pa)
			- aerodynamical resistance (s/m)
			- minimum stomatal resistance (s/m)
			- architectural resistance (s/m)
			- LAI
			- soil moisture stress factor
			- air temperature (to calculate slope of saturated 
			  vapor pressure curve) (C)
			- elevation (m)
              - output: - daily evapotranspiration (mm/day)
  programmer: Bart Nijssen
  date      : August 7, 1995
  references: 
  changes   :
  04-25-03 modified to use pressure in Pa directly, rather than
           converting it within the routine.              KAC
  2004-May-07 Changed
		if (vpd > 0.0 && evap < 0.0)
	      to
		if (vpd >= 0.0 && evap < 0.0)
	      to correctly handle evap when vpd == 0.0.			TJB
  2009-Jun-09 Moved calculation of canopy resistance from penman() into
	      separate function calc_rc().  Also allowed for reference
	      crop case in calc_rc().					TJB
  2009-Jun-09 Removed unnecessary functions quick_penman() and
	      compute_penman_constants().				TJB

********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#define CLOSURE 4000		/** Pa **/
#define RSMAX 5000
#define VPDMINFACTOR 0.1

double calc_rc(double rs, 
	       double net_short, 
	       float  RGL, 
	       double tair, 
	       double vpd, 
	       double lai, 
	       double gsm_inv, 
	       char   ref_crop)
{
  double f;
  double DAYfactor;             /* factor for canopy resistance based on photosynthesis */
  double Tfactor;               /* factor for canopy resistance based on temperature */
  double vpdfactor;             /* factor for canopy resistance based on vpd */
  double rc;

  if (rs == 0) {
    rc = 0;
  }
  else if (lai == 0) {
    rc = HUGE_RESIST;
  }
  else if (ref_crop) {
    /* calculate simple reference canopy resistance in s/m */
    rc = rs/(lai * 0.5);
  }
  else {
    /* calculate resistance factors (Wigmosta et al., 1994) */
    if(rs>0.) {
      f = net_short / RGL;
      DAYfactor = (1. + f)/(f + rs/RSMAX);
    }
    else DAYfactor = 1.;

    Tfactor = .08 * tair - 0.0016 * tair * tair;
    Tfactor = (Tfactor <= 0.0) ? 1e-10 : Tfactor;

    vpdfactor = 1 - vpd/CLOSURE;
    vpdfactor = (vpdfactor < VPDMINFACTOR) ? VPDMINFACTOR : vpdfactor;

    /* calculate canopy resistance in s/m */
    rc = rs/(lai * gsm_inv * Tfactor * vpdfactor) * DAYfactor;
    rc = (rc > RSMAX) ? RSMAX : rc;
  }

  return rc;

}

#undef CLOSURE
#undef RSMAX
#undef VPDMINFACTOR

double penman(double tair,
	      double elevation,
	      double rad,
	      double vpd,
	      double ra,
	      double rc,
	      double rarc)
{
  double evap;                  /* Penman-Monteith evapotranspiration */
  double slope;                 /* slope of saturated vapor pressure curve */
  double r_air;                 /* density of air in kg/m3 */
  double h;                     /* scale height in the atmosphere (m) */
  double lv;                    /* latent heat of vaporization (J/kg) */
  double pz;                    /* surface air pressure */
  double gamma;                 /* psychrometric constant (Pa/C) */

  /* calculate the slope of the saturated vapor pressure curve in Pa/K */
  slope = svp_slope(tair);

  /* calculate scale height based on average temperature in the column */
  h  = 287/9.81 * ((tair + 273.15) + 0.5 * (double)elevation * LAPSE_PM);

  /* use hypsometric equation to calculate p_z, assume that virtual
     temperature is equal air_temp */
  pz = PS_PM * exp(-(double)elevation/h);
  
  /* calculate latent heat of vaporization. Eq. 4.2.1 in Handbook of 
     Hydrology, assume Ts is Tair */
  lv = 2501000 - 2361 * tair;
  
  /* calculate gamma. Eq. 4.2.28. Handbook of Hydrology */
  gamma = 1628.6 * pz/lv;
  
  /* calculate factor to be applied to rc/ra */
  
  /* calculate the air density, using eq. 4.2.4 Handbook of Hydrology */
  r_air = 0.003486 * pz/(275 + tair);
 
  /* calculate the evaporation in mm/day (by not dividing by the density 
     of water (~1000 kg/m3)), the result ends up being in mm instead of m */ 
    
  evap = (slope * rad + r_air * CP_PM * vpd/ra)/
         (lv * (slope + gamma * (1 + (rc + rarc)/ra))) * SEC_PER_DAY;

  if (vpd >= 0.0 && evap < 0.0) 
    evap = 0.0;

  return evap;
}

