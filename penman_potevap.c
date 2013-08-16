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
  changes   : 07-May-04 Changed
				if (vpd > 0.0 && evap < 0.0)
			to
				if (vpd >= 0.0 && evap < 0.0)
			to correctly handle evap when vpd == 0.0.	TJB
  references: 
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: penman.c,v 4.1.2.1 2004/05/10 18:43:37 tbohn Exp $";

#define CLOSURE 4000		/** Pa **/
#define RSMAX 5000
#define VPDMINFACTOR 0.1

double penman_potevap(double rad, 
		      double vpd, 
		      double ra, 
		      double rs, 
		      double rarc, 
		      double lai, 
		      double gsm_inv, 
		      double tair, 
		      double net_short, 
		      float  elevation, 
		      float  RGL)
{
  double evap;                  /* Penman-Monteith evapotranspiration */
  double slope;                 /* slope of saturated vapor pressure curve */
  double rc;                    /* canopy resistance */
  double r_air;                 /* density of air in kg/m3 */
  double h;                     /* scale height in the atmosphere (m) */
  double lv;                    /* latent heat of vaporization (J/kg) */
  double pz;                    /* surface air pressure */
  double gamma;                 /* psychrometric constant (Pa/C) */
  double Tfactor;               /* factor for canopy resistance based on temperature */
  double vpdfactor;             /* factor for canopy resistance based on vpd */
  double DAYfactor;             /* factor for canopy resistance based on photosynthesis */
  double f;

  /* calculate the slope of the saturated vapor pressure curve in Pa/K */
  slope = svp_slope(tair)*1000;

  /* calculate resistance factors (rs is 0, so it doesn't really influence anything... (ALMA potevap definitions) */
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

  //  printf("1 penman_potevap evap=%f slope=%.2f rad=%.2f r_air=%.2f CP_PM=%d vpd=%.2f ra=%.2f lv=%.2f gamma=%.2f \n",
  // evap,slope,rad,r_air,CP_PM,vpd,ra,lv,gamma);

  if (vpd >= 0.0 && evap < 0.0) 
    evap = 0.0;

  return evap;
}

