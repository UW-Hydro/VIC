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
  2013-Jul-25 Added calc_rc_ps(), which computes canopy resistance
	      based on photosynthetic demand.				TJB
  2014-May-05 Moved pre-compile constants CLOSURE, RSMAX, and
	      VPDMINFACTOR to vicNl_def.h.				TJB

********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: penman.c,v 5.3.2.1 2009/06/09 09:54:11 vicadmin Exp $";

double calc_rc(double rs, 
	       double net_short, 
	       float  RGL, 
	       double tair, 
	       double vpd, 
	       double lai, 
	       double gsm_inv, 
	       char   ref_crop,
	       int flag_irr)
{
  extern option_struct options;

  double f;
  double DAYfactor;             /* factor for canopy resistance based on photosynthesis */
  double Tfactor;               /* factor for canopy resistance based on temperature */
  double vpdfactor;             /* factor for canopy resistance based on vpd */
  double rc;

  if (rs == 0) {
    rc = 0;
  }
  else if (lai == 0) {
    rc = RSMAX;
  }
  else if (ref_crop) {
    /* calculate simple reference canopy resistance in s/m */
    rc = rs/(lai * 0.5);
    rc = (rc > RSMAX) ? RSMAX : rc;
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

    /*Set limitations to 1 if irrigated vegetation (both when irrig=true and when irrig=false). 
      Not perfect..., but is included because of FAO method used. Meaning T, vpd and Sw do not limit ET 
      for irrigated crops, to be consistent with FAO method used when preprocessing data. */
    if(flag_irr==1) {
      Tfactor=1.;
      DAYfactor=1.;
      vpdfactor=1.;
      if(options.IRRIGATION) gsm_inv=1.; //only when irrig=true
    }

    /* calculate canopy resistance in s/m */
    rc = rs/(lai * gsm_inv * Tfactor * vpdfactor) * DAYfactor;
    rc = (rc > RSMAX) ? RSMAX : rc;
  }

  return rc;

}


void calc_rc_ps(char    Ctype,
		double  MaxCarboxRate,
		double  MaxETransport,
		double  CO2Specificity,
		double *NscaleFactor,
		double  tair, 
		double  shortwave, 
		double *aPAR,
		double  elevation, 
		double  Catm,
		double *CanopLayerBnd,
		double  lai, 
		double  gsm_inv, 
		double  vpd, 
		double *rsLayer,
		double *rc)
/********************************************************************************
  function  : calc_rc_ps
  purpose   : Calculate canopy resistance rc as a function of photosynthetic
              activity, soil moisture stress, and vapor pressure deficit
  interface : - input : - C3/C4 type (W/m2)
                        - Maximum Carboxlyation Rate 
                        - Maximum Electron Transport Rate 
			- CO2 Specificity
			- Nitrogen scaling factor
			- air temperature
			- net shortwave radiation
			- absorbed photosynthetically active radiation
			- elevation (m)
			- Atmospheric CO2 concentration
			- Array of canopy layer boundaries (fractions of total LAI)
			- Whole-canopy LAI (m**2/m**2)
			- Soil moisture limitation factor gsm_inv
			- Vapor Pressure Deficit (Pa)
              - output: - Canopy-layer-aggregate stomatal resistance (s/m)
                        - Whole-canopy aggregate canopy resistance (s/m)
  programmer: Ted Bohn
  date      : March 7, 2007
  changes   :

  references: 
********************************************************************************/
{
  extern option_struct options;
  double GPP0;                  /* aggregate canopy assimilation (photosynthesis)
                                   in absence of soil moisture stress */
  double Rdark0;                /* aggregate canopy dark respiration in absence of
                                   soil moisture stress */
  double Rphoto0;               /* aggregate canopy photorespiration in absence of
                                   soil moisture stress */
  double Rmaint0;               /* aggregate plant maintenance respiration in absence of
                                   soil moisture stress */
  double Rgrowth0;              /* aggregate plant growth respiration in absence of
                                   soil moisture stress */
  double Raut0;                 /* aggregate plant respiration in absence of
                                   soil moisture stress */
  double NPP0;                  /* aggregate net primary productivity in absence of
                                   soil moisture stress */
  double Ci0;                   /* aggregate canopy leaf-internal CO2 mixing ratio
                                   in absence of soil moisture stress */
  double rc0;                   /* aggregate canopy resistance in absence of
                                   soil moisture stress */
  double rcRatio;
  int    cidx;
  double vpdfactor;             /* factor for canopy resistance based on vpd */

  /* Compute canopy resistance and photosynthetic demand in absence of soil moisture stress */
  canopy_assimilation(Ctype,
                      MaxCarboxRate,
                      MaxETransport,
                      CO2Specificity,
                      NscaleFactor,
                      tair,
                      shortwave,
                      aPAR,
                      elevation,
                      Catm,
                      CanopLayerBnd,
                      lai,
                      "ci",
                      rsLayer,
                      &rc0,
                      &Ci0,
                      &GPP0,
                      &Rdark0,
                      &Rphoto0,
                      &Rmaint0,
                      &Rgrowth0,
                      &Raut0,
                      &NPP0);

  /* calculate vapor pressure deficit factor */
  vpdfactor = 1 - vpd/CLOSURE;
  vpdfactor = (vpdfactor < VPDMINFACTOR) ? VPDMINFACTOR : vpdfactor;

  /* calculate canopy resistance in presence of soil moisture stress */
  *rc = rc0/(gsm_inv*vpdfactor);
  *rc = (*rc > RSMAX) ? RSMAX : *rc;
  rcRatio = *rc/rc0;
  /* this next calculation assumes canopy layers are of equal size */
  for (cidx=0; cidx<options.Ncanopy; cidx++) {
    rsLayer[cidx] *= rcRatio;
    rsLayer[cidx] = (rsLayer[cidx] > RSMAX) ? RSMAX : rsLayer[cidx];
  }

}

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
  h  = 287/9.81 * ((tair + 273.15) + 0.5 * (double)elevation * T_LAPSE);

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

