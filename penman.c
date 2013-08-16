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

double penman(double rad, 
	      double vpd, 
	      double ra, 
	      double rs, 
	      double rarc, 
              double lai, 
	      double gsm_inv, 
	      double tair, 
              double net_short, 
	      float  elevation, 
	      float  RGL,
	      int    flag_irr)
{
  extern option_struct options;

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

  //   printf("penman slope %.3f rad %.3f rair %.3f vpd %.3f ra %.3f rc %.3f teller %.3f nevner %.3f\n",
  // slope,rad,r_air,vpd,ra,rc,slope * rad + r_air * CP_PM * vpd/ra*SEC_PER_DAY,lv * (slope + gamma * (1 + (rc + rarc)/ra))); 
  
  //if(gsm_inv<0.8 || Tfactor<0.8 || vpdfactor<0.8) printf("penman evap %f rc %f lai %f gsm_inv %f Tfactor %f vpdfactor %f DAYfactor %f\n",
  //evap,rc,lai,gsm_inv,Tfactor,vpdfactor,1/DAYfactor);
  //if(lai>5.36) printf("penman evap %f rc %f lai %f gsm_inv %f Tfactor %f vpdfactor %f DAYfactor %f\n",
  //evap,rc,lai,gsm_inv,Tfactor,vpdfactor,1/DAYfactor);


  return evap;
}

double quick_penman(double rad, 
		    double vpd, 
		    double gsm_inv,
		    double CONST_1,
		    double CONST_2,
		    double CONST_3,
		    double CONST_4,
		    double CONST_5)
{
  double evap;			/* Penman-Monteith evapotranspiration */
  double rc;			/* canopy resistance */

  /* calculate canopy resistance in s/m */
  rc = CONST_5 / gsm_inv;
  rc = (rc > RSMAX) ? RSMAX : rc;

  /* calculate the evaporation in mm/day (by not dividing by the density 
     of water (~1000 kg/m3)), the result ends up being in mm instead of m */ 
    
  evap = (CONST_1 * rad + CONST_2) / (CONST_3 + CONST_4 * rc);

/*   if (vpd > 0.0 && evap < 0.0)  */
/*     evap = 0.0; */

  return evap;
}

void compute_penman_constants(double  vpd, 
			      double  ra, 
			      double  rs, 
			      double  rarc, 
			      double  lai, 
			      double  tair, 
			      double  net_short, 
			      float   elevation, 
			      float   RGL,
			      double *CONST_1,
			      double *CONST_2,
			      double *CONST_3,
			      double *CONST_4,
			      double *CONST_5)
{
  double slope;			/* slope of saturated vapor pressure curve */
  double r_air;			/* density of air in kg/m3 */
  double h;			/* scale height in the atmosphere (m) */
  double lv;			/* latent heat of vaporization (J/kg) */
  double pz;			/* surface air pressure */
  double gamma;			/* psychrometric constant (Pa/C) */
  double Tfactor;		/* factor for canopy resistance based on 
				   temperature */
  double vpdfactor;		/* factor for canopy resistance based on vpd */
  double DAYfactor;		/* factor for canopy resistance based on 
				   photosynthesis */
  double f;

  /* calculate the slope of the saturated vapor pressure curve in Pa/K */
  slope = svp_slope(tair)*1000;

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

  /* Set constant values for later use */
  *CONST_1 = slope;
  *CONST_2 = r_air * CP_PM * vpd / ra;
  *CONST_3 = lv * (slope + gamma * (1 + rarc / ra));
  *CONST_4 = lv * gamma / ra;
  *CONST_5 = rs / (lai * Tfactor * vpdfactor) * DAYfactor;
  
}
