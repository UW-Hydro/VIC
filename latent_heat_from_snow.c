#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void latent_heat_from_snow(double  AirDens,
			   double  Density,
			   double  EactAir,
			   double  Lv,
			   double  Press,
			   double  Ra,
			   double  TMean, 
			   double  Vpd,
			   double *LatentHeat,
			   double *LatentHeatSublimation,
			   double *VaporMassFlux,
			   double *BlowingMassFlux,
			   double *SurfaceMassFlux,
			   double delta_t,
			   double Tair,
			   int LastSnow,
			   double SurfaceLiquidWater,
			   double Wind,
			   double *Z0,
			   double Z,
			   double SnowDepth,
			   int overstory,
			   float lag_one,
			   float sigma_slope,
			   float fetch,
			   int iveg,
			   int Nveg,
			   int month) {
/**********************************************************************
  latent_heat_from_snow.c       Laura Bowling           

  Split out of the snowpack energy balance, this subroutine computes
  the latent heat from the snowpack.

  Modifications:
  11-18-02 Modified to handle the effects of blowing snow.       LCB
  16-Jul-04 Moved the calculation of BlowingMassFlux back into this
	    function.  Added "month" to the parameter list and added
	    declaration of *veg_lib, to enable this calculation to
	    occur.  Modified calculations of all sublimation terms
	    to ensure that VaporMassFlux, BlowingMassFlux, and
	    SurfaceMassFlux consistently have units of kg/m2s.	TJB

***********************************************************************/

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  double EsSnow;
  double Ls;

  EsSnow = svp(TMean);

  // SurfaceMassFlux and BlowingMassFlux in kg/m2s

  *SurfaceMassFlux = AirDens * ( EPS / Press ) * ( EactAir - EsSnow ) / Ra;
  
  if ( Vpd == 0.0 && *SurfaceMassFlux < 0.0 ) 
    *SurfaceMassFlux = 0.0;

  if( options.BLOWING && !overstory ) {
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    *BlowingMassFlux = CalcBlowingSnow(delta_t/SECPHOUR, Tair, LastSnow,
				       SurfaceLiquidWater, Wind, Ls,
				       AirDens, Press, EactAir, Z0, Z,
				       SnowDepth, lag_one, sigma_slope,
				       TMean, iveg, Nveg, fetch,
				       veg_lib[iveg].displacement[month],
				       veg_lib[iveg].roughness[month]);
  }
  else
    *BlowingMassFlux = 0.0;
  
  /* Calculate latent heat flux */
  *VaporMassFlux = *SurfaceMassFlux + *BlowingMassFlux;

  if ( TMean >= 0.0 ) {
    /* Melt conditions: use latent heat of vaporization */
    *LatentHeat = Lv * (*VaporMassFlux);
    *LatentHeatSublimation = 0;
  }
  else {
    /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    *LatentHeatSublimation = Ls * (*VaporMassFlux);
    *LatentHeat = 0;
  }

}  
