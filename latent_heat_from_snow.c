#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

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
			   double *VaporMassFlux) {

  double EsSnow;
  double Ls;

  EsSnow = svp(TMean);
  
  *VaporMassFlux = AirDens * ( EPS / Press ) * ( EactAir - EsSnow ) / Ra;
  if ( Vpd == 0.0 && *VaporMassFlux < 0.0 ) 
    *VaporMassFlux = 0.0;
  
  /* Calculate latent heat flux */
 
  if ( TMean >= 0.0 ) {
    /* Melt conditions: use latent heat of vaporization */
    *LatentHeat = Lv * *VaporMassFlux;
    *LatentHeatSublimation = 0;
  }
  else {
    /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    *LatentHeatSublimation = Ls * *VaporMassFlux;
    *LatentHeat = 0;
  }
  *VaporMassFlux /= Density;

}  
