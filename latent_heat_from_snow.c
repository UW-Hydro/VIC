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
			   double *VaporMassFlux,
			   double *BlowingMassFlux,
			   double *SurfaceMassFlux,
			   double Dt,
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
			   int Nveg,
			   int iveg) {
/**********************************************************************
  latent_heat_from_snow.c       Laura Bowling           

  Split out of the snowpack energy balance, this subroutine computes
  the latent heat from the snowpack.

  Modifications:
  11-18-02 Modified to handle the effects of blowing snow.       LCB

***********************************************************************/

  extern option_struct options;
  double EsSnow;
  double Ls;

  EsSnow = svp(TMean);

  // SurfaceMassFlux and BlowingMassFlux in kg*m-2*s

  *SurfaceMassFlux = AirDens * ( EPS / Press ) * ( EactAir - EsSnow ) / Ra;
  
  if ( Vpd == 0.0 && *SurfaceMassFlux < 0.0 ) 
    *SurfaceMassFlux = 0.0;

  //      if( !overstory && options.BLOWING) {
  //     Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
  //	*BlowingMassFlux = CalcBlowingSnow(Dt, Tair, LastSnow, SurfaceLiquidWater, Wind, Ls, AirDens, Press, EactAir, Z0, Z, SnowDepth, lag_one, sigma_slope, TMean, iveg, Nveg, fetch);
  //	*BlowingMassFlux *= 3600./RHO_W;
  // }
  //    else
  //     *BlowingMassFlux = 0.0;
  
  /* Calculate latent heat flux */
   *VaporMassFlux = *SurfaceMassFlux + *BlowingMassFlux*Density/3600.;
   *SurfaceMassFlux *= 3600./Density;
   //  fprintf(stderr, "SurfaceMAssFlux = %f\n",*SurfaceMassFlux);

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

  *VaporMassFlux /= Density;

}  
