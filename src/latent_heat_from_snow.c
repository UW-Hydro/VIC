#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void latent_heat_from_snow(double  AirDens,
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
			   double *SurfaceMassFlux) {
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
  05-Aug-04 Moved the calculation of BlowingMassFlux back out of this
	    function into surface_fluxes(), as part of merge with
	    Laura Bowling's latest code.			TJB

***********************************************************************/

  double EsSnow;
  double Ls;

  EsSnow = svp(TMean);

  // SurfaceMassFlux and BlowingMassFlux in kg/m2s

  *SurfaceMassFlux = AirDens * ( EPS / Press ) * ( EactAir - EsSnow ) / Ra;
  
  if ( Vpd == 0.0 && *SurfaceMassFlux < 0.0 ) 
    *SurfaceMassFlux = 0.0;

  /* Calculate total latent heat flux */

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
