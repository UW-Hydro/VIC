/******************************************************************************
 * @section DESCRIPTION
 *
 * Computes the latent heat from the snowpack.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Computes the latent heat from the snowpack.
 *****************************************************************************/
void
latent_heat_from_snow(double  AirDens,
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
                      double *SurfaceMassFlux)
{
    double EsSnow;
    double Ls;

    EsSnow = svp(TMean);

    // SurfaceMassFlux and BlowingMassFlux in kg/m2s

    *SurfaceMassFlux = AirDens * (CONST_EPS / Press) * (EactAir - EsSnow) / Ra;

    if (Vpd == 0.0 && *SurfaceMassFlux < 0.0) {
        *SurfaceMassFlux = 0.0;
    }

    /* Calculate total latent heat flux */

    *VaporMassFlux = *SurfaceMassFlux + *BlowingMassFlux;

    if (TMean >= 0.0) {
        /* Melt conditions: use latent heat of vaporization */
        *LatentHeat = Lv * (*VaporMassFlux);
        *LatentHeatSublimation = 0;
    }
    else {
        /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
        Ls = calc_latent_heat_of_sublimation(TMean);
        *LatentHeatSublimation = Ls * (*VaporMassFlux);
        *LatentHeat = 0;
    }
}
