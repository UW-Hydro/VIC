/******************************************************************************
 * @section DESCRIPTION
 *
 * calculates nitrogen scaling factors for all canopy layers, following eqns
 * 106 and 107 in Knorr 1997.
 *
 * Note: this should only be applied to veg types that have a canopy, e.g.
 * trees and shrubs, but not grass or tundra vegetation.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate nitrogen scaling factors for all canopy layers.
 *****************************************************************************/
void
calc_Nscale_factors(bool    NscaleFlag,
                    double *CanopLayerBnd,
                    double  LAItotal,
                    double  coszen,
                    double *NscaleFactor)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   k12;
    size_t                   cidx; // canopy layer index

    /* solar zenith angle at local noon */
    if (coszen < param.PHOTO_ZENITHMINPAR) {
        coszen = param.PHOTO_ZENITHMINPAR;
    }

    /* Extinction factor; eqn 119c in Knorr 1997 */
    k12 = 0.5 / coszen;

    /* Condition: LAI > LaiLimit; eqns 107 and 108 in Knorr 1997 */
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
        if (NscaleFlag && LAItotal > param.PHOTO_LAILIMIT && cidx > 0) {
            NscaleFactor[cidx] = exp(-k12 * CanopLayerBnd[cidx - 1] * LAItotal);
            if (NscaleFactor[cidx] < 1e-10) {
                NscaleFactor[cidx] = 1e-10;
            }
        }
        else {
            NscaleFactor[cidx] = 1;
        }
    }
}
