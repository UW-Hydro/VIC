/******************************************************************************
 * @section DESCRIPTION
 *
 * Stability correction
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate atmospheric stability correction for non-neutral
 *           conditions
 *****************************************************************************/
double
StabilityCorrection(double Z,
                    double d,
                    double TSurf,
                    double Tair,
                    double Wind,
                    double Z0)
{
    double Correction;          /* Correction to aerodynamic resistance */
    double Ri;                   /* Richardson's Number */
    double RiCr = 0.2;           /* Critical Richardson's Number */
    double RiLimit;              /* Upper limit for Richardson's Number */

    Correction = 1.0;

    /* Calculate the effect of the atmospheric stability using a Richardson
       Number approach */

    if (TSurf != Tair) {
        /* Non-neutral conditions */

        Ri = CONST_G * (Tair - TSurf) * (Z - d) /
             (((Tair +
                CONST_TKFRZ) + (TSurf + CONST_TKFRZ)) / 2.0 * Wind * Wind);

        RiLimit = (Tair + CONST_TKFRZ) /
                  (((Tair +
                     CONST_TKFRZ) +
                    (TSurf + CONST_TKFRZ)) / 2.0 * (log((Z - d) / Z0) + 5));

        if (Ri > RiLimit) {
            Ri = RiLimit;
        }

        if (Ri > 0.0) {
            Correction = (1 - Ri / RiCr) * (1 - Ri / RiCr);
        }
        else {
            if (Ri < -0.5) {
                Ri = -0.5;
            }

            Correction = sqrt(1 - 16 * Ri);
        }
    }

    return Correction;
}
