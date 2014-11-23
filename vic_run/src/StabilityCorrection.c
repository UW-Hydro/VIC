/*
 * SUMMARY:      StabilityCorrection.c - Calculate the stability correction
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD: Fri Mar  2 18:42:20 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate the stability correction for exchange of sensible
 *               heat between the surface and the atmosphere
 * DESCRIP-END.
 * FUNCTIONS:    StabilityCorrection()
 * COMMENTS:
 */

#include <vic_def.h>
#include <vic_run.h>

/*****************************************************************************
   Function name: StabilityCorrection()

   Purpose      : Calculate atmospheric stability correction for non-neutral
                 conditions

   Required     :
    double Z          - Reference height (m)
    double d          - Displacement height (m)
    double TSurf      - Surface temperature (C)
    double Tair       - Air temperature (C)
    double Wind       - Wind speed (m/s)
    double Z0         - Roughness length (m)

   Returns      :
    double Correction - Multiplier for aerodynamic resistance

   Modifies     : None

   Comments     :
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

        Ri = G * (Tair - TSurf) * (Z - d) /
             (((Tair + 273.15) + (TSurf + 273.15)) / 2.0 * Wind * Wind);

        RiLimit = (Tair + 273.15) /
                  (((Tair +
                     273.15) +
                    (TSurf + 273.15)) / 2.0 * (log((Z - d) / Z0) + 5));

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
