/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate GPP, Raut, and NPP for veg cover with multi-layer canopy.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Calculate GPP, Raut, and NPP for veg cover with multi-layer canopy
 *****************************************************************************/
void
canopy_assimilation(char    Ctype,
                    double  MaxCarboxRate,
                    double  MaxETransport,
                    double  CO2Specificity,
                    double *NscaleFactor,
                    double  Tfoliage,
                    double  SWdown,
                    double *aPAR,
                    double  elevation,
                    double  Catm,
                    double *CanopLayerBnd,
                    double  LAItotal,
                    char   *mode,
                    double *rsLayer,
                    double *rc,
                    double *Ci,
                    double *GPP,
                    double *Rdark,
                    double *Rphoto,
                    double *Rmaint,
                    double *Rgrowth,
                    double *Raut,
                    double *NPP)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   h;
    double                   pz;
    size_t                   cidx;
    double                   dLAI;
    double                  *CiLayer = NULL;
    double                   AgrossLayer;
    double                   RdarkLayer;
    double                   RphotoLayer;
    double                   gc; /* 1/rs */

    /* calculate scale height based on average temperature in the column */
    h = calc_scale_height(Tfoliage, elevation);

    /* use hypsometric equation to calculate p_z, assume that virtual
       temperature is equal air_temp */
    pz = CONST_PSTD * exp(-(double) elevation / h);

    CiLayer = calloc(options.Ncanopy, sizeof(*CiLayer));
    check_alloc_status(CiLayer, "Memory allocation error.");

    if (!strcasecmp(mode, "ci")) {
        /* Assume a default leaf-internal CO2; compute assimilation,
           respiration, and stomatal resistance */

        /* Default leaf-internal CO2 */
        for (cidx = 0; cidx < options.Ncanopy; cidx++) {
            if (Ctype == PHOTO_C3) {
                CiLayer[cidx] = param.PHOTO_FCI1C3 * Catm;
            }
            else if (Ctype == PHOTO_C4) {
                CiLayer[cidx] = param.PHOTO_FCI1C4 * Catm;
            }
        }
        if (Ctype == PHOTO_C3) {
            *Ci = param.PHOTO_FCI1C3 * Catm;
        }
        else if (Ctype == PHOTO_C4) {
            *Ci = param.PHOTO_FCI1C4 * Catm;
        }

        /* Sum over canopy layers */
        *GPP = 0.0;
        *Rdark = 0.0;
        *Rphoto = 0.0;
        gc = 0.0;
        for (cidx = 0; cidx < options.Ncanopy; cidx++) {
            photosynth(Ctype,
                       MaxCarboxRate,
                       MaxETransport,
                       CO2Specificity,
                       NscaleFactor[cidx],
                       Tfoliage,
                       SWdown / param.PHOTO_EPAR, /* note: divide by Epar to convert from W/m2 to mol(photons)/m2s */
                       aPAR[cidx],
                       pz,
                       Catm,
                       mode,
                       &(rsLayer[cidx]),
                       &(CiLayer[cidx]),
                       &RdarkLayer,
                       &RphotoLayer,
                       &AgrossLayer);

            if (cidx > 0) {
                dLAI = LAItotal *
                       (CanopLayerBnd[cidx] - CanopLayerBnd[cidx - 1]);
            }
            else {
                dLAI = LAItotal * CanopLayerBnd[cidx];
            }

            *GPP += AgrossLayer * dLAI;
            *Rdark += RdarkLayer * dLAI;
            *Rphoto += RphotoLayer * dLAI;
            gc += (1 / rsLayer[cidx]) * dLAI;
        }

        if (gc < DBL_EPSILON) {
            gc = DBL_EPSILON;
        }
        *rc = 1 / gc;
        if (*rc > param.HUGE_RESIST) {
            *rc = param.HUGE_RESIST;
        }
    }
    else {
        /* Stomatal resistance given; compute assimilation, respiration, and leaf-internal CO2 */

        /* Sum over canopy layers */
        *GPP = 0.0;
        *Rdark = 0.0;
        *Rphoto = 0.0;
        *Ci = 0.0;
        for (cidx = 0; cidx < options.Ncanopy; cidx++) {
            photosynth(Ctype,
                       MaxCarboxRate,
                       MaxETransport,
                       CO2Specificity,
                       NscaleFactor[cidx],
                       Tfoliage,
                       SWdown / param.PHOTO_EPAR,
                       aPAR[cidx],
                       pz,
                       Catm,
                       mode,
                       &(rsLayer[cidx]),
                       &(CiLayer[cidx]),
                       &RdarkLayer,
                       &RphotoLayer,
                       &AgrossLayer);

            if (cidx > 0) {
                dLAI = LAItotal *
                       (CanopLayerBnd[cidx] - CanopLayerBnd[cidx - 1]);
            }
            else {
                dLAI = LAItotal * CanopLayerBnd[cidx];
            }

            *GPP += AgrossLayer * dLAI;
            *Rdark += RdarkLayer * dLAI;
            *Rphoto += RphotoLayer * dLAI;
            *Ci += CiLayer[cidx] * dLAI;
        }
    }

    /* Compute whole-plant respiration terms and NPP */
    *Rmaint = *Rdark / param.PHOTO_FRLEAF;
    *Rgrowth = (param.PHOTO_FRGROWTH / (1 + param.PHOTO_FRGROWTH)) *
               ((*GPP) - (*Rmaint));
    *Raut = *Rmaint + *Rgrowth;
    *NPP = *GPP - *Raut;

    free((char*) CiLayer);
}
