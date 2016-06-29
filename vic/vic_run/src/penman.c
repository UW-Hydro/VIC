/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine calculates evapotranspiration using the Penman-Monteith
 * approach.
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
 * @brief    Calculate canopy resistance
 *****************************************************************************/
double
calc_rc(double rs,
        double net_short,
        double RGL,
        double tair,
        double vpd,
        double lai,
        double gsm_inv,
        char   ref_crop)
{
    extern parameters_struct param;

    double                   f;
    double                   DAYfactor; /* factor for canopy resistance based on photosynthesis */
    double                   Tfactor; /* factor for canopy resistance based on temperature */
    double                   vpdfactor; /* factor for canopy resistance based on vpd */
    double                   rc;

    if (rs == 0) {
        rc = 0;
    }
    else if (lai == 0) {
        rc = param.CANOPY_RSMAX;
    }
    else if (ref_crop) {
        /* calculate simple reference canopy resistance in s/m */
        rc = rs / (lai * 0.5);
        rc = (rc > param.CANOPY_RSMAX) ? param.CANOPY_RSMAX : rc;
    }
    else {
        /* calculate resistance factors (Wigmosta et al., 1994) */
        if (rs > 0.) {
            f = net_short / RGL;
            DAYfactor = (1. + f) / (f + rs / param.CANOPY_RSMAX);
        }
        else {
            DAYfactor = 1.;
        }

        Tfactor = .08 * tair - 0.0016 * tair * tair;
        Tfactor = (Tfactor <= 0.0) ? 1e-10 : Tfactor;

        vpdfactor = 1 - vpd / param.CANOPY_CLOSURE;
        vpdfactor =
            (vpdfactor <
             param.CANOPY_VPDMINFACTOR) ? param.CANOPY_VPDMINFACTOR : vpdfactor;

        /* calculate canopy resistance in s/m */
        rc = rs / (lai * gsm_inv * Tfactor * vpdfactor) * DAYfactor;
        rc = (rc > param.CANOPY_RSMAX) ? param.CANOPY_RSMAX : rc;
    }

    return rc;
}

/******************************************************************************
 * @brief    Calculate canopy resistance rc as a function of photosynthetic.
 *****************************************************************************/
void
calc_rc_ps(char    Ctype,
           double  MaxCarboxRate,
           double  MaxETransport,
           double  CO2Specificity,
           double *NscaleFactor,
           double  tair,
           double  shortwave,
           double *aPAR,
           double  elevation,
           double  Catm,
           double *CanopLayerBnd,
           double  lai,
           double  gsm_inv,
           double  vpd,
           double *rsLayer,
           double *rc)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   GPP0; /* aggregate canopy assimilation (photosynthesis)
                                      in absence of soil moisture stress */
    double                   Rdark0; /* aggregate canopy dark respiration in absence of
                                        soil moisture stress */
    double                   Rphoto0; /* aggregate canopy photorespiration in absence of
                                         soil moisture stress */
    double                   Rmaint0; /* aggregate plant maintenance respiration in absence of
                                         soil moisture stress */
    double                   Rgrowth0; /* aggregate plant growth respiration in absence of
                                          soil moisture stress */
    double                   Raut0; /* aggregate plant respiration in absence of
                                       soil moisture stress */
    double                   NPP0; /* aggregate net primary productivity in absence of
                                      soil moisture stress */
    double                   Ci0; /* aggregate canopy leaf-internal CO2 mixing ratio
                                     in absence of soil moisture stress */
    double                   rc0; /* aggregate canopy resistance in absence of
                                     soil moisture stress */
    double                   rcRatio;
    double                   vpdfactor; /* factor for canopy resistance based on vpd */
    size_t                   cidx;

    /* Compute canopy resistance and photosynthetic demand in absence of soil moisture stress */
    canopy_assimilation(Ctype,
                        MaxCarboxRate,
                        MaxETransport,
                        CO2Specificity,
                        NscaleFactor,
                        tair,
                        shortwave,
                        aPAR,
                        elevation,
                        Catm,
                        CanopLayerBnd,
                        lai,
                        "ci",
                        rsLayer,
                        &rc0,
                        &Ci0,
                        &GPP0,
                        &Rdark0,
                        &Rphoto0,
                        &Rmaint0,
                        &Rgrowth0,
                        &Raut0,
                        &NPP0);

    /* calculate vapor pressure deficit factor */
    vpdfactor = 1 - vpd / param.CANOPY_CLOSURE;
    vpdfactor =
        (vpdfactor <
         param.CANOPY_VPDMINFACTOR) ? param.CANOPY_VPDMINFACTOR : vpdfactor;

    /* calculate canopy resistance in presence of soil moisture stress */
    *rc = rc0 / (gsm_inv * vpdfactor);
    *rc = (*rc > param.CANOPY_RSMAX) ? param.CANOPY_RSMAX : *rc;
    rcRatio = *rc / rc0;
    /* this next calculation assumes canopy layers are of equal size */
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
        rsLayer[cidx] *= rcRatio;
        rsLayer[cidx] =
            (rsLayer[cidx] >
             param.CANOPY_RSMAX) ? param.CANOPY_RSMAX : rsLayer[cidx];
    }
}

/******************************************************************************
 * @brief    Calculate daily evapotranspiration using the combination equation.
 *****************************************************************************/
double
penman(double tair,
       double elevation,
       double rad,
       double vpd,
       double ra,
       double rc,
       double rarc)
{
    double evap;                /* Penman-Monteith evapotranspiration */
    double slope;               /* slope of saturated vapor pressure curve */
    double r_air;               /* density of air in kg/m3 */
    double h;                   /* scale height in the atmosphere (m) */
    double lv;                  /* latent heat of vaporization (J/kg) */
    double pz;                  /* surface air pressure */
    double gamma;               /* psychrometric constant (Pa/C) */

    /* calculate the slope of the saturated vapor pressure curve in Pa/K */
    slope = svp_slope(tair);

    /* calculate scale height based on average temperature in the column */
    h = calc_scale_height(tair, elevation);

    /* use hypsometric equation to calculate p_z, assume that virtual
       temperature is equal air_temp */
    pz = CONST_PSTD * exp(-elevation / h);

    /* calculate latent heat of vaporization. Eq. 4.2.1 in Handbook of
       Hydrology, assume Ts is Tair */
    lv = calc_latent_heat_of_vaporization(tair);

    /* calculate gamma. Eq. 4.2.28. Handbook of Hydrology */
    gamma = 1628.6 * pz / lv;

    /* calculate factor to be applied to rc/ra */

    /* calculate the air density, using eq. 4.2.4 Handbook of Hydrology */
    r_air = 0.003486 * pz / (275 + tair);

    /* calculate the evaporation in mm/day (by not dividing by the density
       of water (~1000 kg/m3)), the result ends up being in mm instead of m */

    evap = (slope * rad + r_air * CONST_CPMAIR * vpd / ra) /
           (lv * (slope + gamma * (1 + (rc + rarc) / ra))) * CONST_CDAY;

    if (vpd >= 0.0 && evap < 0.0) {
        evap = 0.0;
    }

    return evap;
}
