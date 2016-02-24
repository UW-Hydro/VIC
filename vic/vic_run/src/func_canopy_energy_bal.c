/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine iterates to determine the temperature of the canopy, and solve
 * the resulting fluxes between the canopy and the atmosphere and the canopy
 * and the ground.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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
 * @brief    Calculate the canopy energy balance.
 *****************************************************************************/
double
func_canopy_energy_bal(double  Tfoliage,
                       va_list ap)
{
    extern option_struct     options;
    extern parameters_struct param;

    /* General Model Parameters */
    double                   delta_t;
    double                   elevation;

    double                  *Wmax;
    double                  *Wcr;
    double                  *Wpwp;
    double                  *frost_fract;

    /* Atmopheric Condition and Forcings */
    double                   AirDens;
    double                   EactAir;
    double                   Press;
    double                   Le;
    double                   Tcanopy;
    double                   Vpd;
    double                   shortwave;
    double                   Catm;
    double                  *dryFrac;

    double                  *Evap;
    double                  *Ra;
    double                  *Ra_used;
    double                   Rainfall;
    double                  *Wind;

    /* Vegetation Terms */
    int                      veg_class;

    double                  *displacement;
    double                  *ref_height;
    double                  *roughness;

    double                  *root;
    double                  *CanopLayerBnd;

    /* Water Flux Terms */
    double                   IntRain;
    double                   IntSnow;

    double                  *Wdew;

    layer_data_struct       *layer;
    veg_var_struct          *veg_var;

    /* Energy Flux Terms */
    double                   LongOverIn;
    double                   LongUnderOut;
    double                   NetShortOver;

    double                  *AdvectedEnergy;
    double                  *LatentHeat;
    double                  *LatentHeatSub;
    double                  *LongOverOut;
    double                  *NetLongOver;
    double                  *NetRadiation;
    double                  *RefreezeEnergy;
    double                  *SensibleHeat;
    double                  *VaporMassFlux;

    /* Internal Variables */
    double                   EsSnow;
    double                   Ls;
    double                   RestTerm;
    double                   prec;

    /** Read variables from variable length argument list **/

    /* General Model Parameters */
    delta_t = (double) va_arg(ap, double);
    elevation = (double) va_arg(ap, double);

    Wmax = (double *) va_arg(ap, double *);
    Wcr = (double *) va_arg(ap, double *);
    Wpwp = (double *) va_arg(ap, double *);
    frost_fract = (double *) va_arg(ap, double *);

    /* Atmopheric Condition and Forcings */
    AirDens = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Le = (double) va_arg(ap, double);
    Tcanopy = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    shortwave = (double) va_arg(ap, double);
    Catm = (double) va_arg(ap, double);
    dryFrac = (double *) va_arg(ap, double *);

    Evap = (double *) va_arg(ap, double *);
    Ra = (double *) va_arg(ap, double *);
    Ra_used = (double *) va_arg(ap, double *);
    Rainfall = (double) va_arg(ap, double);
    Wind = (double *) va_arg(ap, double *);

    /* Vegetation Terms */
    veg_class = (unsigned int) va_arg(ap, unsigned int);

    displacement = (double *) va_arg(ap, double *);
    ref_height = (double *) va_arg(ap, double *);
    roughness = (double *) va_arg(ap, double *);

    root = (double *) va_arg(ap, double *);
    CanopLayerBnd = (double *) va_arg(ap, double *);

    /* Water Flux Terms */
    IntRain = (double) va_arg(ap, double);
    IntSnow = (double) va_arg(ap, double);

    Wdew = (double *) va_arg(ap, double *);

    layer = (layer_data_struct *) va_arg(ap, layer_data_struct *);
    veg_var = (veg_var_struct *) va_arg(ap, veg_var_struct *);

    /* Energy Flux Terms */
    LongOverIn = (double) va_arg(ap, double);
    LongUnderOut = (double) va_arg(ap, double);
    NetShortOver = (double) va_arg(ap, double);

    AdvectedEnergy = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    LongOverOut = (double *) va_arg(ap, double *);
    NetLongOver = (double *) va_arg(ap, double *);
    NetRadiation = (double *) va_arg(ap, double *);
    RefreezeEnergy = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    VaporMassFlux = (double *) va_arg(ap, double *);

    /* Calculate the net radiation at the canopy surface, using the canopy
       temperature.  The outgoing longwave is subtracted twice, because the
       canopy radiates in two directions */

    *LongOverOut = calc_outgoing_longwave(Tfoliage + CONST_TKFRZ,
                                          param.EMISS_VEG);
    *NetRadiation = NetShortOver + LongOverIn + LongUnderOut -
                    2 * (*LongOverOut);

    *NetLongOver = LongOverIn - (*LongOverOut);

    if (IntSnow > 0) {
        Ra_used[0] = Ra[0];
        Ra_used[1] = Ra[1];

        /** Added multiplication by 10 to incorporate change in canopy resistance due
            to smoothing by intercepted snow **/
        if (options.AERO_RESIST_CANSNOW == AR_406 ||
            options.AERO_RESIST_CANSNOW == AR_406_LS ||
            options.AERO_RESIST_CANSNOW == AR_406_FULL) {
            Ra_used[1] *= 10.;
        }

        /** Calculate the vapor mass flux between intercepted snow in
            the canopy and the surrounding air mass **/

        EsSnow = svp(Tfoliage);

        /* Apply stability correction to aerodynamic resistance */
        if (options.AERO_RESIST_CANSNOW == AR_410) {
            if (Wind[1] > 0.0) {
                Ra_used[1] /= StabilityCorrection(ref_height[1],
                                                  displacement[1], Tfoliage,
                                                  Tcanopy, Wind[1],
                                                  roughness[1]);
            }
            else {
                Ra_used[1] = param.HUGE_RESIST;
            }
        }

        *VaporMassFlux = AirDens * (CONST_EPS / Press) * (EactAir - EsSnow) /
                         Ra_used[1] / CONST_RHOFW;

        if (Vpd == 0.0 && *VaporMassFlux < 0.0) {
            *VaporMassFlux = 0.0;
        }

        /* Calculate the latent heat flux */

        Ls = calc_latent_heat_of_sublimation(Tfoliage);
        *LatentHeatSub = Ls * *VaporMassFlux * CONST_RHOFW;
        *LatentHeat = 0;
        *Evap = 0;
        veg_var->throughfall = 0;

        if (options.AERO_RESIST_CANSNOW == AR_406) {
            Ra_used[1] /= 10;
        }
    }
    else {
        if (options.AERO_RESIST_CANSNOW == AR_406_FULL ||
            options.AERO_RESIST_CANSNOW == AR_410) {
            Ra_used[0] = Ra[0];
            Ra_used[1] = Ra[1];
        }
        else {
            Ra_used[0] = Ra[0];
            Ra_used[1] = Ra[0];
        }

        *Wdew = IntRain * MM_PER_M;
        prec = Rainfall * MM_PER_M;
        *Evap = canopy_evap(layer, veg_var, false,
                            veg_class, Wdew, delta_t, *NetRadiation,
                            Vpd, NetShortOver, Tcanopy, Ra_used[1],
                            elevation, prec, Wmax, Wcr, Wpwp, frost_fract,
                            root, dryFrac, shortwave, Catm, CanopLayerBnd);
        *Wdew /= MM_PER_M;

        *LatentHeat = Le * *Evap * CONST_RHOFW;
        *LatentHeatSub = 0;
    }

    /* Calculate the sensible heat flux */

    *SensibleHeat = calc_sensible_heat(AirDens, Tcanopy, Tfoliage, Ra_used[1]);

    /* Calculate the advected energy */

    *AdvectedEnergy = (CONST_CPFW * Tcanopy * Rainfall) / (delta_t);

    /* Calculate the amount of energy available for refreezing */

    RestTerm = *SensibleHeat + *LatentHeat + *LatentHeatSub + *NetRadiation +
               *AdvectedEnergy;

    if (IntSnow > 0) {
        /* Intercepted snow present, check if excess energy can be used to
           melt or refreeze it */

        *RefreezeEnergy = (IntRain * CONST_LATICE * CONST_RHOFW) / (delta_t);

        if (Tfoliage == 0.0 && RestTerm > -(*RefreezeEnergy)) {
            *RefreezeEnergy = -RestTerm; /* available energy input over cold content
                                            used to melt, i.e. Qrf is negative value
                                            (energy out of pack)*/
            RestTerm = 0.0;
        }
        else {
            RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
        }
    }
    else {
        *RefreezeEnergy = 0;
    }

    return (RestTerm);
}
