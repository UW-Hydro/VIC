/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate snow accumulation and melt for the lake model
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
 * @brief    Calculate snow accumulation and melt using an energy balance
 *           approach for a two layer snow model
 *****************************************************************************/
int
ice_melt(double            z2,
         double            aero_resist,
         double           *aero_resist_used,
         double            Le,
         snow_data_struct *snow,
         lake_var_struct  *lake,
         double            delta_t,
         double            displacement,
         double            Z0,
         double            surf_atten,
         double            rainfall,
         double            snowfall,
         double            wind,
         double            Tcutoff,
         double            air_temp,
         double            net_short,
         double            longwave,
         double            density,
         double            pressure,
         double            vpd,
         double            vp,
         double           *melt,
         double           *save_advection,
         double           *save_deltaCC,
         double           *save_SnowFlux,
         double           *save_latent,
         double           *save_sensible,
         double           *save_Qnet,
         double           *save_refreeze_energy,
         double           *save_LWnet)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   DeltaPackCC; /* Change in cold content of the pack */
    double                   DeltaPackSwq; /* Change in snow water equivalent of the pack (m) */
    double                   InitialSwq; /* Initial snow water equivalent (m) */
    double                   InitialIce;
    double                   MassBalanceError; /* Mass balance error (m) */
    double                   MaxLiquidWater; /* Maximum liquid water content of pack (m) */
    double                   OldTSurf; /* Old snow surface temperature (C) */
    double                   Qnet; /* Net energy exchange at the surface (W/m2) */
    double                   PackRefreezeEnergy; /* refreeze/melt energy in pack layer (W/m2) */
    double                   RefreezeEnergy; /* refreeze energy (W/m2) */
    double                   RefrozenWater; /* Amount of refrozen water (m) */
    double                   SnowFallCC; /* Cold content of new snowfall (J) */
    double                   SurfaceCC;
    double                   PackCC;
    double                   SurfaceSwq;
    double                   PackSwq;
    double                   PackIce;
    double                   SnowMelt; /* Amount of snow melt during time interval (m water equivalent) */
    double                   IceMelt;
    double                   LWnet;
    double                   avgcond;
    double                   SWconducted;
    double                   SnowIce;
    double                   LakeIce;
    double                   Ice;
    double                   SnowFall;
    double                   RainFall;
    double                   vapor_flux;
    double                   blowing_flux;
    double                   surface_flux;
    double                   advection;
    double                   deltaCC;
    double                   SnowFlux; /* thermal flux through snowpack from ground */
    double                   latent_heat;
    double                   latent_heat_sub;
    double                   sensible_heat;
    double                   Ls;
    double                   melt_energy = 0.;

    SnowFall = snowfall / MM_PER_M; /* convert to m */
    RainFall = rainfall / MM_PER_M; /* convert to m */
    IceMelt = 0.0;
    RefrozenWater = 0.0;

    InitialSwq = snow->swq;
    OldTSurf = snow->surf_temp;

    /* Initialize snowpack variables */
    SnowIce = snow->swq - snow->pack_water - snow->surf_water;
    LakeIce = lake->ice_water_eq / lake->areai;      /* meters of water equivalent based on average ice thickness. */
    InitialIce = LakeIce;
    Ice = SnowIce + LakeIce;

    /* Reconstruct snow pack */
    if (Ice > param.SNOW_MAX_SURFACE_SWE) {
        SurfaceSwq = param.SNOW_MAX_SURFACE_SWE;
    }
    else {
        SurfaceSwq = Ice;
    }
    if (SurfaceSwq <= SnowIce) {
        PackSwq = SnowIce - SurfaceSwq;
        PackIce = LakeIce;
    }
    else {
        PackSwq = 0.;
        PackIce = Ice - SurfaceSwq;
    }

    /* Calculate cold contents */
    SurfaceCC = CONST_VCPICE_WQ * SurfaceSwq * snow->surf_temp;
    PackCC = CONST_VCPICE_WQ * (PackSwq + PackIce) * snow->pack_temp;
    if (air_temp > 0.0) {
        SnowFallCC = 0.0;
    }
    else {
        SnowFallCC = CONST_VCPICE_WQ * SnowFall * air_temp;
    }

    /* Distribute fresh snowfall */

    /* Surface layer was not already full, snow will exceed space.
       What happens to snow if SurfaceSwq = MAX_SURFACE_SWE?*/
    if (SnowFall > (param.SNOW_MAX_SURFACE_SWE - SurfaceSwq)) {
        DeltaPackSwq = SurfaceSwq + SnowFall - param.SNOW_MAX_SURFACE_SWE;
        if (DeltaPackSwq > SurfaceSwq) {
            DeltaPackCC = SurfaceCC +
                          (SnowFall -
                           param.SNOW_MAX_SURFACE_SWE) / SnowFall * SnowFallCC;
        }
        else {
            DeltaPackCC = DeltaPackSwq / SurfaceSwq * SurfaceCC;
        }
        SurfaceSwq = param.SNOW_MAX_SURFACE_SWE;
        SurfaceCC += SnowFallCC - DeltaPackCC;
        PackSwq += DeltaPackSwq;
        PackCC += DeltaPackCC;
    }
    else {
        SurfaceSwq += SnowFall;
        SurfaceCC += SnowFallCC;
        DeltaPackCC = 0;
    }
    if (SurfaceSwq > 0.0) {
        snow->surf_temp = SurfaceCC / (CONST_VCPICE_WQ * SurfaceSwq);
    }
    else {
        snow->surf_temp = 0.0;
    }
    if (PackSwq + PackIce > 0.0) {
        snow->pack_temp = PackCC /
                          (CONST_VCPICE_WQ * (PackSwq + PackIce));
    }
    else {
        snow->pack_temp = 0.0;
    }

    /* Adjust ice and snow->surf_water */
    SnowIce += SnowFall;
    Ice += SnowFall;
    snow->surf_water += RainFall;

    icerad(net_short, lake->hice, SnowIce * CONST_RHOFW / param.LAKE_RHOSNOW,
           &avgcond,
           &SWconducted, &deltaCC);

    /* Calculate blowing snow sublimation (m/timestep) */
    // Currently have hard-wired parameters that are approximated for ice-covered area:
    // lag-one autocorrelation = 0.95, sigma_slope = .005 (both appropriate for
    // flat terrain. Fetch = 2000 m (i.e. unlimited fetch), roughness and displacement
    // calculated assuming 10 cm high protrusions on frozen ponds.

    if (options.BLOWING && snow->swq > 0.) {
        Ls = calc_latent_heat_of_sublimation(snow->surf_temp);
        snow->blowing_flux = CalcBlowingSnow(delta_t, air_temp,
                                             snow->last_snow, snow->surf_water,
                                             wind, Ls, density,
                                             vp, Z0,
                                             z2, snow->depth, .95, 0.005,
                                             snow->surf_temp, 0, 1, 100.,
                                             .067, .0123, &snow->transport);
        if ((int)snow->blowing_flux == ERROR) {
            log_err("Error calculating blowing snow flux");
        }

        snow->blowing_flux *= delta_t / CONST_RHOFW;
    }
    else {
        snow->blowing_flux = 0.0;
    }

    /* Store sublimation terms in temporary variables */
    vapor_flux = snow->vapor_flux;
    blowing_flux = snow->blowing_flux;
    surface_flux = snow->surface_flux;

    /* Calculate the surface energy balance for snow_temp = 0.0 */

    Qnet = CalcIcePackEnergyBalance((double) 0.0, delta_t, aero_resist,
                                    aero_resist_used, z2, Z0, wind, net_short,
                                    longwave, density, Le, air_temp,
                                    pressure * PA_PER_KPA, vpd * PA_PER_KPA,
                                    vp * PA_PER_KPA,
                                    RainFall, snow->surf_water, &RefreezeEnergy,
                                    &vapor_flux, &blowing_flux, &surface_flux,
                                    &advection, Tcutoff, avgcond, SWconducted,
                                    &SnowFlux, &latent_heat, &latent_heat_sub,
                                    &sensible_heat, &LWnet);

    snow->vapor_flux = vapor_flux;
    snow->surface_flux = surface_flux;
    save_refreeze_energy[0] = RefreezeEnergy;

    /* Check that snow swq exceeds minimum value for model stability */

    /* If Qnet == 0.0, then set the surface temperature to 0.0 */
    if (Qnet == 0.0) {
        snow->surf_temp = 0.0;
        if (RefreezeEnergy >= 0.0) {                 /* Surface is freezing. */
            RefrozenWater = RefreezeEnergy /
                            (CONST_LATICE *
                             CONST_RHOFW) * delta_t;
            if (RefrozenWater > snow->surf_water) {
                RefrozenWater = snow->surf_water;
                RefreezeEnergy = RefrozenWater * CONST_LATICE * CONST_RHOFW /
                                 (delta_t);
            }
            melt_energy += RefreezeEnergy;
            SurfaceSwq += RefrozenWater;
            SnowIce += RefrozenWater;
            Ice += RefrozenWater;
            snow->surf_water -= RefrozenWater;
            if (snow->surf_water < 0.0) {
                snow->surf_water = 0.0;
            }
            SnowMelt = 0.0;
        }
        else {
            /* Calculate snow melt if refreeze energy is negative */
            SnowMelt = fabs(RefreezeEnergy) /
                       (CONST_LATICE * CONST_RHOFW) * delta_t;
            melt_energy += RefreezeEnergy;
        }

        /* Adjust snow->surf_water for vapor_flux */
        if (snow->surf_water < -(snow->vapor_flux)) {
            // if vapor_flux exceeds stored water, we not only need to
            // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
            snow->blowing_flux *= -(snow->surf_water) / snow->vapor_flux;
            snow->vapor_flux = -(snow->surf_water);
            snow->surface_flux = -(snow->surf_water) - snow->blowing_flux;
            snow->surf_water = 0.0;
        }
        else {
            snow->surf_water += snow->vapor_flux;
        }

        /* If SnowMelt < Ice, there was incomplete melting of the snow/ice */
        if (SnowMelt < Ice) {
            /* Subtract melt from snow pack and lake ice */

            /* Since we redefine the snow surface layer to always encompass the topmost
               portion of the snow, we essentially reduce the bottom snow layer first.
               This is only for accounting purposes and may not reflect where in the
               pack the melt actually occurs. */
            if (SnowMelt <= PackSwq) { /* Only melt part of the pack (bottom) layer. */
                snow->surf_water += SnowMelt;
                PackSwq -= SnowMelt;
                Ice -= SnowMelt;
                SnowIce -= SnowMelt;
            }
            else if (SnowMelt <= SnowIce) { /* Melt all of pack layer and part of surface layer. */
                snow->surf_water += SnowMelt + snow->pack_water;
                snow->pack_water = 0.0;
                SurfaceSwq -= (SnowMelt - PackSwq);
                PackSwq = 0.0;
                SnowIce -= SnowMelt;
                Ice -= SnowMelt;
            }
            else { /* Melt snow pack completely and also part of the ice */
                snow->surf_water += SnowIce + snow->pack_water;
                snow->pack_water = 0.0;
                PackSwq = 0.0;
                Ice -= SnowMelt;
                LakeIce -= SnowMelt - SnowIce;
                IceMelt = SnowMelt - SnowIce;
                if (SurfaceSwq > SnowMelt) {
                    SurfaceSwq -= SnowMelt;
                }
                else {
                    SurfaceSwq = 0.0;
                    PackIce -= (SnowMelt - SurfaceSwq - PackSwq);
                }
                SnowIce = 0.0;
            }
        }
        /* Else, SnowMelt > TotalIce and there was complete melting of the snow and ice */
        else {
            snow->surf_water += SnowIce + snow->pack_water;
            snow->pack_water = 0.0;
            PackSwq = 0.0;
            SurfaceSwq = 0.0;
            SnowIce = 0.0;
            SnowMelt = Ice;
            IceMelt = LakeIce;
            LakeIce = 0.0;
            PackIce = 0.0;
            Ice = 0.0;
            snow->surf_temp = 0.0;
            snow->pack_temp = 0.0;
            /* readjust melt energy to account for melt only of available snow */
            melt_energy -= RefreezeEnergy;
            RefreezeEnergy = RefreezeEnergy / fabs(RefreezeEnergy) * SnowMelt *
                             CONST_LATICE * CONST_RHOFW / (delta_t);
            melt_energy += RefreezeEnergy;
        }
    }
    /* Else, IceEnergyBalance(T=0.0) <= 0.0 */
    else {
        /* Calculate surface layer temperature using "Brent method" */
        if (SurfaceSwq > param.SNOW_MIN_SWQ_EB_THRES) {
            snow->surf_temp =
                root_brent((double) (snow->surf_temp - param.SNOW_DT),
                           (double) (snow->surf_temp + param.SNOW_DT),
                           IceEnergyBalance, delta_t,
                           aero_resist, aero_resist_used, z2, Z0,
                           wind, net_short, longwave, density,
                           Le, air_temp, pressure * PA_PER_KPA,
                           vpd * PA_PER_KPA, vp * PA_PER_KPA,
                           RainFall,
                           snow->surf_water, &RefreezeEnergy,
                           &vapor_flux, &blowing_flux,
                           &surface_flux, &advection, Tcutoff,
                           avgcond, SWconducted, &SnowFlux,
                           &latent_heat, &latent_heat_sub,
                           &sensible_heat, &LWnet);

            if (snow->surf_temp <= -998) {
                if (options.TFALLBACK) {
                    snow->surf_temp = OldTSurf;
                    snow->surf_temp_fbflag = 1;
                    snow->surf_temp_fbcount++;
                }
                else {
                    ErrorIcePackEnergyBalance(snow->surf_temp, delta_t,
                                              aero_resist,
                                              aero_resist_used, z2,
                                              displacement, Z0, wind, net_short,
                                              longwave, density, Le, air_temp,
                                              pressure * PA_PER_KPA,
                                              vpd * PA_PER_KPA,
                                              vp * PA_PER_KPA,
                                              RainFall, SurfaceSwq,
                                              snow->surf_water, OldTSurf,
                                              &RefreezeEnergy,
                                              &vapor_flux, &blowing_flux,
                                              &surface_flux,
                                              &advection, deltaCC, Tcutoff,
                                              avgcond, SWconducted,
                                              snow->swq * CONST_RHOFW / param.LAKE_RHOSNOW,
                                              param.LAKE_RHOSNOW, surf_atten,
                                              &SnowFlux, &latent_heat,
                                              &latent_heat_sub,
                                              &sensible_heat, &LWnet);
                    return(ERROR);
                }
            }
        }
        else {
            snow->surf_temp = 999;
        }
        if (snow->surf_temp > -998 && snow->surf_temp < 999) {
            Qnet = CalcIcePackEnergyBalance(snow->surf_temp, delta_t,
                                            aero_resist, aero_resist_used, z2,
                                            Z0, wind, net_short, longwave,
                                            density, Le, air_temp,
                                            pressure * PA_PER_KPA,
                                            vpd * PA_PER_KPA,
                                            vp * PA_PER_KPA, RainFall,
                                            snow->surf_water, &RefreezeEnergy,
                                            &vapor_flux, &blowing_flux,
                                            &surface_flux, &advection, Tcutoff,
                                            avgcond, SWconducted, &SnowFlux,
                                            &latent_heat, &latent_heat_sub,
                                            &sensible_heat, &LWnet);

            snow->vapor_flux = vapor_flux;
            snow->surface_flux = surface_flux;
            save_refreeze_energy[0] = RefreezeEnergy;

            /* since we iterated, the surface layer is below freezing and no snowmelt */

            SnowMelt = 0.0;
            IceMelt = 0.0;

            /* Since updated snow_temp < 0.0, all of the liquid water in the surface
               layer has been frozen */

            SnowIce += snow->surf_water;
            Ice += snow->surf_water;
            melt_energy += snow->surf_water * CONST_LATICE * CONST_RHOFW /
                           (delta_t);
            RefrozenWater = snow->surf_water;
            snow->surf_water = 0.0;

            /* Adjust SurfaceSwq for vapor_flux */
            if (SurfaceSwq < -(snow->vapor_flux)) {
                // if vapor_flux exceeds stored snow/ice, we not only need to
                // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
                if (SurfaceSwq > SnowIce) {
                    snow->blowing_flux *= -(SurfaceSwq) / snow->vapor_flux;
                    snow->vapor_flux = -SurfaceSwq;
                    snow->surface_flux = -SurfaceSwq - snow->blowing_flux;
                    LakeIce -= SurfaceSwq - SnowIce;
                    Ice = PackIce;
                    SnowIce = 0.0;
                }
                else {
                    snow->blowing_flux *= -(SurfaceSwq) / snow->vapor_flux;
                    snow->vapor_flux = -SurfaceSwq;
                    snow->surface_flux = -SurfaceSwq - snow->blowing_flux;
                    SurfaceSwq = 0.0;
                    Ice = PackSwq + PackIce;
                }
            }
            else {
                SurfaceSwq += snow->vapor_flux;
                if (SnowIce > -(snow->vapor_flux)) {
                    SnowIce += snow->vapor_flux;
                }
                else {
                    LakeIce += (snow->vapor_flux + SnowIce);
                    SnowIce = 0.;
                }
                Ice += snow->vapor_flux;
            }
        }
        else {
            snow->surf_temp = 999;
        }
    }

    /* Done with iteration etc, now Update the liquid water content of the
       surface layer */

    if (SnowIce > SurfaceSwq) {
        MaxLiquidWater = param.SNOW_LIQUID_WATER_CAPACITY * SurfaceSwq;
    }
    else {
        MaxLiquidWater = param.SNOW_LIQUID_WATER_CAPACITY * SnowIce;
    }
    if (snow->surf_water > MaxLiquidWater) {
        melt[0] = snow->surf_water - MaxLiquidWater;
        snow->surf_water = MaxLiquidWater;
    }
    else {
        melt[0] = 0.0;
    }

    /* Refreeze liquid water in the pack.
       variable 'RefreezeEnergy' is the heat released to the snow pack
       if all liquid water were refrozen.
       if RefreezeEnergy < PackCC then all water IS refrozen
       PackCC always <=0.0
       WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does
       not involve energy transported to the pixel.  Instead heat from the snow
       pack is used to refreeze water */
    snow->pack_water += melt[0]; /* add surface layer outflow to pack liquid water*/
    PackRefreezeEnergy = snow->pack_water * CONST_LATICE * CONST_RHOFW;

    /* calculate energy released to freeze*/
    if (PackCC < -PackRefreezeEnergy) { /* cold content not fully depleted*/
        PackSwq += snow->pack_water;        /* refreeze all water and update*/
        Ice += snow->pack_water;
        SnowIce += snow->pack_water;
        snow->pack_water = 0.0;
        if (PackSwq + PackIce > 0.0) {
            PackCC =
                (PackSwq +
                 PackIce) * CONST_VCPICE_WQ * snow->pack_temp +
                PackRefreezeEnergy;
            snow->pack_temp = PackCC / (CONST_VCPICE_WQ *
                                        (PackSwq + PackIce));
            if (snow->pack_temp > 0.) {
                snow->pack_temp = 0.;
            }
        }
        else {
            snow->pack_temp = 0.0;
        }
    }
    else {
        /* cold content has been either exactly satisfied or exceeded. If
           PackCC = refreeze then pack is ripe and all pack water is
           refrozen, else if energy released in refreezing exceeds PackCC
           then exactly the right amount of water is refrozen to satify PackCC.
           The refrozen water is added to PackSwq and Ice */
        snow->pack_temp = 0.0;
        DeltaPackSwq = -PackCC / (CONST_LATICE * CONST_RHOFW);
        snow->pack_water -= DeltaPackSwq;
        PackSwq += DeltaPackSwq;
        Ice += DeltaPackSwq;
        SnowIce += DeltaPackSwq;
    }

    /* Update the liquid water content of the pack */
    MaxLiquidWater = param.SNOW_LIQUID_WATER_CAPACITY * PackSwq;
    if (snow->pack_water > MaxLiquidWater) {
        melt[0] = snow->pack_water - MaxLiquidWater;
        snow->pack_water = MaxLiquidWater;
    }
    else {
        melt[0] = 0.0;
    }

    /* Update snow properties */
    Ice = PackIce + PackSwq + SurfaceSwq;
    if (Ice > param.SNOW_MAX_SURFACE_SWE) {
        SurfaceCC = CONST_VCPICE_WQ * snow->surf_temp * SurfaceSwq;
        PackCC = CONST_VCPICE_WQ * snow->pack_temp * (PackSwq + PackIce);
        if (SurfaceSwq > param.SNOW_MAX_SURFACE_SWE) {
            PackCC += SurfaceCC *
                      (SurfaceSwq - param.SNOW_MAX_SURFACE_SWE) / SurfaceSwq;
            SurfaceCC -= SurfaceCC *
                         (SurfaceSwq - param.SNOW_MAX_SURFACE_SWE) / SurfaceSwq;
            PackSwq += SurfaceSwq - param.SNOW_MAX_SURFACE_SWE;
            SurfaceSwq -= SurfaceSwq - param.SNOW_MAX_SURFACE_SWE;
        }
        else if (SurfaceSwq < param.SNOW_MAX_SURFACE_SWE) {
            PackCC -= PackCC *
                      (param.SNOW_MAX_SURFACE_SWE -
                       SurfaceSwq) / (PackSwq + PackIce);
            SurfaceCC += PackCC *
                         (param.SNOW_MAX_SURFACE_SWE -
                          SurfaceSwq) / (PackSwq + PackIce);
            PackSwq -= param.SNOW_MAX_SURFACE_SWE - SurfaceSwq;
            SurfaceSwq += param.SNOW_MAX_SURFACE_SWE - SurfaceSwq;
        }
        snow->pack_temp = PackCC / (CONST_VCPICE_WQ *
                                    (PackSwq + PackIce));
        snow->surf_temp = SurfaceCC / (CONST_VCPICE_WQ * SurfaceSwq);
    }
    else {
        PackSwq = 0.0;
        PackCC = 0.0;
        PackIce = 0.0;
        snow->pack_temp = 0.0;
    }
    snow->swq = SnowIce + snow->surf_water + snow->pack_water;
    lake->ice_water_eq = LakeIce * lake->areai;
    lake->volume -= (InitialIce - LakeIce - IceMelt) * lake->areai;
    if (lake->ice_water_eq <= 0.0) {
        lake->ice_water_eq = 0.0;
    }
    if (!options.SPATIAL_SNOW) {
        if (snow->swq > 0) {
            snow->coverage = 1.;
        }
        else {
            snow->coverage = 0.;
        }
    }

    /* Mass balance test */
    MassBalanceError = (InitialSwq - snow->swq) + (InitialIce - LakeIce) +
                       (RainFall +
                        SnowFall) - IceMelt - melt[0] + snow->vapor_flux;

    melt[0] *= MM_PER_M; /* converts back to mm */
    snow->mass_error = MassBalanceError;
    snow->coldcontent = SurfaceCC;
    snow->vapor_flux *= -1.;
    *save_LWnet = LWnet;
    *save_advection = advection;
    *save_deltaCC = deltaCC;
    *save_SnowFlux = SnowFlux;
    *save_latent = latent_heat + latent_heat_sub;
    *save_sensible = sensible_heat;
    *save_refreeze_energy = RefreezeEnergy;
    *save_Qnet = Qnet;

    return (0);
}

/******************************************************************************
 * @brief    Dummy function to make a direct call to IceEnergyBalance()
 *           possible.
 *****************************************************************************/
double
CalcIcePackEnergyBalance(double Tsurf,
                         ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    double  Qnet;                /* Net energy exchange at the IcePack snow
                                    surface (W/m^2) */

    va_start(ap, Tsurf);
    Qnet = IceEnergyBalance(Tsurf, ap);
    va_end(ap);

    return Qnet;
}

/******************************************************************************
 * @brief    Dummy function to make a direct call to
 *           ErrorPrintIcePackEnergyBalance() possible.
 *****************************************************************************/
double
ErrorIcePackEnergyBalance(double Tsurf,
                          ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    double  Qnet;                /* Net energy exchange at the IcePack snow
                                    surface (W/m^2) */

    va_start(ap, Tsurf);
    Qnet = ErrorPrintIcePackEnergyBalance(Tsurf, ap);
    va_end(ap);

    return Qnet;
}

/******************************************************************************
 * @brief    Print ice pack energy balance terms
 *****************************************************************************/
double
ErrorPrintIcePackEnergyBalance(double  TSurf,
                               va_list ap)
{
    double  Dt;                  /* Model time step (seconds) */
    double  Ra;                  /* Aerodynamic resistance (s/m) */
    double *Ra_used;             /* Aerodynamic resistance (s/m) after stability correction */
    double  Z;                   /* Reference height (m) */
    double  Displacement;        /* Displacement height (m) */
    double  Z0;                  /* surface roughness height (m) */
    double  Wind;                /* Wind speed (m/s) */
    double  ShortRad;            /* Net incident shortwave radiation (W/m2) */
    double  LongRadIn;           /* Incoming longwave radiation (W/m2) */
    double  AirDens;             /* Density of air (kg/m3) */
    double  Lv;                  /* Latent heat of vaporization (J/kg3) */
    double  Tair;                /* Air temperature (C) */
    double  Press;               /* Air pressure (Pa) */
    double  Vpd;                /* Vapor pressure deficit (Pa) */
    double  EactAir;             /* Actual vapor pressure of air (Pa) */
    double  Rain;                /* Rain fall (m/timestep) */
    double  SweSurfaceLayer;     /* Snow water equivalent in surface layer (m)
                                  */
    double  SurfaceLiquidWater;  /* Liquid water in the surface layer (m) */
    double  OldTSurf;            /* Surface temperature during previous time
                                    step */
    double *RefreezeEnergy;      /* Refreeze energy (W/m2) */
    double *vapor_flux;          /* Total mass flux of water vapor to or from
                                    snow (m/timestep) */
    double *blowing_flux;        /* Mass flux of water vapor to or from
                                    blowing snow (m/timestep) */
    double *surface_flux;        /* Mass flux of water vapor to or from
                                    snow pack (m/timestep) */
    double *AdvectedEnergy;      /* Energy advected by precipitation (W/m2) */
    double  DeltaColdContent;    /* Change in cold content (W/m2) */
    double  Tfreeze;
    double  AvgCond;
    double  SWconducted;
    double  SnowDepth;
    double  SnowDensity;
    double  SurfAttenuation;

    /* end of list of arguments in variable argument list */
    double *GroundFlux;
    double *LatentHeat;         /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;      /* Latent heat exchange at surface (W/m2) due to sublimation */
    double *SensibleHeat;       /* Sensible heat exchange at surface (W/m2) */
    double *LWnet;

    /* initialize variables */
    Dt = (double) va_arg(ap, double);
    Ra = (double) va_arg(ap, double);
    Ra_used = (double *) va_arg(ap, double *);
    Z = (double) va_arg(ap, double);
    Displacement = (double) va_arg(ap, double);
    Z0 = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);
    ShortRad = (double) va_arg(ap, double);
    LongRadIn = (double) va_arg(ap, double);
    AirDens = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Tair = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    SweSurfaceLayer = (double) va_arg(ap, double);
    SurfaceLiquidWater = (double) va_arg(ap, double);
    OldTSurf = (double) va_arg(ap, double);
    RefreezeEnergy = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);
    blowing_flux = (double *) va_arg(ap, double *);
    surface_flux = (double *) va_arg(ap, double *);
    AdvectedEnergy = (double *) va_arg(ap, double *);
    DeltaColdContent = (double) va_arg(ap, double);
    Tfreeze = (double) va_arg(ap, double);
    AvgCond = (double) va_arg(ap, double);
    SWconducted = (double) va_arg(ap, double);
    SnowDepth = (double) va_arg(ap, double);
    SnowDensity = (double) va_arg(ap, double);
    SurfAttenuation = (double) va_arg(ap, double);
    GroundFlux = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    LWnet = (double *) va_arg(ap, double *);

    /* print variables */
    log_warn("ice_melt failed to converge to a solution in root_brent.  "
             "Variable values will be dumped to the screen, check for invalid "
             "values.");

    fprintf(LOG_DEST, "Dt = %f\n", Dt);
    fprintf(LOG_DEST, "Ra = %f\n", Ra);
    fprintf(LOG_DEST, "Ra_used = %f\n", *Ra_used);
    fprintf(LOG_DEST, "Z = %f\n", Z);
    fprintf(LOG_DEST, "Displacement = %f\n", Displacement);
    fprintf(LOG_DEST, "Z0 = %f\n", Z0);
    fprintf(LOG_DEST, "Wind = %f\n", Wind);
    fprintf(LOG_DEST, "ShortRad = %f\n", ShortRad);
    fprintf(LOG_DEST, "LongRadIn = %f\n", LongRadIn);
    fprintf(LOG_DEST, "AirDens = %f\n", AirDens);
    fprintf(LOG_DEST, "Lv = %f\n", Lv);
    fprintf(LOG_DEST, "Tair = %f\n", Tair);
    fprintf(LOG_DEST, "Press = %f\n", Press);
    fprintf(LOG_DEST, "Vpd = %f\n", Vpd);
    fprintf(LOG_DEST, "EactAir = %f\n", EactAir);
    fprintf(LOG_DEST, "Rain = %f\n", Rain);
    fprintf(LOG_DEST, "SweSurfaceLayer = %f\n", SweSurfaceLayer);
    fprintf(LOG_DEST, "SurfaceLiquidWater = %f\n", SurfaceLiquidWater);
    fprintf(LOG_DEST, "TSurf = %f\n", TSurf);
    fprintf(LOG_DEST, "OldTSurf = %f\n", OldTSurf);
    fprintf(LOG_DEST, "RefreezeEnergy = %f\n", RefreezeEnergy[0]);
    fprintf(LOG_DEST, "vapor_flux = %f\n", *vapor_flux);
    fprintf(LOG_DEST, "blowing_flux = %f\n", *blowing_flux);
    fprintf(LOG_DEST, "surface_flux = %f\n", *surface_flux);
    fprintf(LOG_DEST, "AdvectedEnergy = %f\n", AdvectedEnergy[0]);
    fprintf(LOG_DEST, "DeltaColdContent = %f\n", DeltaColdContent);
    fprintf(LOG_DEST, "Tfreeze = %f\n", Tfreeze);
    fprintf(LOG_DEST, "AvgCond = %f\n", AvgCond);
    fprintf(LOG_DEST, "SWconducted = %f\n", SWconducted);
    fprintf(LOG_DEST, "SnowDepth = %f\n", SnowDepth);
    fprintf(LOG_DEST, "SnowDensity = %f\n", SnowDensity);
    fprintf(LOG_DEST, "SurfAttenuation = %f\n", SurfAttenuation);
    fprintf(LOG_DEST, "GroundFlux = %f\n", GroundFlux[0]);
    fprintf(LOG_DEST, "LatentHeat = %f\n", LatentHeat[0]);
    fprintf(LOG_DEST, "LatentHeatSub = %f\n", LatentHeatSub[0]);
    fprintf(LOG_DEST, "SensibleHeat = %f\n", SensibleHeat[0]);
    fprintf(LOG_DEST, "LWnet = %f\n", *LWnet);

    log_warn("Finished dumping snow_melt variables.\n"
             "Try increasing SNOW_DT to get model to complete cell.\n"
             "Then check output for instabilities.");

    return(ERROR);
}
