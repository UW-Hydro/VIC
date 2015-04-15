#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

int
snow_melt(double            latent,
          double            net_short_snow,       // net SW at absorbed by snow
          double            t_canopy,
          double            t_grnd,
          double           *z0,    // roughness
          double            aero_resist,       // aerodynamic resistance
          double           *aero_resist_used,        // stability-corrected aerodynamic resistance
          double            air_temp,       // air temperature
          double            coverage,  // snowpack cover fraction
          double            delta_t,   // time step in secs
          double            density,       // atmospheric density
          double            displacement,   // surface displacement
          double            grnd_flux,       // ground heat flux
          double            long_snow_in,   // incoming longwave radiation
          double            pressure,
          double            rainfall,
          double            snowfall,
          double            vp,
          double            vpd,
          double            wind,
          double            z2,
          double           *net_long_snow,
          double           *old_tsurf,
          double           *melt,
          double           *save_qnet,
          double           *save_advected_sensible,
          double           *save_advection,
          double           *save_deltaCC,
          double           *save_grnd_flux,
          double           *save_latent,
          double           *save_latent_sub,
          double           *save_refreeze_energy,
          double           *save_sensible,
          int               unstable_snow,
          int               rec,
          int               iveg,
          int               band,
          snow_data_struct *snow,
          soil_con_struct  *soil_con)
{
    extern option_struct options;
    double               error;
    double               delta_pack_glacier_cc; /* Change in cold content of the pack */
    double               delta_pack_glacier_we; /* Change in snow water equivalent of the
                                                   pack (m) */
    double               surface_pack_glacier_we;    /* Ice content of snow pack (m)*/
    double               initial_surface_pack_we; /* Initial snow water equivalent (m) */
    double               initial_glacier_we; /* Initial ice water equivalent (m) */
    double               mass_balance_error; /* Mass balance error (m) */
    double               max_liquid_water; /* Maximum liquid water content of pack (m) */
    double               surface_pack_cc; /* Cold content of snow pack (J) */
    double               pack_glacier_we; /* Snow pack snow water equivalent (m) */
    double               surface_pack_glacier_depth; /* Snow pack depth (m) */
    double               surface_pack_glacier_density; /* Snow pack density (kg/m3) */
    double               qnet;   /* Net energy exchange at the surface (W/m2) */
    double               refreeze_energy; /* refreeze/melt energy in surface layer (W/m2) */
    double               pack_glacier_refreeze_energy; /* refreeze/melt energy in pack layer (W/m2) */
    double               refrozen_water; /* Amount of refrozen water (m) */
    double               snowfall_cc; /* Cold content of new snowfall (J) */
    double               snowmelt; /* Amount of snow melt during time interval
                                      (m water equivalent) */
    double               surface_cc; /* Cold content of snow pack (J) */
    double               surface_we; /* Surface layer snow water equivalent (m) */
    double               snowfall_m;
    double               rainfall_m;
    double               advection;
    double               delta_cc;
    double               latent_heat;
    double               latent_heat_sub;
    double               sensible_heat;
    double               advected_sensible_heat;
    double               melt_energy = 0.;
    double               delswe;
    double               deliwe;
    double               tmp;
    char                 errorstring[MAXSTRING];

    snowfall_m = snowfall / (double) MMPERMETER; /* convet to m */
    rainfall_m = rainfall / (double) MMPERMETER; /* convet to m */

    initial_surface_pack_we = snow->swq;
    initial_glacier_we = snow->iwq;
    snowmelt = 0.0;
    (*old_tsurf) = snow->surf_temp;

    /* Initialize snowpack variables */
    surface_pack_glacier_we = snow->swq + snow->iwq - snow->pack_water -
                              snow->surf_water;

    /* Reconstruct snow pack */
    if (surface_pack_glacier_we > MAX_SURFACE_SWE) {
        surface_we = MAX_SURFACE_SWE;
        pack_glacier_we = surface_pack_glacier_we - surface_we;
    }
    else {
        surface_we = surface_pack_glacier_we;
        pack_glacier_we = 0.;
    }

    if (snow->icedepth > 0.) {
        surface_pack_glacier_depth = snow->depth + snow->icedepth;
        surface_pack_glacier_density = (snow->density * snow->depth +
                                        ice_density *
                                        snow->icedepth) /
                                       surface_pack_glacier_depth;
    }
    else {
        surface_pack_glacier_depth = snow->depth;
        surface_pack_glacier_density = snow->density;
    }

    /* Calculate cold contents */
    surface_cc = calc_cold_content(surface_we, snow->surf_temp);
    tmp = pack_glacier_we - snow->iwq;

    if (tmp < 0.) {
        surface_pack_cc = 0.;
    }
    else {
        surface_pack_cc = calc_cold_content(tmp, snow->pack_temp);
    }

    if (air_temp > 0.0) {
        snowfall_cc = 0.0;
    }
    else {
        snowfall_cc = calc_cold_content(snowfall_m, air_temp);
    }

    /* Distribute fresh snowfall */
    if (snowfall_m > (MAX_SURFACE_SWE - surface_we) &&
        (MAX_SURFACE_SWE - surface_we) > SMALL) {
        delta_pack_glacier_we = surface_we + snowfall_m - MAX_SURFACE_SWE;
        if (delta_pack_glacier_we > surface_we) {
            delta_pack_glacier_cc = surface_cc +
                                    (snowfall_m -
                                     MAX_SURFACE_SWE) / snowfall_m *
                                    snowfall_cc;
        }
        else {
            delta_pack_glacier_cc = delta_pack_glacier_we / surface_we *
                                    surface_cc;
        }
        surface_we = MAX_SURFACE_SWE;
        surface_cc += snowfall_cc - delta_pack_glacier_cc;
        pack_glacier_we += delta_pack_glacier_we;
        surface_pack_cc += delta_pack_glacier_cc;
    }
    else {
        surface_we += snowfall_m;
        surface_cc += snowfall_cc;
    }
    if (surface_we > 0.0) {
        snow->surf_temp = calc_temp_from_cc(surface_we, surface_cc);
    }
    else {
        snow->surf_temp = 0.0;
    }
    if (pack_glacier_we > 0.0) {
        snow->pack_temp = calc_temp_from_cc(pack_glacier_we, surface_pack_cc);
    }
    else {
        snow->pack_temp = 0.0;
    }

    /* Adjust ice and snow->surf_water */
    surface_pack_glacier_we += snowfall_m;
    snow->surf_water += rainfall_m;

    /* Calculate the surface energy balance for snow_temp = 0.0 */
    qnet = CalcSnowPackEnergyBalance((double) 0.0, delta_t, aero_resist,
                                     aero_resist_used,
                                     displacement, z2, z0,
                                     density, vp, long_snow_in, latent,
                                     pressure, rainfall_m, net_short_snow, vpd,
                                     wind, (*old_tsurf), coverage,
                                     surface_pack_glacier_depth,
                                     surface_pack_glacier_density,
                                     snow->surf_water, surface_we,
                                     t_canopy, t_grnd,
                                     &advection, &advected_sensible_heat,
                                     &delta_cc,
                                     &grnd_flux, &latent_heat,
                                     &latent_heat_sub, net_long_snow,
                                     &refreeze_energy, &sensible_heat,
                                     &snow->vapor_flux, &snow->blowing_flux,
                                     &snow->surface_flux);

    // If we can solve the snow SEB seperate from the ground surface
    if (!unstable_snow) {
        /* If qnet == 0.0, then set the surface temperature to 0.0 */
        if (qnet == 0.0) {
            snow->surf_temp = 0.0;
            if (refreeze_energy >= 0.0) {
                refrozen_water = refreeze_energy / (Lf * RHO_W) * delta_t;
                if (refrozen_water > snow->surf_water) {
                    refrozen_water = snow->surf_water;
                    refreeze_energy = refrozen_water * Lf * RHO_W / (delta_t);
                }
                melt_energy += refreeze_energy;
                surface_we += refrozen_water;
                surface_pack_glacier_we += refrozen_water;
                snow->surf_water -= refrozen_water;
                if (snow->surf_water < 0.0) {
                    snow->surf_water = 0.0;
                }
                snowmelt = 0.0;
            }
            else {
                /* Calculate snow melt */
                snowmelt = fabs(refreeze_energy) / (Lf * RHO_W) * delta_t;
                melt_energy += refreeze_energy;
            }

            /* Adjust snow->surf_water for vapor_flux */
            if (snow->surf_water < -(snow->vapor_flux)) {
                // if vapor_flux exceeds surf_water, we not only need to
                // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
                snow->blowing_flux *= -(snow->surf_water / snow->vapor_flux);
                snow->vapor_flux = -(snow->surf_water);
                snow->surface_flux = -(snow->surf_water) - snow->blowing_flux;
                snow->surf_water = 0.0;
            }
            else {
                snow->surf_water += snow->vapor_flux;
            }

            /* If snowmelt < surface_pack_glacier_we, there was incomplete melting of the pack */
            if (snowmelt < surface_pack_glacier_we) {
                if (snowmelt <= pack_glacier_we) {
                    snow->surf_water += snowmelt;
                    pack_glacier_we -= snowmelt;
                    surface_pack_glacier_we -= snowmelt;
                }
                else {
                    snow->surf_water += snowmelt + snow->pack_water;
                    snow->pack_water = 0.0;
                    pack_glacier_we = 0.0;
                    surface_pack_glacier_we -= snowmelt;
                    surface_we = surface_pack_glacier_we;
                }
            }
            /* Else, snowmelt > surface_pack_glacier_we and there was complete melting of the pack */
            else {
                snowmelt = surface_pack_glacier_we;
                snow->surf_water += surface_pack_glacier_we;
                surface_we = 0.0;
                snow->surf_temp = 0.0;
                pack_glacier_we = 0.0;
                snow->pack_temp = 0.0;
                surface_pack_glacier_we = 0.0;
                /* readjust melt energy to account for melt only of available snow */
                melt_energy -= refreeze_energy;
                refreeze_energy = refreeze_energy / fabs(refreeze_energy) *
                                  snowmelt * Lf * RHO_W / (delta_t);
                melt_energy += refreeze_energy;
            }
        }
        /* Else, SnowPackEnergyBalance(T=0.0) <= 0.0 */
        else {
            /* Calculate surface layer temperature using "Brent method" */
            if (surface_we > MIN_SWQ_EB_THRES) {
                snow->surf_temp = root_brent(
                    (double) (snow->surf_temp - SNOW_DT),
                    (double) (snow->surf_temp + SNOW_DT),
                    errorstring, SnowPackEnergyBalance,
                    delta_t, aero_resist, aero_resist_used,
                    displacement, z2, z0,
                    density, vp, long_snow_in, latent, pressure,
                    rainfall_m, net_short_snow, vpd,
                    wind, (*old_tsurf), coverage,
                    surface_pack_glacier_depth, surface_pack_glacier_density,
                    snow->surf_water, surface_we,
                    t_canopy, t_grnd,
                    &advection, &advected_sensible_heat,
                    &delta_cc,
                    &grnd_flux, &latent_heat,
                    &latent_heat_sub, net_long_snow,
                    &refreeze_energy, &sensible_heat,
                    &snow->vapor_flux, &snow->blowing_flux,
                    &snow->surface_flux);

                if (snow->surf_temp <= -998) {
                    if (options.TFALLBACK) {
                        snow->surf_temp = *old_tsurf;
                        snow->surf_temp_fbflag = 1;
                        snow->surf_temp_fbcount++;
                    }
                    else {
                        error = ErrorSnowPackEnergyBalance(snow->surf_temp, rec,
                                                           iveg, band,
                                                           delta_t, aero_resist,
                                                           aero_resist_used,
                                                           displacement, z2, z0,
                                                           density, vp,
                                                           long_snow_in, latent,
                                                           pressure,
                                                           rainfall_m,
                                                           net_short_snow, vpd,
                                                           wind, (*old_tsurf),
                                                           coverage,
                                                           surface_pack_glacier_depth,
                                                           surface_pack_glacier_density,
                                                           snow->surf_water,
                                                           surface_we,
                                                           t_canopy, t_grnd,
                                                           &advection,
                                                           &advected_sensible_heat,
                                                           &delta_cc,
                                                           &grnd_flux, &latent_heat,
                                                           &latent_heat_sub,
                                                           net_long_snow, &refreeze_energy,
                                                           &sensible_heat, &snow->vapor_flux,
                                                           &snow->blowing_flux,
                                                           &snow->surface_flux);
                        return(ERROR);
                    }
                }
            }
            else {
                /* Thin snowpack must be solved in conjunction with ground
                   surface energy balance */
                snow->surf_temp = 999;
            }
            if (snow->surf_temp > -998 && snow->surf_temp < 999) {
                qnet = CalcSnowPackEnergyBalance(snow->surf_temp,
                                                 delta_t, aero_resist,
                                                 aero_resist_used,
                                                 displacement, z2, z0,
                                                 density, vp, long_snow_in,
                                                 latent,
                                                 pressure,
                                                 rainfall_m, net_short_snow,
                                                 vpd,
                                                 wind, (*old_tsurf), coverage,
                                                 surface_pack_glacier_depth,
                                                 surface_pack_glacier_density,
                                                 snow->surf_water, surface_we,
                                                 t_canopy, t_grnd,
                                                 &advection,
                                                 &advected_sensible_heat,
                                                 &delta_cc,
                                                 &grnd_flux, &latent_heat,
                                                 &latent_heat_sub,
                                                 net_long_snow,
                                                 &refreeze_energy,
                                                 &sensible_heat,
                                                 &snow->vapor_flux,
                                                 &snow->blowing_flux,
                                                 &snow->surface_flux);

                /* since we iterated, the surface layer is below freezing and no snowmelt */
                snowmelt = 0.0;

                /* Since updated snow_temp < 0.0, all of the liquid water in the surface
                   layer has been frozen */
                surface_we += snow->surf_water;
                surface_pack_glacier_we += snow->surf_water;
                snow->surf_water = 0.0;
                melt_energy += snow->surf_water * Lf * RHO_W / (delta_t);

                /* Adjust surface_we for vapor flux */
                if (surface_we < -(snow->vapor_flux)) {
                    // if vapor_flux exceeds surface_we, we not only need to
                    // re-scale vapor_flux, we need to re-scale surface_flux and blowing_flux
                    snow->blowing_flux *= -(surface_we / snow->vapor_flux);
                    snow->vapor_flux = -surface_we;
                    snow->surface_flux = -surface_we - snow->blowing_flux;
                    surface_we = 0.0;
                    surface_pack_glacier_we = pack_glacier_we;
                }
                else {
                    surface_we += snow->vapor_flux;
                    surface_pack_glacier_we += snow->vapor_flux;
                }
            }
        }
    }
    else {
        /* Snow solution is unstable as independent layer */
        snow->surf_temp = 999;
    }

    /* Done with iteration etc, now Update the liquid water content of the
       surface layer */
    max_liquid_water = LIQUID_WATER_CAPACITY * surface_we;
    if (snow->surf_water > max_liquid_water) {
        melt[0] = snow->surf_water - max_liquid_water;
        snow->surf_water = max_liquid_water;
    }
    else {
        melt[0] = 0.0;
    }

    /* Refreeze liquid water in the pack.
       variable 'refreeze_energy' is the heat released to the snow pack
       if all liquid water were refrozen.
       if refreeze_energy < surface_pack_cc then all water IS refrozen
       surface_pack_cc always <=0.0

       WORK IN PROGRESS: This energy is NOT added to MeltEnergy, since this does
       not involve energy transported to the pixel.  Instead heat from the snow
       pack is used to refreeze water */
    snow->pack_water += melt[0]; /* add surface layer outflow to pack
                                        liquid water*/
    pack_glacier_refreeze_energy = snow->pack_water * Lf * RHO_W;

    /* calculate energy released to freeze*/
    if (surface_pack_cc < -pack_glacier_refreeze_energy) { /* cold content not fully depleted*/
        pack_glacier_we += snow->pack_water;  /* refreeze all water and update*/
        surface_pack_glacier_we += snow->pack_water;
        snow->pack_water = 0.0;
        if (pack_glacier_we > 0.0) {
            surface_pack_cc += pack_glacier_refreeze_energy;
            snow->pack_temp =
                calc_temp_from_cc(pack_glacier_we, surface_pack_cc);
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
           surface_pack_cc = refreeze then pack is ripe and all pack water is
           refrozen, else if energy released in refreezing exceeds surface_pack_cc
           then exactly the right amount of water is refrozen to satify surface_pack_cc.
           The refrozen water is added to pack_glacier_we and surface_pack_glacier_we */
        snow->pack_temp = 0.0;
        delta_pack_glacier_we = -surface_pack_cc / (Lf * RHO_W);
        snow->pack_water -= delta_pack_glacier_we;
        pack_glacier_we += delta_pack_glacier_we;
        surface_pack_glacier_we += delta_pack_glacier_we;
    }

    /* Update the liquid water content of the pack */
    max_liquid_water = LIQUID_WATER_CAPACITY * pack_glacier_we;
    if (snow->pack_water > max_liquid_water) {
        melt[0] = snow->pack_water - max_liquid_water;
        snow->pack_water = max_liquid_water;
    }
    else {
        melt[0] = 0.0;
    }

    /* Update snow properties */
    surface_pack_glacier_we = pack_glacier_we + surface_we;

    // Update glacier properties and remove glacier ice from pack
    if (initial_glacier_we > 0.) {
        if (surface_pack_glacier_we >= initial_glacier_we) {
            // No melting of the glacier ice so iwq didn't change
            snow->iwq = initial_glacier_we;
            surface_pack_glacier_we -= initial_glacier_we;
            snow->glmelt = 0.0;
        }
        else {
            // Complete melting of snow and some glacier
            snow->iwq = surface_pack_glacier_we;  // What's left is glacier
            surface_pack_glacier_we = 0.0;  // No ice to be apportioned later
            snow->glmelt = initial_glacier_we - surface_pack_glacier_we;  // Lost ice
        }

        // Handle case where the vapor flux deposits directly onto the glacier.
        // If we began without snow and there was no snowfall, all the mass
        // will remain part of the glacier.
        if ((initial_surface_pack_we == 0) && (snowfall_m == 0)) {
            snow->iwq += surface_pack_glacier_we;
            surface_pack_glacier_we = 0.;
        }
    }  // After removing the glacier ice from the "surface_pack_glacier_we" pack, the rest should
       // be just like the standard VIC snow_melt() function

    if (surface_pack_glacier_we > MAX_SURFACE_SWE) {
        surface_cc = calc_cold_content(surface_we, snow->surf_temp);
        surface_pack_cc = calc_cold_content(pack_glacier_we, snow->pack_temp);
        if (surface_we > MAX_SURFACE_SWE) {
            surface_pack_cc += surface_cc *
                               (surface_we - MAX_SURFACE_SWE) / surface_we;
            surface_cc -= surface_cc *
                          (surface_we - MAX_SURFACE_SWE) / surface_we;
            pack_glacier_we += surface_we - MAX_SURFACE_SWE;
            surface_we -= surface_we - MAX_SURFACE_SWE;
        }
        else if (surface_we < MAX_SURFACE_SWE) {
            surface_pack_cc -= surface_pack_cc *
                               (MAX_SURFACE_SWE - surface_we) / pack_glacier_we;
            surface_cc += surface_pack_cc *
                          (MAX_SURFACE_SWE - surface_we) / pack_glacier_we;
            pack_glacier_we -= MAX_SURFACE_SWE - surface_we;
            surface_we += MAX_SURFACE_SWE - surface_we;
        }
        snow->pack_temp = calc_temp_from_cc(pack_glacier_we, surface_pack_cc);
        snow->surf_temp = calc_temp_from_cc(surface_we, surface_cc);
    }
    else {
        pack_glacier_we = 0.0;
        surface_pack_cc = 0.0;
        snow->pack_temp = 0.0;
    }

    if (surface_pack_glacier_we > 0.) {
        // Reconstitute the swq variable
        snow->swq = surface_pack_glacier_we + snow->pack_water +
                    snow->surf_water;
    }
    else {
        // surface_pack_glacier_we was totally depleted so remove any remaining liquid water
        snow->swq = 0.;
        melt[0] += snow->pack_water + snow->surf_water;
        snow->pack_water = 0.;
        snow->surf_water = 0.;
    }
    snow->icedepth = calc_snow_depth(snow->iwq, ice_density);

    if (snow->swq <= 0.0) {
        snow->swq = 0.0;
        snow->surf_temp = 0.0;
        snow->pack_temp = 0.0;
    }

    /* Snow solution is unstable as independent layer */
    if (snow->surf_temp == 999) {
        snow->swq += snow->iwq;  // Include ice in snowpack for now
        snow->iwq = 0.0;
        snow->depth = surface_pack_glacier_depth;
        snow->icedepth = 0.0;
        snow->density = surface_pack_glacier_density;
    }

    /* Mass balance test */
    mass_balance_error =
        (initial_surface_pack_we - snow->swq) +
        (initial_glacier_we - snow->iwq) +
        (rainfall_m + snowfall_m) - melt[0] + snow->vapor_flux;

    // glacier mass balance calculation
    delswe = snow->swq - snow->swqold; // change in swe in previous time step
    snow->swqold = snow->swq;
    deliwe = snow->iwq - snow->iwqold; // change in iwe in previous time step
    snow->iwqold = snow->iwq;
    snow->bn = delswe + deliwe; // glacier mass balance

    melt[0] *= MMPERMETER;               /* converts back to mm */

    // add glacier outflow to melt and glmelt
    melt[0] += snow->gl_overflow;
    snow->glmelt += snow->gl_overflow;

    // store final values
    snow->mass_error = mass_balance_error;
    snow->surf_coldcontent = surface_cc;
    snow->pack_coldcontent = surface_pack_cc;
    snow->vapor_flux *= -1.;
    *save_advection = advection;
    *save_deltaCC = delta_cc;
    *save_grnd_flux = grnd_flux;
    *save_latent = latent_heat;
    *save_latent_sub = latent_heat_sub;
    *save_sensible = sensible_heat;
    *save_advected_sensible = advected_sensible_heat;
    *save_refreeze_energy = refreeze_energy;
    *save_qnet = qnet;

    return (0);
}

/*****************************************************************************
   Function name: CalcSnowPackEnergyBalance()

   Purpose      : Dummy function to make a direct call to
                 SnowEnergyBalance() possible.

   Required     :
    double TSurf - SnowPack surface temperature (C)
    other arguments required by SnowPackEnergyBalance()

   Returns      :
    double qnet - Net energy exchange at the SnowPack snow surface (W/m^2)

   Modifies     : none

   Comments     : function is local to this module
*****************************************************************************/
double
CalcSnowPackEnergyBalance(double Tsurf,
                          ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    double  qnet;                /* Net energy exchange at the SnowPack snow
                                    surface (W/m^2) */

    va_start(ap, Tsurf);
    qnet = SnowPackEnergyBalance(Tsurf, ap);
    va_end(ap);

    return qnet;
}

double
ErrorSnowPackEnergyBalance(double Tsurf,
                           ...)
{
    va_list ap;                 /* Used in traversing variable argument list
                                 */
    double  qnet;                /* Net energy exchange at the SnowPack snow
                                    surface (W/m^2) */

    va_start(ap, Tsurf);
    qnet = ErrorPrintSnowPackEnergyBalance(Tsurf, ap);
    va_end(ap);

    return qnet;
}

double
ErrorPrintSnowPackEnergyBalance(double  TSurf,
                                va_list ap)
{
    /* Define Variable Argument List */

    /* General Model Parameters */
    int    rec, iveg, band;
    double Dt;                    /* Model time step (sec) */

    /* Vegetation Parameters */
    double Ra;                    /* Aerodynamic resistance (s/m) */
    double Displacement;          /* Displacement height (m) */
    double Z;                     /* Reference height (m) */
    double z0;                    /* surface roughness height (m) */

    /* Atmospheric Forcing Variables */
    double AirDens;               /* Density of air (kg/m3) */
    double EactAir;               /* Actual vapor pressure of air (Pa) */
    double long_snow_in;            /* Incoming longwave radiation (W/m2) */
    double Lv;                    /* Latent heat of vaporization (J/kg3) */
    double Press;                 /* Air pressure (Pa) */
    double Rain;                  /* Rain fall (m/timestep) */
    double ShortRad;              /* Net incident shortwave radiation
                                     (W/m2) */
    double Vpd;           /* Vapor pressure deficit (Pa) */
    double Wind;                  /* Wind speed (m/s) */

    /* Snowpack Variables */
    double old_tsurf;              /* Surface temperature during previous time
                                      step */
    double SnowCoverFract;        /* Fraction of area covered by snow */
    double SnowDensity;           /* Density of snowpack (kg/m^3) */
    double SurfaceLiquidWater;    /* Liquid water in the surface layer (m) */
    double SweSurfaceLayer;       /* Snow water equivalent in surface layer
                                     (m) */

    /* Energy Balance Components */
    double  Tair;                 /* Canopy surface temperature (C) */
    double  TGrnd;                /* Ground surface temperature (C) */

    double *AdvectedEnergy;       /* Energy advected by precipitation (W/m2) */
    double *AdvectedSensibleHeat; /* Sensible heat advected from snow-free
                                     area into snow covered area (W/m^2) */
    double *DeltaColdContent;     /* Change in cold content of surface
                                     layer (W/m2) */
    double *DeltaPackColdContent; /* Change in sold content of pack
                                     layer (W/m^2) */
    double *GroundFlux;       /* Ground Heat Flux (W/m2) */
    double *LatentHeat;       /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub;        /* Latent heat of sub exchange at
                                     surface (W/m2) */
    double *net_long_snow;          /* Net longwave radiation at snowpack
                                       surface (W/m^2) */
    double *refreeze_energy;       /* Refreeze energy (W/m2) */
    double *SensibleHeat;     /* Sensible heat exchange at surface
                                 (W/m2) */
    double *VaporMassFlux;        /* Mass flux of water vapor to or from the
                                     intercepted snow */
    double *BlowingMassFlux;        /* Mass flux of water vapor to or from the
                                       intercepted snow */
    double *SurfaceMassFlux;        /* Mass flux of water vapor to or from the
                                       intercepted snow */

    char   *errorstring;

    /* Read Variable Argument List */

    /* General Model Parameters */
    rec = (int) va_arg(ap, int);
    iveg = (int) va_arg(ap, int);
    band = (int) va_arg(ap, int);
    Dt = (double) va_arg(ap, double);

    /* Vegetation Parameters */
    Ra = (double) va_arg(ap, double);
    Displacement = (double) va_arg(ap, double);
    Z = (double) va_arg(ap, double);
    z0 = (double) va_arg(ap, double);

    /* Atmospheric Forcing Variables */
    AirDens = (double) va_arg(ap, double);
    EactAir = (double) va_arg(ap, double);
    long_snow_in = (double) va_arg(ap, double);
    Lv = (double) va_arg(ap, double);
    Press = (double) va_arg(ap, double);
    Rain = (double) va_arg(ap, double);
    ShortRad = (double) va_arg(ap, double);
    Vpd = (double) va_arg(ap, double);
    Wind = (double) va_arg(ap, double);

    /* Snowpack Variables */
    old_tsurf = (double) va_arg(ap, double);
    SnowCoverFract = (double) va_arg(ap, double);
    SnowDensity = (double) va_arg(ap, double);
    SurfaceLiquidWater = (double) va_arg(ap, double);
    SweSurfaceLayer = (double) va_arg(ap, double);

    /* Energy Balance Components */
    Tair = (double) va_arg(ap, double);
    TGrnd = (double) va_arg(ap, double);

    AdvectedEnergy = (double *) va_arg(ap, double *);
    AdvectedSensibleHeat = (double *)va_arg(ap, double *);
    DeltaColdContent = (double *) va_arg(ap, double *);
    DeltaPackColdContent = (double *) va_arg(ap, double *);
    GroundFlux = (double *) va_arg(ap, double *);
    LatentHeat = (double *) va_arg(ap, double *);
    LatentHeatSub = (double *) va_arg(ap, double *);
    net_long_snow = (double *) va_arg(ap, double *);
    refreeze_energy = (double *) va_arg(ap, double *);
    SensibleHeat = (double *) va_arg(ap, double *);
    VaporMassFlux = (double *) va_arg(ap, double *);
    BlowingMassFlux = (double *) va_arg(ap, double *);
    SurfaceMassFlux = (double *) va_arg(ap, double *);
    errorstring = (char *) va_arg(ap, char *);

    /* print variables */
    fprintf(stderr, "%s", errorstring);
    fprintf(stderr,
            "ERROR: snow_melt failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

    /* general model terms */
    fprintf(stderr, "rec = %i\n", rec);
    fprintf(stderr, "iveg = %i\n", iveg);
    fprintf(stderr, "band = %i\n", band);
    fprintf(stderr, "Dt = %f\n", Dt);

    /* land surface parameters */
    fprintf(stderr, "Ra = %f\n", Ra);
    fprintf(stderr, "Displacement = %f\n", Displacement);
    fprintf(stderr, "Z = %f\n", Z);
    fprintf(stderr, "z0 = %f\n", z0);

    /* meteorological terms */
    fprintf(stderr, "AirDens = %f\n", AirDens);
    fprintf(stderr, "EactAir = %f\n", EactAir);
    fprintf(stderr, "long_snow_in = %f\n", long_snow_in);
    fprintf(stderr, "Lv = %f\n", Lv);
    fprintf(stderr, "Press = %f\n", Press);
    fprintf(stderr, "Rain = %f\n", Rain);
    fprintf(stderr, "ShortRad = %f\n", ShortRad);
    fprintf(stderr, "Vpd = %f\n", Vpd);
    fprintf(stderr, "Wind = %f\n", Wind);

    /* snow pack terms */
    fprintf(stderr, "old_tsurf = %f\n", old_tsurf);
    fprintf(stderr, "SnowCoverFract = %f\n", SnowCoverFract);
    fprintf(stderr, "SnowDensity = %f\n", SnowDensity);
    fprintf(stderr, "SurfaceLiquidWater = %f\n", SurfaceLiquidWater);
    fprintf(stderr, "SweSurfaceLayer = %f\n", SweSurfaceLayer);
    fprintf(stderr, "Tair = %f\n", Tair);
    fprintf(stderr, "TGrnd = %f\n", TGrnd);
    fprintf(stderr, "AdvectedEnergy = %f\n", AdvectedEnergy[0]);
    fprintf(stderr, "AdvectedSensibleHeat = %f\n", AdvectedSensibleHeat[0]);
    fprintf(stderr, "DeltaColdContent = %f\n", DeltaColdContent[0]);
    fprintf(stderr, "DeltaPackColdContent = %f\n", DeltaPackColdContent[0]);
    fprintf(stderr, "GroundFlux = %f\n", GroundFlux[0]);
    fprintf(stderr, "LatentHeat = %f\n", LatentHeat[0]);
    fprintf(stderr, "LatentHeatSub = %f\n", LatentHeatSub[0]);
    fprintf(stderr, "net_long_snow = %f\n", net_long_snow[0]);
    fprintf(stderr, "refreeze_energy = %f\n", refreeze_energy[0]);
    fprintf(stderr, "SensibleHeat = %f\n", SensibleHeat[0]);
    fprintf(stderr, "VaporMassFlux = %f\n", VaporMassFlux[0]);
    fprintf(stderr, "BlowingMassFlux = %f\n", BlowingMassFlux[0]);
    fprintf(stderr, "SurfaceMassFlux = %f\n", SurfaceMassFlux[0]);

    fprintf(stderr,
            "Finished dumping snow_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThencheck output for instabilities.\n");

    return(ERROR);
}
