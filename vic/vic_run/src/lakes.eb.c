/******************************************************************************
 * @section DESCRIPTION
 *
 * This set of routines handel the energy balance of lakes
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
 * @brief    This subroutine solves the energy budget for open water bodies.
 *****************************************************************************/
int
solve_lake(double           snowfall,
           double           rainfall,
           double           tair,
           double           wind,
           double           vp,
           double           shortin,
           double           longin,
           double           vpd,
           double           pressure,
           double           air_density,
           lake_var_struct *lake,
           soil_con_struct  soil_con,
           double           dt,
           double           wind_h,
           dmy_struct       dmy,
           double           fracprv)
{
    extern parameters_struct param;

    double                   LWnetw, LWneti;
    double                   sw_water, sw_ice;
    double                   T[MAX_LAKE_NODES]; /* temp of the water column, open fraction. */
    double                   Ti[MAX_LAKE_NODES]; /* temp of the water column, ice fraction. */
    double                   water_density[MAX_LAKE_NODES],
                             water_cp[MAX_LAKE_NODES];
    double                   albi, albw;
    double                   tempalbs;
    double                   Tcutoff; /* Lake freezing temperature (K). */
    double                   Qhw, Qhi;
    double                   Qew, Qei;
    double                   qw;
    int                      freezeflag;
    int                      mixdepth;
    double                   new_ice_area; /* Ice area formed by freezing in the open water portion. */
    int                      k, i, ErrorFlag;
    double                   Le;
    double                   sumjoulb;
    double                   temphi;
    double                   new_ice_height;
    double                   windw, windi;
    double                   energy_ice_formation;
    double                   energy_ice_melt_bot;
    double                   Qnet_ice;
    double                   energy_out_bottom;
    double                   energy_out_bottom_ice;
    double                   qf;
    double                   inputs, outputs, internal, phasechange;
    double                   new_ice_water_eq;
    double                   temp_refreeze_energy;
    energy_bal_struct       *lake_energy;
    snow_data_struct        *lake_snow;

    /**********************************************************************
    * 1. Initialize variables.
    **********************************************************************/

    lake->sarea_save = lake->sarea;
    lake->volume_save = lake->volume;
    lake->swe_save = lake->swe;

    lake_energy = &(lake->energy);
    lake_snow = &(lake->snow);

    lake_energy->advection = 0.0;
    lake_energy->deltaCC = 0.0;
    lake_energy->grnd_flux = 0.0;
    lake_energy->snow_flux = 0.0;
    lake->snowmlt = 0.0;
    qw = 0.0;
    new_ice_height = new_ice_area = new_ice_water_eq = 0.0;
    lake->evapw = 0.0;
    energy_ice_formation = 0.0;
    energy_out_bottom = energy_out_bottom_ice = 0.0;
    energy_ice_melt_bot = 0.0;
    lake_snow->vapor_flux = 0.0;
    lake->vapor_flux = 0.0;
    lake_energy->Tsurf = lake->temp[0];
    temp_refreeze_energy = 0.0;
    lake->ice_throughfall = 0.0;

    if (lake->activenod > 0 || lake->areai > 0.0) {
        /* --------------------------------------------------------------------
         * Calculate the water freezing point.
         * ------------------------------------------------------------------- */

        rhoinit(&Tcutoff, pressure);      /* degrees C */

        /* --------------------------------------------------------------------
         * Initialize liquid water and ice temperature profiles to be used
         * later in the module.
         * -------------------------------------------------------------------- */

        for (k = 0; k < lake->activenod; k++) {
            T[k] = lake->temp[k];
            Ti[k] = lake->temp[k];
            water_density[k] = calc_density(T[k]);
            water_cp[k] = specheat(T[k]);
        }

        energycalc(lake->temp, &sumjoulb, lake->activenod, lake->dz,
                   lake->surfdz, lake->surface, water_cp, water_density);

        /**********************************************************************
        * 2. Calculate added precipitation and total snow height.
        **********************************************************************/

        // Convert swq from m/(lake area) to m/(ice area)
        if (lake_snow->swq > 0.0) {
            if (fracprv > 0.0) {
                lake_snow->swq /= fracprv;
            }
            else if (fracprv == 0.0) {
                lake->ice_throughfall += (lake->sarea) * (lake_snow->swq);
                lake_snow->swq = 0.0;
            }
        }

        if (fracprv >= 1.0) { /* areai is relevant */
            // If there is no snow, add the rain over ice directly to the lake.
            if (lake_snow->swq <= 0.0 && rainfall > 0.0) {
                lake->ice_throughfall += (rainfall / MM_PER_M) * lake->areai;
                rainfall = 0.0;
            }
        }
        else if (fracprv > param.LAKE_FRACLIM && fracprv < 1.0) { /* sarea is relevant */
            /* Precip over open water directly increases lake volume. */
            lake->ice_throughfall +=
                ((snowfall / MM_PER_M + rainfall /
                  MM_PER_M) * (1 - fracprv) * lake->sarea);

            // If there is no snow, add the rain over ice directly to the lake.
            if (lake_snow->swq <= 0.0 && rainfall > 0.0) {
                lake->ice_throughfall +=
                    (rainfall / MM_PER_M) * fracprv * lake->sarea;
                rainfall = 0.0; /* Because do not want it added to snow->surf_water */
            }
        }
        else {
            lake->ice_throughfall +=
                ((rainfall + snowfall) / MM_PER_M) * lake->sarea;
            rainfall = 0.0;
            snowfall = 0.0;
        }

        /**********************************************************************
        * 3. Calculate incoming solar radiation over water and ice.
        **********************************************************************/

        /* --------------------------------------------------------------------
         * Calculate the albedo of the lake for ice, snow and liquid water.
         * To be consistent with VIC snow model, tempalbs = NEW_SNOW_ALB if snow
         * is falling, but aging routing does not get reset.
         * -------------------------------------------------------------------- */

        alblake(Tcutoff, tair, &lake->SAlbedo, &tempalbs, &albi, &albw,
                snowfall, lake_snow->coldcontent, dt, &lake_snow->last_snow,
                lake_snow->swq, &lake_snow->MELTING, dmy.day_in_year,
                soil_con.lat);

        /* --------------------------------------------------------------------
         * Calculate the incoming solar radiaton for both the ice fraction
         * and the liquid water fraction of the lake.
         * -------------------------------------------------------------------- */

        if (lake_snow->swq > param.LAKE_SNOWCRIT * param.LAKE_RHOSNOW /
            CONST_RHOFW) {
            sw_ice = shortin * (1. - tempalbs);
            lake_energy->AlbedoLake = (fracprv) * tempalbs +
                                      (1. - fracprv) * albw;
        }
        else if (lake_snow->swq > 0. && lake_snow->swq <= param.LAKE_SNOWCRIT *
                 param.LAKE_RHOSNOW /
                 CONST_RHOFW) {
            sw_ice = shortin * (1. - (albi + tempalbs) / 2.);
            lake_energy->AlbedoLake = (fracprv) *
                                      (albi +
                                       tempalbs) / 2. + (1. - fracprv) * albw;
        }
        else if (fracprv > 0. && lake_snow->swq <= 0.) {
            sw_ice = shortin * (1. - albi);
            lake_energy->AlbedoLake = (fracprv) * albi + (1. - fracprv) * albw;
        }
        else {
            sw_ice = 0.0;
            lake_energy->AlbedoLake = albw;
        }
        lake_energy->AlbedoUnder = lake_energy->AlbedoLake;

        sw_water = shortin * (1. - albw);

        /**********************************************************************
        * 4. Calculate initial energy balance over ice-free water.
        **********************************************************************/

        if ((1. - fracprv) > DBL_EPSILON && lake->activenod > 0) {
            freezeflag = 1; /* Calculation for water, not ice. */
            windw = wind *
                    log((2. + param.LAKE_ZWATER) / param.LAKE_ZWATER) / log(
                wind_h / param.LAKE_ZWATER);

            ErrorFlag = water_energy_balance(lake->activenod, lake->surface,
                                             &lake->evapw,
                                             dt, lake->dz,
                                             lake->surfdz,
                                             soil_con.lat, Tcutoff,
                                             tair, windw,
                                             pressure, vp, air_density, longin,
                                             sw_water, wind_h, &Qhw, &Qew,
                                             &LWnetw, T,
                                             water_density,
                                             &lake_energy->deltaH,
                                             &energy_ice_formation, fracprv,
                                             &new_ice_area, water_cp,
                                             &new_ice_height,
                                             &energy_out_bottom,
                                             &new_ice_water_eq,
                                             lake->volume - lake->ice_water_eq);
            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            /* --------------------------------------------------------------------
             * Do the convective mixing of the lake water.
             * -------------------------------------------------------------------- */

            mixdepth = 0;  /* Set to zero for this time step. */
            tracer_mixer(T, &mixdepth, lake->surface,
                         lake->activenod, lake->dz, lake->surfdz, water_cp);

            lake_energy->AtmosLatent = (1. - fracprv) * Qew;
            lake_energy->AtmosSensible = (1. - fracprv) * Qhw;
            lake_energy->NetLongAtmos = (1. - fracprv) * LWnetw;
            lake_energy->NetShortAtmos = (1. - fracprv) * sw_water;
            lake_energy->refreeze_energy = energy_ice_formation *
                                           (1. - fracprv);
            lake_energy->deltaH *= (1. - fracprv);
            lake_energy->grnd_flux = -1. * (energy_out_bottom * (1. - fracprv));
            lake_energy->Tsurf = (1. - fracprv) * T[0];
        }      /* End of water fraction calculations. */
        else {
            // ice covers 100% of lake, reset open water fluxes
            mixdepth = 0;
            LWnetw = 0;
            Qew = 0;
            Qhw = 0;

            lake_energy->AtmosLatent = 0.0;
            lake_energy->AtmosSensible = 0.0;
            lake_energy->NetLongAtmos = 0.0;
            lake_energy->NetShortAtmos = 0.0;
            lake_energy->refreeze_energy = 0.0;
            lake_energy->deltaH = 0.0;
            lake_energy->grnd_flux = 0.0;
            lake_energy->Tsurf = 0.0;
        }

        /**********************************************************************
        *  6. Calculate initial energy balance over ice.
        **********************************************************************/

        windi = (wind * log((2. + soil_con.snow_rough) / soil_con.snow_rough) /
                 log(wind_h / soil_con.snow_rough));
        if (windi < 1.0) {
            windi = 1.0;
        }

        if (fracprv >= param.LAKE_FRACLIM) {
            freezeflag = 0;   /* Calculation for ice. */
            Le = calc_latent_heat_of_sublimation(tair); /* ice*/

            lake->aero_resist =
                (log((2. + soil_con.snow_rough) / soil_con.snow_rough) *
                 log(wind_h /
                     soil_con.snow_rough) /
                 (CONST_KARMAN * CONST_KARMAN)) / windi;

            /* Calculate snow/ice temperature and change in ice thickness from
               surface melting. */
            ErrorFlag = ice_melt(wind_h + soil_con.snow_rough,
                                 lake->aero_resist, &(lake->aero_resist),
                                 Le, lake_snow, lake, dt, 0.0,
                                 soil_con.snow_rough, 1.0,
                                 rainfall, snowfall, windi, Tcutoff, tair,
                                 sw_ice,
                                 longin, air_density, pressure, vpd, vp,
                                 &lake->snowmlt,
                                 &lake_energy->advection, &lake_energy->deltaCC,
                                 &lake_energy->snow_flux, &Qei, &Qhi, &Qnet_ice,
                                 &temp_refreeze_energy, &LWneti);
            if (ErrorFlag == ERROR) {
                return (ERROR);
            }

            lake_energy->refreeze_energy += temp_refreeze_energy * fracprv;
            lake->tempi = lake_snow->surf_temp;

            /**********************************************************************
            *  7. Adjust temperatures of water column in ice fraction.
            **********************************************************************/

            /* --------------------------------------------------------------------
             * Calculate inputs to temp_area..
             * -------------------------------------------------------------------- */

            if (lake->activenod > 0) {
                ErrorFlag = water_under_ice(freezeflag, sw_ice, wind, Ti,
                                            water_density,
                                            soil_con.lat,
                                            lake->activenod, lake->dz,
                                            lake->surfdz,
                                            Tcutoff, &qw, lake->surface,
                                            &temphi, water_cp,
                                            mixdepth, lake->hice,
                                            lake_snow->swq * CONST_RHOFW / param.LAKE_RHOSNOW,
                                            dt, &energy_out_bottom_ice);
                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }
            else {
                temphi = -sumjoulb;
            }

            /**********************************************************************
            *   8.  Calculate change in ice thickness and fraction
            *    within fraction that already has ice.
            **********************************************************************/
            /* Check to see if ice has already melted (from the top) in this time step. */
            if (lake->ice_water_eq > 0.0) {
                ErrorFlag = lakeice(sw_ice,
                                    fracprv, dt, lake_energy->snow_flux, qw,
                                    &energy_ice_melt_bot,
                                    lake_energy->deltaCC, &qf,
                                    &lake->ice_water_eq,
                                    lake->volume - new_ice_water_eq,
                                    lake->surface[0]);
                if (ErrorFlag == ERROR) {
                    return (ERROR);
                }
            }

            lake_energy->AtmosLatent += fracprv * Qei;
            lake_energy->advection *= fracprv;
            lake_energy->AtmosSensible += fracprv * Qhi;
            lake_energy->NetLongAtmos += fracprv * LWneti;
            lake_energy->NetShortAtmos += fracprv * sw_ice;
            lake_energy->deltaH += fracprv * temphi;
            lake_energy->grnd_flux += -1. * (energy_out_bottom_ice * fracprv);
            lake_energy->refreeze_energy += energy_ice_melt_bot * fracprv;
            lake_energy->Tsurf += fracprv * lake_snow->surf_temp;
        }
        else {
            /* No Lake Ice Fraction */
            LWneti = 0;
            Qei = 0.;
            Qhi = 0;
            qf = 0.0;
            temphi = 0.0;
            lake_energy->refreeze_energy = 0.0;
            if (fracprv > 0.0) {
                energy_ice_melt_bot =
                    (lake->hice * CONST_RHOICE +
                     (snowfall /
                      MM_PER_M) *
                     CONST_RHOFW) * CONST_LATICE / dt;
                lake->areai = 0.0;
                lake->hice = 0.0;
                lake->ice_water_eq = 0.0;
            }
            else {
                energy_ice_melt_bot = 0.0;
                lake->areai = 0.0;
                lake->hice = 0.0;
                lake->ice_water_eq = 0.0;
            }
            lake->aero_resist =
                (log((2. + param.LAKE_ZWATER) / param.LAKE_ZWATER) *
                 log(wind_h / param.LAKE_ZWATER) /
                 (CONST_KARMAN * CONST_KARMAN)) / windi;
        }
        lake->soil.aero_resist[0] = lake->aero_resist;

        /**********************************************************************
        * 9. Average water temperature.
        **********************************************************************/

        if (lake->activenod > 0) {
            // Average ice-covered and non-ice water columns.
            colavg(lake->temp, T, Ti, fracprv, lake->density, lake->activenod,
                   lake->dz, lake->surfdz);

            // Calculate depth average temperature of the lake
            lake->tempavg = 0.0;
            for (i = 0; i < lake->activenod; i++) {
                lake->tempavg += lake->temp[i] / lake->activenod;
            }
        }
        else {
            lake->tempavg = -99;
        }

        /**********************************************************************
        * 10. Calculate the final water heat content and energy balance.
        **********************************************************************/

        /* Incoming energy. */
        inputs = (sw_ice + LWneti + lake_energy->advection + Qhi + Qei);
        outputs = energy_out_bottom_ice;
        internal = temphi;
        phasechange = -1 * (lake_energy->refreeze_energy) - 1. *
                      energy_ice_melt_bot;

        lake_energy->error = inputs - outputs - internal - phasechange;

        lake_energy->snow_flux = 0.0;

        // Sign convention
        lake_energy->deltaH *= -1;

        lake_energy->error = (lake_energy->NetShortAtmos +
                              lake_energy->NetLongAtmos +
                              lake_energy->AtmosSensible +
                              lake_energy->AtmosLatent +
                              lake_energy->deltaH +
                              lake_energy->grnd_flux +
                              lake_energy->refreeze_energy +
                              lake_energy->advection);

        temphi = 0.0;

        /**********************************************************************
        * 11. Final accounting for passing variables back to VIC.
        **********************************************************************/

        // Adjust lake_snow variables to represent storage and flux over entire lake
        lake_snow->swq *= fracprv;
        lake_snow->surf_water *= fracprv;
        lake_snow->pack_water *= fracprv;
        lake_snow->depth = lake_snow->swq * CONST_RHOFW / param.LAKE_RHOSNOW;
        lake_snow->coldcontent *= fracprv;
        lake_snow->vapor_flux *= fracprv;
        lake_snow->blowing_flux *= fracprv;
        lake_snow->surface_flux *= fracprv;
        lake_snow->melt = lake->snowmlt * fracprv;    // in mm

        /* Lake data structure terms */
        lake->evapw *=
            ((1. - fracprv) * dt) / MM_PER_M * lake->sarea;            // in m3
        lake->vapor_flux = lake->snow.vapor_flux * lake->sarea; // in m3
        lake->swe = lake_snow->swq * lake->sarea; // in m3
        lake->snowmlt *= fracprv / MM_PER_M * lake->sarea; // in m3
        lake->pack_water = lake_snow->pack_water * lake->sarea; // in m3
        lake->surf_water = lake_snow->surf_water * lake->sarea; // in m3
        lake->sdepth = lake_snow->depth * lake->sarea;
        lake->pack_temp = lake_snow->pack_temp;
        lake->surf_temp = lake_snow->surf_temp;

        /* Update ice area to include new ice growth in water fraction. */
        lake->new_ice_area = lake->areai;
        if (new_ice_area > 0.0) {
            lake->new_ice_area += new_ice_area;
            lake->ice_water_eq += new_ice_water_eq;
        }
        if (lake->ice_water_eq > 0.0 && lake->new_ice_area > 0.0) {
            lake->hice =
                (lake->ice_water_eq /
                 lake->new_ice_area) * CONST_RHOFW / CONST_RHOICE;
        }
        else {
            lake->hice = 0.0;
        }

        /* Change area of ice-covered fraction if ice has thinned. */
        if (lake->hice <= 0.0) {
            lake->new_ice_area = 0.0;
            lake->hice = 0.0;
        }
        else if (lake->hice < param.LAKE_FRACMIN) {
            lake->new_ice_area =
                (lake->new_ice_area * lake->hice) / param.LAKE_FRACMIN;
            lake->hice = param.LAKE_FRACMIN;
        }

        if (lake->snow.swq > 0) {
            lake->snow.coverage = lake->new_ice_area / lake->sarea;
        }
        else {
            lake->snow.coverage = 0;
        }
    } /* End of if activenods > 0 */


    return (0);
}   /* End of solve_lake function. */

/******************************************************************************
 * @brief    Calculate the partitioning of the energy balance into latent and
 *           sensible heat fluxes.
 *****************************************************************************/
void
latsens(double  Tsurf,
        double  Tcutk,
        double  hice,
        double  tair,
        double  wind,
        double  pressure,
        double  vp,
        double  air_density,
        double *evap,
        double *qsen,
        double  wind_h)
{
    extern parameters_struct param;

    double                   dragcoeff;
    double                   eog, delT;
    double                   qair, qlake; /* Specific humidity of the atmosphere and lake, respectively. */
    double                   delq; /* Difference in absolute humidity between the lake */
    /* surface and higher up (-). */

/**********************************************************************
* Calculate the drag coefficient.
**********************************************************************/

    if (hice > 0.) {
        dragcoeff = lkdrag(Tsurf, tair + CONST_TKFRZ, wind, param.LAKE_ZSNOW,
                           wind_h);
    }
    else {
        dragcoeff = lkdrag(Tsurf, tair + CONST_TKFRZ, wind, param.LAKE_ZWATER,
                           wind_h);
    }

/**********************************************************************
* Determine the coefficients to be used in the calculation of
* the vapor pressure at lake level depending on whether the lake is
* covered with ice or not.
**********************************************************************/

    if ((hice <= 0.) && (Tsurf > Tcutk)) {
/* --------------------------------------------------------------------
 * Lake is not covered with ice (Handbook of Hydrology eq. 4.2.2).
 * eog in kPa.
 * -------------------------------------------------------------------- */
        eog = .611 *
              exp(17.269 *
                  (Tsurf - CONST_TKFRZ) / (Tsurf + 237.3 - CONST_TKFRZ));
    }
    else {
/* --------------------------------------------------------------------
 * Lake is covered with ice.
 * -------------------------------------------------------------------- */

        eog = .611 * exp(21.874 * (Tsurf - CONST_TKFRZ) / (Tsurf - 7.66));
    }


/**********************************************************************
* Calculate the specific humidity at lake level and the difference
* between this humidity and the air humidity at measurement height.
**********************************************************************/

    qlake = CONST_EPS * (eog / (pressure - 0.378 * eog));
    qair = CONST_EPS * (vp / (pressure - 0.378 * vp));
    delq = qair - qlake;

    /**********************************************************************
    * Calculate the evaporation rate.  Eq. 4 in Hostetler (1991),
    * after Brutsaert (1982).  evap in mm/s
    **********************************************************************/

    *evap = -1 * dragcoeff * wind * air_density * delq;

/**********************************************************************
* Calculate the difference in lake surface temperature and the
* air temperature at measurement height and consequently the
* sensible heat flux.  Hostetler (1991) eq 5, in W/m^2.
**********************************************************************/

    delT = tair + CONST_TKFRZ - Tsurf;
    *qsen = dragcoeff * wind * air_density * CONST_CPMAIR;
    *qsen = *qsen * delT;
}

/******************************************************************************
 * @brief    Calculate the albedo of snow, ice and water of the lake.
 *****************************************************************************/
void
alblake(double         Tcutoff,
        double         Tair,
        double        *snowalbedo,
        double        *albs,
        double        *albi,
        double        *albw,
        double         newsnow,
        double         coldcontent,
        double         dt,
        unsigned      *last_snow,
        double         swq,
        bool          *MELTING,
        unsigned short day_in_year,
        double         latitude)
{
    extern parameters_struct param;

    double                   albgl, albgs;

    if ((Tair - Tcutoff) > 0.0) {
        if ((Tair - Tcutoff) < 20.) {
            albgl = 0.4 - 0.011 * (Tair - Tcutoff);
            albgs = 0.6 - 0.0245 * (Tair - Tcutoff);
        }
        else {
            albgl = 0.4 - 0.011 * 20.;
            albgs = 0.6 - 0.0245 * 20.;
        }
    }
    else {
        albgl = 0.4;
        albgs = 0.6;
    }

    *albi = 0.5 * albgs + 0.5 * albgl;

    // update number of days since last significant snowfall
    if (newsnow > param.SNOW_TRACESNOW) {
        *last_snow = 1;
    }
    else if (swq == 0) {
        *last_snow = 0;
    }
    else {
        *last_snow += 1;
    }

    /** Record if snowpack is melting this time step **/
    if (swq > 0.0) {
        if (coldcontent >= 0 && (
                (latitude >= 0 && (day_in_year > 60 && // ~ March 1
                                   day_in_year < 273)) || // ~ October 1
                (latitude < 0 && (day_in_year < 60 || // ~ March 1
                                  day_in_year > 273)) // ~ October 1
                )) {
            *MELTING = true;
        }
        else {
            *MELTING = false;
        }
    }
    else {
        *MELTING = false;
    }

    if (*MELTING && newsnow > param.SNOW_TRACESNOW) {
        *MELTING = false;
    }

    // compute snow surface albedo
    if (swq > 0.0) {
        *snowalbedo = snow_albedo(newsnow, swq, *snowalbedo, coldcontent,
                                  dt, *last_snow, *MELTING);
    }
    else if (swq == 0.0 && newsnow > 0.0) {
        *snowalbedo = param.SNOW_NEW_SNOW_ALB;
    }
    else {
        *snowalbedo = 0.0;
    }

    if (newsnow > 0.0) {
        *albs = param.SNOW_NEW_SNOW_ALB;
    }
    else {
        *albs = *snowalbedo;
    }

    *albw = 0.15;
}

/******************************************************************************
 * @brief    Calculate the temperature of the lake water at the different node
 *           levels and calculate the water density.
 *****************************************************************************/
void
colavg(double *finaltemp,
       double *T,
       double *Ti,
       double  lakeprv,
       double *density,
       int     numnod,
       double  dz,
       double  surfdz)
{
    int    j;
    double water_densityw, water_densityi;
    double temp;
    double z;

    for (j = 0; j < numnod; j++) {
/* --------------------------------------------------------------------
 * Calculate the densities of the ice and water fractions.
 * -------------------------------------------------------------------- */

        water_densityw = calc_density(T[j]);
        water_densityi = calc_density(Ti[j]);
        water_densityw = water_densityw + CONST_RHOFW;
        water_densityi = water_densityi + CONST_RHOFW;

/* --------------------------------------------------------------------
 * Calculate the depth differences.
 * -------------------------------------------------------------------- */

        z = dz;
        if (j == 0) {
            z = surfdz;
        }

/* --------------------------------------------------------------------
 * Calculate the lake (water) temperature as a weight of ice and water
 * temperatures.
 * -------------------------------------------------------------------- */

        temp = ((1. - lakeprv) * T[j] * z * water_densityw +
                lakeprv * Ti[j] * z * water_densityi) /
               ((1. -
                 lakeprv) * z * water_densityw + lakeprv * z * water_densityi);

        finaltemp[j] = temp;

        /* --------------------------------------------------------------------
         * Recalculate the water density at the different nodes.
         * -------------------------------------------------------------------- */

        density[j] = calc_density(finaltemp[j]);
    }
}

/******************************************************************************
 * @brief    Calculate the water density.
 *****************************************************************************/
double
calc_density(double ts)
{
    double rhow, t, rhostps;

    t = ts;

/* --------------------------------------------------------------------
 * Calculate the water density as a function of temperature.
 * -------------------------------------------------------------------- */

    rhow = 999.842594 + .06793952 * t - .009095290 * t * t +
           .0001001685 * t * t * t - .000001120083 * t * t * t * t +
           .000000006536332 * t * t * t * t * t;


/* --------------------------------------------------------------------
 * Return the difference between 1000 kg/m3 and the calculated
 * density, not the actual value.
 * -------------------------------------------------------------------- */

    rhostps = rhow - CONST_RHOFW;
    return (rhostps);
}

/******************************************************************************
 * @brief    Calculate the eddy diffusivity.
 *****************************************************************************/
void
eddy(int     freezeflag,
     double  wind,
     double *water_density,
     double *de,
     double  lat,
     int     numnod,
     double  dz,
     double  surfdz)
{
    extern parameters_struct param;

    double                   ks, N2, ws, radmax, Po;
    double                   dpdz, rad;
    double                   zhalf[MAX_LAKE_NODES];
    int                      k;
    double                   z;
    double                   Ri; /* Richardson's number. */

/**********************************************************************
* Calculate the density at all nodes and the distance between all nodes.
**********************************************************************/

    for (k = 0; k < numnod; k++) {
        zhalf[k] = dz;
    }

/**********************************************************************
* Calculate the distance between the first node and the surface node.
**********************************************************************/

    zhalf[0] = (surfdz + dz) * 0.5;

/**********************************************************************
* If there is ice only molecular diffusivity is taken in account,
* no eddy diffusivity.
**********************************************************************/

    if (freezeflag != 1) {
        for (k = 0; k < numnod; k++) {
            de[k] = param.LAKE_DM;
        }
    }

    /**********************************************************************
    * Avoid too low wind speeds for computational stability.
    **********************************************************************/

    else {
        if (wind < 1.0) {
            wind = 1.0;
        }

        /**********************************************************************
        * Determine the latitudinaly dependent parameter of the Ekman
        * profile, the surface value of the friction velocity and the neutral
        * value of the Prandtl number (Hostetler and Bartlein eq. 6 and 7).
        **********************************************************************/

        ks = 6.6 * pow(sin((double) fabs(lat) * CONST_PI / 180.), 0.5) * pow(
            wind,
            -1.84);
        ws = 0.0012 * wind;
        Po = 1.0;

        /**********************************************************************
        * Determine the maximum value of the second term in the calculation
        * of the Richardson number for computational stability.
        **********************************************************************/

        radmax = 4.e4;

        for (k = 0; k < numnod - 1; k++) {
            /**********************************************************************
            * Calculate the eddy diffusivity for each node.
            **********************************************************************/

            /* --------------------------------------------------------------------
             * Calculate the derivative of density with depth, the Brunt-Vaisala
             * frequency and the depth of the node (Hostetler and Bartlein eq. 9).
             * -------------------------------------------------------------------- */

            dpdz = (water_density[k + 1] - water_density[k]) / zhalf[k];
            N2 = (dpdz / (1.e3 + water_density[k])) * 9.8;
            z = surfdz + ((double) k) * dz;

            /* --------------------------------------------------------------------
             * Calculate the second term in the calculation of the Richardson
             * number, make sure this number does not get too large for
             * computational stability.  Also make sure this term does not go
             * below 1 to avoid negative Richardson numbers (Hostetler and Bartlein eq. 8).
             * -------------------------------------------------------------------- */

            if ((z * exp(ks * z) / ws) > 1.e8) {
                rad = radmax;
            }
            else {
                rad = 1. + 40. * N2 *
                      (CONST_KARMAN *
                       z) * (CONST_KARMAN * z) / (ws * ws * exp(-2. * ks * z));

                if (rad > radmax) {
                    rad = radmax;
                }

                if (rad < 1.0) {
                    rad = 1.0;
                }
            }

            /* --------------------------------------------------------------------
             * Calculate the Richardson number and the eddy diffusivity.
             * (Hostetler and Bartlein eq. 8 and 5.
             * -------------------------------------------------------------------- */

            Ri = (-1.0 + sqrt(rad)) / 20.0;
            de[k] = param.LAKE_DM +
                    (CONST_KARMAN * ws * z / Po) * exp(-ks * z) /
                    (1.0 + 37.0 * Ri * Ri);
        }

        /* --------------------------------------------------------------------
         * The eddy diffusivity of the last node is assumed to equal the
         * eddy diffusivity of the second last node.
         * -------------------------------------------------------------------- */

        de[numnod - 1] = de[numnod - 2];
    }
}

/******************************************************************************
 * @brief    Calculate the form of new ice in the lake as long as the
 *           fractional coverage of the ice is not 1 yet.
 *****************************************************************************/
void
iceform(double *qfusion,
        double *T,
        double  Tcutoff,
        double  fracprv,
        double *areaadd,
        int     numnod,
        double  dt,
        double  dz,
        double  surfdz,
        double *cp,
        double *surface,
        double *new_ice_height,
        double *new_ice_water_eq,
        double  lvolume)
{
    extern parameters_struct param;

    double                   sum, extra;
    int                      j;
    double                   di;

/**********************************************************************
* Calculate the newly added ice to the ice layer.
**********************************************************************/
    *qfusion = 0.0;
    sum = 0.0;
    *new_ice_water_eq = 0.0;

    for (j = 0; j < numnod; j++) {
        if (T[j] < Tcutoff) {
/* --------------------------------------------------------------------
 * Calculate the ice growth (extra in J).
 * -------------------------------------------------------------------- */

            if (j == 0) {
                extra = ((Tcutoff - T[j]) * surfdz * CONST_RHOFW *
                         cp[j] *
                         (1.0 - fracprv) * (surface[j] + surface[j + 1]) / 2.);
            }
            else if (j == numnod - 1) {
                extra = ((Tcutoff - T[j]) * dz * CONST_RHOFW *
                         cp[j] * (1.0 - fracprv) * surface[j]);
            }
            else {
                extra = ((Tcutoff - T[j]) * dz * CONST_RHOFW *
                         cp[j] *
                         (1.0 - fracprv) * (surface[j] + surface[j + 1]) / 2.);
            }

/* --------------------------------------------------------------------
 * If ice, the temperature of the water is the freezing temperature.
 * -------------------------------------------------------------------- */

            T[j] = Tcutoff;

            sum += extra;
        }
    }

/**********************************************************************
* Calculate the heat flux absorbed into the ice (in W/m2) and the
* thickness of the new ice.
**********************************************************************/

    *new_ice_water_eq = sum / (CONST_RHOFW * CONST_LATICE);

    if (lvolume > *new_ice_water_eq) {
        *qfusion = (sum / (dt * surface[0] * (1.0 - fracprv))); /*W/m2*/

        di = sum / (CONST_LATICE * CONST_RHOICE);    /* m^3 of ice formed */
    }
    else if (lvolume > 0.0) {
        *new_ice_water_eq = lvolume;

        di = *new_ice_water_eq * CONST_RHOFW / CONST_RHOICE;

        // NEED TO CHANGE ICE TEMPERATURE TO ACCOUNT FOR EXTRA qfusion
        *qfusion =
            (*new_ice_water_eq * CONST_RHOFW /
             CONST_RHOICE) / (dt * surface[0] * (1.0 - fracprv));                    /*W/m2*/
    }
    else {
        *new_ice_water_eq = 0.0;

        di = 0.0;

        *qfusion = 0.0; /*W/m2*/
    }

    /**********************************************************************
    * Calculate the added fractional coverage of ice, make sure the
    * total fractional coverage does not exceed 1.
    **********************************************************************/

    *new_ice_height = param.LAKE_FRACMIN;
    *areaadd = di / (param.LAKE_FRACMIN);

    if (*areaadd > (1.0 - fracprv) * surface[0]) {
        *new_ice_height = di / ((1. - fracprv) * surface[0]);
        *areaadd = (1. - fracprv) * surface[0];
    }
}

/******************************************************************************
 * @brief    Calculate the radiation balance over ice.
 *****************************************************************************/
void
icerad(double  sw,
       double  hi,
       double  hs,
       double *avgcond,
       double *SWnet,
       double *SW_under_ice)
{
    extern parameters_struct param;

    double                   a, b, c, d;

/**********************************************************************
* Calculate the thermal conductivity of the combined snow and ice
* layer.
**********************************************************************/

    *avgcond =
        (hs * param.LAKE_CONDI + hi *
         param.LAKE_CONDS) / (param.LAKE_CONDI * param.LAKE_CONDS);

/**********************************************************************
* Calculate the incoming radiation at different levels in the
* ice-snow layer.
**********************************************************************/

/* --------------------------------------------------------------------
 * Calculation of constants. Patterson and Hamblin eq, 7
 * -------------------------------------------------------------------- */

    a = -1. *
        (1. -
         exp(-param.LAKE_LAMSSW * hs)) / (param.LAKE_CONDS * param.LAKE_LAMSSW);
    b = -1. *
        exp(-param.LAKE_LAMSSW *
            hs) *
        (1 -
         exp(-param.LAKE_LAMISW * hi)) / (param.LAKE_CONDI * param.LAKE_LAMISW);
    c = -1. *
        (1. -
         exp(-param.LAKE_LAMSLW * hs)) / (param.LAKE_CONDS * param.LAKE_LAMSLW);
    d = -1. *
        exp(-param.LAKE_LAMSLW *
            hs) *
        (1 -
         exp(-param.LAKE_LAMILW * hi)) / (param.LAKE_CONDI * param.LAKE_LAMILW);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow pack. RHS of Patterson and Hamblin, eq. 7
 * -------------------------------------------------------------------- */

    *SWnet = sw * param.LAKE_A1 * (a + b) + sw * param.LAKE_A2 * (c + d);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow/ice layer.
 * Patterson and Hamblin eq. 8 (qf-qo)
 * -------------------------------------------------------------------- */

    *SW_under_ice =
        (param.LAKE_A1 * sw *
         (1 - exp(-(param.LAKE_LAMSSW * hs + param.LAKE_LAMISW * hi))) +
         param.LAKE_A2 * sw *
         (1 - exp(-(param.LAKE_LAMSLW * hs + param.LAKE_LAMILW * hi))));
}

/******************************************************************************
 * @brief    Calculate the growth and decrease in the lake ice cover.
 *****************************************************************************/
int
lakeice(double  sw_ice,
        double  fracice,
        double  dt,
        double  snowflux,
        double  qw,
        double *energy_ice_melt_bot,
        double  SWabsorbed,
        double *qf,
        double *ice_water_eq,
        double  volume,
        double  sarea)
{
    double dibot;          /* change in ice surface at the bottom of the pack. */
    double new_water_eq;

    /**********************************************************************
    * Calculate fluxes at the base of the ice.
    **********************************************************************/

    /* --------------------------------------------------------------------
     * Flux of heat in the ice at the ice/water interface (P & H, eq. 8)
     * -------------------------------------------------------------------- */

    *qf = snowflux + sw_ice - SWabsorbed;

    /* --------------------------------------------------------------------
     * Amount of heat used to melt the ice (positive means freezing).
     * -------------------------------------------------------------------- */

    *energy_ice_melt_bot = *qf - qw; // from Hostetler 1991

    /* --------------------------------------------------------------------
     * Calculate the growth of the ice pack at the bottom (in meters).
     * -------------------------------------------------------------------- */

    dibot =
        (*energy_ice_melt_bot /
         (CONST_RHOICE * CONST_LATICE)) * dt;

    /* --------------------------------------------------------------------
     * Calculate the water equivalent of the ice volume (in cubic meters).
     * Freezing occurs over surface area of water, not ice, if ice area exceeds
     * water area.
     * -------------------------------------------------------------------- */

    new_water_eq = dibot * sarea * (fracice) * CONST_RHOICE / CONST_RHOFW;

    /**********************************************************************
    * Calculate the new height of the ice pack and the fractional ice
    * cover.
    **********************************************************************/
    if (dibot > 0.0) { /*Freezing; check if enough unfrozen water is available. */
        if (volume - *ice_water_eq >= new_water_eq) {
            *ice_water_eq += new_water_eq;
        }
        else { /* Freezing is restricted by available water. */
            if ((volume - *ice_water_eq) > 0.0) {
                dibot =
                    (volume -
                     *ice_water_eq) /
                    (sarea * (fracice) * CONST_RHOICE / CONST_RHOFW);
                *ice_water_eq = volume;
            }
            else {
                dibot = 0.0;
            }
        }
    }
    else { /* Melt */
        *ice_water_eq += new_water_eq;

        // check that ice not completely melted in the current time step
        if (*ice_water_eq <= 0.0) {
            *ice_water_eq = 0.0;
        }
    }
    // *energy_ice_melt_bottom is not currently adjusted if there is not enough water to freeze or ice to melt.
    // This energy should go to cooling the ice and warming the water, respectively.

    return (0);
}

/******************************************************************************
 * @brief    Calculate the lake drag coefficient.
 *****************************************************************************/
double
lkdrag(double Tsurf,
       double Tair,
       double wind,
       double roughness,
       double Z1)
{
    double cdrn, ribn, ribd, rib;
    double cdr, cdrmin;

/**********************************************************************
* Calculate the Richardson number.
**********************************************************************/

    cdrn =
        (CONST_KARMAN /
         log(Z1 / roughness)) * (CONST_KARMAN / log(Z1 / roughness));                   /* dimensionless */

    ribn = Z1 * CONST_G * (1. - Tsurf / Tair);  /* m2/s2 */

    if ((Tsurf / Tair) <= 1.0) {
        ribd = wind * wind + 0.1 * 0.1;
    }
    else {
        ribd = wind * wind + 1.0 * 1.0;
    }

    rib = ribn / ribd;              /* dimensionless */

/**********************************************************************
* Calculate the drag coefficient using the Richardson number.
**********************************************************************/

    if (rib < 0.) {
        cdr = cdrn * (1.0 + 24.5 * sqrt(-cdrn * rib));
    }
    else {
        cdr = cdrn / (1.0 + 11.5 * rib);
    }

    if ((.25 * cdrn) > 6.e-4) {
        cdrmin = .25 * cdrn;
    }
    else {
        cdrmin = 6.e-4;
    }

    if (cdr < cdrmin) {
        cdr = cdrmin;
    }

    return(cdr);
}

/******************************************************************************
 * @brief    Calculate the temperature at which water freezes depending on
 *           salinity and air pressure.
 *****************************************************************************/
void
rhoinit(double *tfsp,
        double  pressure)
{
    double salinity;

/**********************************************************************
* Salinity is assumed to be zero.
**********************************************************************/

    salinity = 0.;

/**********************************************************************
* Calculate the lake freezing temperature (C). Pressure in bars.
**********************************************************************/

    *tfsp = (-0.0575 * salinity + 1.710523e-3 * pow(salinity, 1.5) -
             2.154996e-4 * salinity * salinity - 7.53e-3 * (pressure) /
             BAR_PER_KPA);
}

/******************************************************************************
 * @brief    Calculate the specific heat of the water depending on water
 *           temperature.  Salinity is assumed to be zero.
 *****************************************************************************/
double
specheat(double t)
{
    double cpt;                 /* Specific heat (J/Kg K) */

    cpt = 4217.4 - 3.720283 * t + 0.1412855 * t * t - 2.654387e-3 * t * t * t +
          2.093236e-5 * t * t * t * t;

    return cpt;
}

/******************************************************************************
 * @brief    Calculate the water temperature for different levels in the lake.
 *****************************************************************************/
void
temp_area(double  sw_visible,
          double  sw_nir,
          double  surface_force,
          double *T,
          double *Tnew,
          double *water_density,
          double *de,
          double  dt,
          double *surface,
          int     numnod,
          double  dz,
          double  surfdz,
          double *temph,
          double *cp,
          double *energy_out_bottom)
{
    extern parameters_struct param;

    double                   z[MAX_LAKE_NODES], zhalf[MAX_LAKE_NODES];
    double                   a[MAX_LAKE_NODES], b[MAX_LAKE_NODES],
                             c[MAX_LAKE_NODES];
    double                   d[MAX_LAKE_NODES];

    int                      k;
    double                   surface_1, surface_2, surface_avg, T1;
    double                   cnextra;
    double                   top, bot; /* The depth of the top and the bottom of the current water layer. */
    double                   energyinput;
    double                   joulenew;
    double                   energymixed;
    double                   term1, term2;

/**********************************************************************
* Calculate the distance between the centers of the surface and first
* lake layers.
**********************************************************************/

/**********************************************************************
* Initialize the water density at all the nodes and the depth of all
* and distance between all nodes.
**********************************************************************/

    for (k = 0; k < numnod; k++) {
        if (k == 0) {
            z[k] = surfdz;
        }
        else {
            z[k] = dz;
        }
        zhalf[k] = dz;
    }
    if (numnod > 1) {
        zhalf[0] = 0.5 * (z[0] + z[1]);
    }
    else {
        zhalf[0] = 0.5 * z[0];
    }
    energyinput = 0.0;
    energymixed = 0.0;

/**********************************************************************
* Calculate the right hand side vector in the tridiagonal matrix system
* of equations.
**********************************************************************/

    surface_1 = surface[0];
    surface_2 = surface[1];
    surface_avg = (surface_1 + surface_2) / 2.;

    T1 =
        (sw_visible *
         (1 * surface_1 - surface_2 * exp(-param.LAKE_LAMWSW * surfdz)) +
         sw_nir *
         (1 * surface_1 - surface_2 *
          exp(-param.LAKE_LAMWLW * surfdz))) / surface_avg +
        (surface_force * surface_1) / surface_avg;  /* W/m2 */

    energyinput += T1 * surface_avg;
    cnextra = 0.5 *
              (surface_2 /
               surface_avg) * (de[0] / zhalf[0]) * ((T[1] - T[0]) / z[0]);

    energymixed += cnextra;

    *temph = 0.0;

    if (numnod == 1) {
        Tnew[0] = T[0] +
                  (T1 * dt) / ((1.e3 + water_density[0]) * cp[0] * z[0]);
    }
    else {
        /* --------------------------------------------------------------------
         * First calculate d for the surface layer of the lake.
         * -------------------------------------------------------------------- */

        d[0] = T[0] +
               (T1 * dt) /
               ((1.e3 +
                 water_density[0]) * cp[0] *
                z[0]) + cnextra * dt;

        *energy_out_bottom =
            (surface_1 -
             surface_2) * (sw_visible * exp(-param.LAKE_LAMWSW * surfdz) +
                           sw_nir *
                           exp(-param.LAKE_LAMWLW * surfdz));


        /* --------------------------------------------------------------------
         * Calculate d for the remainder of the column.
         * --------------------------------------------------------------------*/

        /* ....................................................................
         * All nodes but the deepest node.
         * ....................................................................*/

        for (k = 1; k < numnod - 1; k++) {
            top = (surfdz + (k - 1) * dz);
            bot = (surfdz + (k) * dz);

            surface_1 = surface[k];
            surface_2 = surface[k + 1];
            surface_avg = (surface[k] + surface[k + 1]) / 2.;


            T1 =
                (sw_visible *
                 (surface_1 *
                  exp(-param.LAKE_LAMWSW *
                      top) - surface_2 * exp(-param.LAKE_LAMWSW * bot)) +
                 sw_nir *
                 (surface_1 *
                  exp(-param.LAKE_LAMWLW *
                      top) - surface_2 *
                  exp(-param.LAKE_LAMWLW * bot))) / surface_avg;

            energyinput += T1 * surface_avg;
            term1 = 0.5 *
                    (1. /
                     surface_avg) *
                    ((de[k] /
                      zhalf[k]) * ((T[k + 1] - T[k]) / z[k])) * surface_2;
            term2 = 0.5 *
                    (-1. /
                     surface_avg) *
                    ((de[k -
                         1] /
                      zhalf[k - 1]) * ((T[k] - T[k - 1]) / z[k])) * surface_1;

            cnextra = term1 + term2;

            energymixed += (term1 + term2);
            d[k] = T[k] +
                   (T1 * dt) / ((1.e3 + water_density[k]) * cp[k] * z[k]) +
                   cnextra * dt;

            *energy_out_bottom +=
                (surface_1 -
                 surface_2) * (sw_visible * exp(-param.LAKE_LAMWSW * bot) +
                               sw_nir *
                               exp(-param.LAKE_LAMWLW * bot));
        }

        /* ....................................................................
         * Calculation for the deepest node.
         * ....................................................................*/
        k = numnod - 1;
        surface_1 = surface[k];
        surface_2 = surface[k];
        surface_avg = surface[k];

        top = (surfdz + (k - 1) * dz);
        bot = (surfdz + (k) * dz);

        T1 =
            (sw_visible *
             (surface_1 *
              exp(-param.LAKE_LAMWSW *
                  top) - surface_2 * exp(-param.LAKE_LAMWSW * bot)) +
             sw_nir *
             (surface_1 *
              exp(-param.LAKE_LAMWLW *
                  top) - surface_2 *
              exp(-param.LAKE_LAMWLW * bot))) / surface_avg;

        energyinput += T1 * surface_avg;
        energyinput /= surface[0];

        cnextra = 0.5 *
                  (-1. * surface_1 /
                   surface_avg) *
                  ((de[k - 1] / zhalf[k - 1]) * ((T[k] - T[k - 1]) / z[k]));

        *energy_out_bottom = 0.;
        *energy_out_bottom += surface_2 *
                              (sw_visible *
                               exp(-param.LAKE_LAMWSW *
                                   bot) + sw_nir *
                               exp(-param.LAKE_LAMWLW * bot));
        *energy_out_bottom /= surface[0];

        energymixed += cnextra;

        d[k] = T[k] +
               (T1 * dt) / ((1.e3 + water_density[k]) * cp[k] * z[k]) +
               cnextra * dt;

        /**********************************************************************
        * Calculate arrays for tridiagonal matrix.
        **********************************************************************/

        /* --------------------------------------------------------------------
         * Top node of the column.
         * --------------------------------------------------------------------*/


        surface_2 = surface[1];
        surface_avg = (surface[0] + surface[1]) / 2.;

        b[0] = -0.5 * (de[0] / zhalf[0]) *
               (dt / z[0]) * surface_2 / surface_avg;
        a[0] = 1. - b[0];

        /* --------------------------------------------------------------------
         * Second to second last node of the column.
         * --------------------------------------------------------------------*/

        for (k = 1; k < numnod - 1; k++) {
            surface_1 = surface[k];
            surface_2 = surface[k + 1];
            surface_avg = (surface[k] + surface[k + 1]) / 2.;

            b[k] = -0.5 * (de[k] / zhalf[k]) *
                   (dt / z[k]) * surface_2 / surface_avg;
            c[k] = -0.5 * (de[k - 1] / zhalf[k - 1]) *
                   (dt / z[k]) * surface_1 / surface_avg;
            a[k] = 1. - b[k] - c[k];
        }

        /* --------------------------------------------------------------------
         * Deepest node of the column.
         * --------------------------------------------------------------------*/

        surface_1 = surface[numnod - 1];
        surface_avg = surface[numnod - 1];
        c[numnod - 1] = -0.5 * (de[numnod - 1] / zhalf[numnod - 1]) *
                        (dt /
                         z[numnod - 1]) * surface_1 / surface_avg;
        a[numnod - 1] = 1. - c[numnod - 1];

        /**********************************************************************
        * Solve the tridiagonal matrix.
        **********************************************************************/

        tridia(numnod, c, a, b, d, Tnew);
    }

    /**********************************************************************
    * Adjust energy fluxes for change in density -> should be fixed by
    * moving to lagrangian scheme
    **********************************************************************/

    energycalc(Tnew, &joulenew, numnod, dz, surfdz, surface, cp, water_density);

    *temph = 0.0;

    *temph = joulenew;
}

/******************************************************************************
 * @brief    Simulate the convective mixing in the lake.
 *****************************************************************************/
void
tracer_mixer(double *T,
             int    *mixdepth,
             double *surface,
             int     numnod,
             double  dz,
             double  surfdz,
             double *cp)
{
    int    k, j, m;         /* Counter variables. */
    int    mixprev;
    double avet, avev;
    double heatcon; /* Heat content of the surface layer (formerly vol). */
    double Tav, densnew;
    double rho_max;
    double water_density[MAX_LAKE_NODES];

    for (k = 0; k < numnod; k++) {
        water_density[k] = calc_density(T[k]);
    }

/**********************************************************************
* Initialize the top depth of local instability.
**********************************************************************/

    mixprev = 0;

    for (k = 0; k < numnod - 1; k++) {
/**********************************************************************
* Check for instability at each slice in water column.
**********************************************************************/
        avet = 0.0;
        avev = 0.0;

        if (water_density[k] > water_density[k + 1]) {
/* --------------------------------------------------------------------
 * If there is instability apply the mixing scheme.
 * --------------------------------------------------------------------*/

            if (mixprev == 0 && (k + 1) > *mixdepth) {
/* ....................................................................
 * Correct the top depth of local instability.
 * ....................................................................*/
                *mixdepth = k + 1;
            }

/*----------------------------------------------------------------------
 * Mix from mixprev to k+1
 *----------------------------------------------------------------------*/

            for (m = mixprev; m <= k + 1; m++) {
/* --------------------------------------------------------------------
 * Apply the mixing scheme from the previous depth to the instability
 * up to this node.
 * --------------------------------------------------------------------*/

                if (m == 0) {
                    /* Calculate the heat content and volume of the surface layer. */
                    heatcon = surfdz *
                              (1.e3 + water_density[m]) * cp[m] * surface[m];
                }
                else {
                    /* Calculate the heat content and volume of all layers but
                       the surface layer. */
                    heatcon = dz *
                              (1.e3 + water_density[m]) * cp[m] * surface[m];
                }

                /* Calculate the volumetric weighed average lake temperature. */
                avet = avet + T[m] * heatcon;
                avev = avev + heatcon;
            }

            Tav = avet / avev;

/* --------------------------------------------------------------------
 * Calculate the density of the surface layer.
 * --------------------------------------------------------------------*/

            densnew = calc_density(Tav);

/* --------------------------------------------------------------------
 * Calculate the maximum density up to the local level.
 * --------------------------------------------------------------------*/

            rho_max = 0.0;

            for (j = 0; j < mixprev; j++) {
                if ((CONST_RHOFW + water_density[j]) > rho_max) {
                    rho_max = CONST_RHOFW + water_density[j];
                }
            }

/* --------------------------------------------------------------------
 * Adjust temperatures and density in the mixed part of column.
 * --------------------------------------------------------------------*/

            for (j = mixprev; j <= k + 1; j++) {
                T[j] = Tav;
                water_density[j] = densnew;
            }

/* --------------------------------------------------------------------
 * Check to make sure that the mixing has not generated new instabilities
 * above the previous depth to the instabilities.
 * --------------------------------------------------------------------*/

            if (rho_max > (CONST_RHOFW + densnew)) {
                /* If there are still instabilities iterate again..*/
                mixprev = 0;
                k = -1;
            }
        }
        else {
/**********************************************************************
* If there are no instabilities up to now then the depth to the
* instability has to be increased by 1 node.
**********************************************************************/

            mixprev = k + 1;
        }
    }

/**********************************************************************
* Recalculate the water density.
**********************************************************************/

    for (k = 0; k < numnod; k++) {
        water_density[k] = calc_density(T[k]);
    }
}

/******************************************************************************
 * @brief    Solve a tridiagonal system of equations.
 *
 * @note     History : Based on a streamlined version of the old NCAR ULIB
 *               subroutine TRDI used in the PHOENIX climate model of Schneider
 *               and Thompson (J.G.R., 1981). Revised by Starley Thompson to
 *               solve multiple systems and vectorize well on the CRAY-1. Later
 *               revised to include a PARAMETER statement to define loop limits
 *               and thus enable Cray short vector loops.
 *           Algorithm:  LU decomposition followed by solution.  NOTE: This
 *               subroutine executes satisfactorily if the input matrix is
 *               diagonally dominant and non-singular.  The diagonal elements
 *               are used to pivot, and no tests are made to determine
 *               singularity. If a singular or numerically singular matrix is
 *               used as input a divide by zero or floating point overflow will
 *               result.
 *****************************************************************************/
void
tridia(int     ne,
       double *a,
       double *b,
       double *c,
       double *y,
       double *x)
{
    double alpha[MAX_LAKE_NODES], gamma[MAX_LAKE_NODES];  /* Work arrays dimensioned (nd,ne).*/

    int    nm1, i;

    nm1 = ne - 1;

/**********************************************************************
* Obtain the LU decompositions.
**********************************************************************/

    alpha[0] = 1. / b[0];
    gamma[0] = c[0] * alpha[0];


    for (i = 1; i < nm1; i++) {
        alpha[i] = 1. / (b[i] - a[i] * gamma[i - 1]);
        gamma[i] = c[i] * alpha[i];
    }

/**********************************************************************
* Solve the system.
**********************************************************************/

    x[0] = y[0] * alpha[0];

    for (i = 1; i < nm1; i++) {
        x[i] = (y[i] - a[i] * x[i - 1]) * alpha[i];
    }


    x[nm1] = (y[nm1] - a[nm1] * x[nm1 - 1]) /
             (b[nm1] - a[nm1] * gamma[nm1 - 1]);

    for (i = nm1 - 1; i >= 0; i--) {
        x[i] = x[i] - gamma[i] * x[i + 1];
    }
}

/******************************************************************************
 * @brief    Calculate the thermal energy in the column.
 *****************************************************************************/
void
energycalc(double *finaltemp,
           double *sumjoule,
           int     numnod,
           double  dz,
           double  surfdz,
           double *surface,
           double *cp,
           double *density)
{
    double energy;
    int    k;

    *sumjoule = 0.0;

    for (k = 0; k < numnod; k++) {
        if (k == 0) {
            energy =
                (finaltemp[k] +
                 CONST_TKFRZ) * surfdz *
                (CONST_RHOFW +
                 density[k]) * cp[k] * (surface[k] + surface[k + 1]) / 2.;
        }
        else if (k == numnod - 1) {
            energy =
                (finaltemp[k] +
                 CONST_TKFRZ) * dz *
                (CONST_RHOFW + density[k]) * cp[k] * surface[k] / 2.;
        }
        else {
            energy =
                (finaltemp[k] +
                 CONST_TKFRZ) * dz *
                (CONST_RHOFW +
                 density[k]) * cp[k] * (surface[k] + surface[k + 1]) / 2.;
        }

        *sumjoule += energy;
    }
}

/******************************************************************************
 * @brief    This routine calculates the water balance of the lake.
 *****************************************************************************/
int
water_balance(lake_var_struct *lake,
              lake_con_struct  lake_con,
              double           dt,
              all_vars_struct *all_vars,
              int              iveg,
              int              band,
              double           lakefrac,
              soil_con_struct  soil_con,
              veg_con_struct   veg_con)
{
    extern option_struct       options;
    extern parameters_struct   param;
    extern global_param_struct global_param;

    int                        isave_n;
    double                     inflow_volume;
    double                     surfacearea, ldepth;
    double                     i_dbl;
    double                     index;
    size_t                     j, k, frost_area;
    double                     Tnew[MAX_LAKE_NODES];
    double                     circum;
    double                     baseflow_out_mm;
    double                     newfraction, Recharge;
    double                     abovegrnd_storage;
    int                        ErrorFlag;
    cell_data_struct         **cell;
    veg_var_struct           **veg_var;
    snow_data_struct         **snow;
    energy_bal_struct        **energy;
    size_t                     lindex;
    double                     frac;
    double                     Dsmax, resid_moist, liq, rel_moist;
    double                    *frost_fract;
    double                     volume_save;
    double                    *delta_moist = NULL;
    double                    *moist = NULL;
    double                     max_newfraction;

    cell = all_vars->cell;
    veg_var = all_vars->veg_var;
    snow = all_vars->snow;
    energy = all_vars->energy;

    frost_fract = soil_con.frost_fract;

    delta_moist = calloc(options.Nlayer, sizeof(*delta_moist));
    check_alloc_status(delta_moist, "Memory allocation error.");
    moist = calloc(options.Nlayer, sizeof(*moist));
    check_alloc_status(moist, "Memory allocation error.");

    /**********************************************************************
    * 1. Preliminary stuff
    **********************************************************************/

    isave_n = lake->activenod; /* save initial no. of nodes for later */

    inflow_volume = lake->runoff_in + lake->baseflow_in + lake->channel_in;

    /**********************************************************************
    * 2. calculate change in lake level for lake outflow calculation
    *     grid cell runoff (m3/TS for non-lake area)
    *     snow meltwater (mm)
    *     open water evaporation (mm)
    *     precip (added in solvelake)
    **********************************************************************/

    // Add runoff from rest of grid cell and wetland to lake, remove evaporation
    // (precip was added in solve_lake, to allow for snow interception)
    // Evaporation is not allowed to exceed the liquid volume of the lake (after incoming runoff & baseflow are added)
    if (fabs(lake->evapw) > DBL_EPSILON && lake->evapw >
        ((lake->volume -
          lake->ice_water_eq) + lake->ice_throughfall + inflow_volume +
         lake->snowmlt)) {
        lake->evapw =
            (lake->volume -
             lake->ice_water_eq) + lake->ice_throughfall + inflow_volume +
            lake->snowmlt;
        lake->volume = lake->ice_water_eq;
    }
    else {
        lake->volume +=
            (lake->ice_throughfall + inflow_volume + lake->snowmlt -
             lake->evapw);
    }

    // Estimate new surface area of liquid water for recharge calculations
    volume_save = lake->volume;
    ErrorFlag = get_depth(lake_con, lake->volume - lake->ice_water_eq, &ldepth);
    if (ErrorFlag == ERROR) {
        log_err("Error calculating depth: volume = %f, depth = %e",
                lake->volume, ldepth);
    }
    ErrorFlag = get_sarea(lake_con, ldepth, &surfacearea);
    if (ErrorFlag == ERROR) {
        log_err("Error calculating area: depth = %f, sarea = %e",
                ldepth, surfacearea);
    }

    // Estimate the new lake fraction (before recharge)
    if (lake->new_ice_area > surfacearea) {
        surfacearea = lake->new_ice_area;
    }
    newfraction = surfacearea / lake_con.basin[0];

    // Save this estimate of the new lake fraction for use later
    max_newfraction = newfraction;

    /**********************************************************************
    * 3. calculate recharge to wetland
    **********************************************************************/

    // Based on the above initial estimate of lake area, the lake
    // will either inundate some of the wetland or will recede and
    // expose new wetland.  In the case of inundation, we will
    // take all above-ground moisture storage (snowpack, etc) and
    // give it to the lake, but at the same time we will take water
    // from the lake and fill the inundated soil to saturation.
    // In the case of a receding lake, newly-exposed wetland is
    // assumed to have saturated soil and 0 above-ground storage.

    // The redistribution of moisture within the wetland will happen
    // during the final step of this function.

    // Note:
    // This does not account for phase changes and temperature
    // changes, so will generate some energy balance errors.
    // In addition, the moisture exchange with the wetland
    // will make our initial estimate of the new lake area
    // somewhat inaccurate (the initial estimate will be an
    // upper bound on the new lake area, resulting in some of
    // the wetland being "splashed" by the lake, causing the
    // snow to run off and the soil to moisten, but not being
    // inundated).  However, lake dimensions will be recalculated
    // after runoff and baseflow are subtracted from the lake.

    lake->recharge = 0.0;
    for (j = 0; j < options.Nlayer; j++) {
        delta_moist[j] = 0; // mm over (1-lakefrac)
    }

    if (max_newfraction > lakefrac) {
        // Lake must fill soil to saturation in the newly-flooded area
        for (j = 0; j < options.Nlayer; j++) {
            delta_moist[j] +=
                (soil_con.max_moist[j] -
                 cell[iveg][band].layer[j].moist) *
                (max_newfraction - lakefrac) / (1 - lakefrac);                                                           // mm over (1-lakefrac)
        }
        for (j = 0; j < options.Nlayer; j++) {
            lake->recharge += (delta_moist[j]) / MM_PER_M *
                              (1 - lakefrac) * lake_con.basin[0];                    // m^3
        }

        // Above-ground storage in newly-flooded area is liberated and goes to lake
        abovegrnd_storage =
            (veg_var[iveg][band].Wdew / MM_PER_M +
             snow[iveg][band].snow_canopy +
             snow[iveg][band].swq) *
            (max_newfraction - lakefrac) * lake_con.basin[0];
        lake->recharge -= abovegrnd_storage;

        // Fill the soil to saturation if possible in inundated area
        // Subtract the recharge (which may be negative) from the lake
        if (lake->volume - lake->ice_water_eq > lake->recharge) { // enough liquid water to support recharge
            lake->volume -= lake->recharge;
        }
        else { // not enough liquid water to support recharge; fill soil as much as allowed by available liquid water in lake and above-ground storage in newly-flooded area; lake will recede back from this point after recharge is taken out of it
            lake->recharge = lake->volume - lake->ice_water_eq;
            lake->volume = lake->ice_water_eq;

            Recharge = MM_PER_M * lake->recharge /
                       ((max_newfraction -
                         lakefrac) *
                        lake_con.basin[0]) +
                       (veg_var[iveg][band].Wdew +
                        snow[iveg][band].snow_canopy * MM_PER_M +
                        snow[iveg][band].swq * MM_PER_M);                                                                                                                               // mm over area that has been flooded

            for (j = 0; j < options.Nlayer; j++) {
                if (Recharge >
                    (soil_con.max_moist[j] - cell[iveg][band].layer[j].moist)) {
                    Recharge -=
                        (soil_con.max_moist[j] -
                         cell[iveg][band].layer[j].moist);
                    delta_moist[j] =
                        (soil_con.max_moist[j] -
                         cell[iveg][band].layer[j].moist) *
                        (max_newfraction - lakefrac) / (1 - lakefrac);                                                      // mm over (1-lakefrac)
                }
                else {
                    delta_moist[j] = Recharge *
                                     (max_newfraction -
                                      lakefrac) / (1 - lakefrac);            // mm over (1-lakefrac)
                    Recharge = 0.0;
                }
            }
        }
    }

    /**********************************************************************
    * 4. Calculate outflow from lake.  Runoff estimate is based on the
    *    equation for flow over a broad crested weir.  Baseflow estimate
    *    is the ARNO formulation, using the ice content of the adjacent
    *    wetland.  Outgoing runoff and baseflow are in m3.
    **********************************************************************/

    Dsmax = soil_con.Dsmax / global_param.model_steps_per_day;
    lindex = options.Nlayer - 1;
    liq = 0;
    for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
        liq +=
            (soil_con.max_moist[lindex] -
             cell[iveg][band].layer[lindex].ice[frost_area]) *
            frost_fract[frost_area];
    }
    resid_moist = soil_con.resid_moist[lindex] * soil_con.depth[lindex] *
                  MM_PER_M;

    /** Compute relative moisture **/
    rel_moist =
        (liq - resid_moist) / (soil_con.max_moist[lindex] - resid_moist);

    /** Compute baseflow as function of relative moisture **/
    frac = Dsmax * soil_con.Ds / soil_con.Ws;
    baseflow_out_mm = frac * rel_moist;
    if (rel_moist > soil_con.Ws) {
        frac = (rel_moist - soil_con.Ws) / (1 - soil_con.Ws);
        baseflow_out_mm += Dsmax * (1 - soil_con.Ds / soil_con.Ws) * pow(frac,
                                                                         soil_con.c);
    }
    if (baseflow_out_mm < 0) {
        baseflow_out_mm = 0;
    }

    // extract baseflow volume from the lake m^3
    // baseflow will only come from under the liquid portion of the lake
    ErrorFlag = get_depth(lake_con, lake->volume - lake->ice_water_eq, &ldepth);
    if (ErrorFlag == ERROR) {
        log_err("Error calculating depth: volume = %f, depth = %e",
                lake->volume, ldepth);
    }
    ErrorFlag = get_sarea(lake_con, ldepth, &surfacearea);
    if (ErrorFlag == ERROR) {
        log_err("Error calculating area: depth = %f, sarea = %e",
                ldepth, surfacearea);
    }
    lake->baseflow_out = baseflow_out_mm * surfacearea / MM_PER_M;
    if (lake->volume - lake->ice_water_eq >= lake->baseflow_out) {
        lake->volume -= lake->baseflow_out;
    }
    else {
        lake->baseflow_out = lake->volume - lake->ice_water_eq;
        lake->volume -= lake->baseflow_out;
    }

    // Find new lake depth for runoff calculations
    ErrorFlag = get_depth(lake_con, lake->volume - lake->ice_water_eq, &ldepth);
    if (ErrorFlag == ERROR) {
        log_err("Error calculating depth: volume = %f, depth = %e",
                lake->volume, ldepth);
    }

    // Compute runoff volume in m^3 and extract runoff volume from lake
    if (ldepth <= lake_con.mindepth) {
        lake->runoff_out = 0.0;
    }
    else {
        circum = 2 * CONST_PI * pow(surfacearea / CONST_PI, 0.5);
        lake->runoff_out = lake_con.wfrac * circum * dt *
                           1.6 * pow(ldepth - lake_con.mindepth, 1.5);
        if ((lake->volume - lake->ice_water_eq) >= lake->runoff_out) {
            /*liquid water is available */
            if ((lake->volume - lake->runoff_out) < lake_con.minvolume) {
                lake->runoff_out = lake->volume - lake_con.minvolume;
            }
            lake->volume -= lake->runoff_out;
        }
        else {
            lake->runoff_out = lake->volume - lake->ice_water_eq;
            if ((lake->volume - lake->runoff_out) < lake_con.minvolume) {
                lake->runoff_out = lake->volume - lake_con.minvolume;
            }
            lake->volume -= lake->runoff_out;
        }
    }

    // Check that lake volume does not exceed our earlier estimate.
    // This will prevent runaway lake growth for the case in which
    // the lake recharge to wetland is negative and large enough to
    // more than compensate for runoff and baseflow out of the lake.
    if (lake->volume > volume_save) {
        lake->runoff_out += lake->volume - volume_save;
        lake->volume = volume_save;
    }

    // check that lake volume does not exceed its maximum
    if (lake->volume - lake_con.maxvolume > DBL_EPSILON) {
        if (lake->ice_water_eq > lake_con.maxvolume) {
            lake->runoff_out += (lake->volume - lake->ice_water_eq);
            lake->volume = lake->ice_water_eq;
        }
        else {
            lake->runoff_out += (lake->volume - lake_con.maxvolume);
            lake->volume = lake_con.maxvolume;
        }
    }
    else if (lake->volume < DBL_EPSILON) {
        lake->volume = 0.0;
    }

    /**********************************************************************/
    /* End of runoff calculation */
    /**********************************************************************/

    // Recalculate lake depth to define surface[] for next time step
    // Here, we only want depth of liquid water (below ice bottom), since surface[] array only applies to liquid water
    ErrorFlag =
        get_depth(lake_con, lake->volume - lake->ice_water_eq, &(lake->ldepth));
    if (ErrorFlag == ERROR) {
        log_err("Error calculating depth: volume = %f, depth = %e",
                lake->volume, lake->ldepth);
    }

    /**********************************************************************
    *  4. Adjust the activenodes and lake area array.
    **********************************************************************/

    compute_derived_lake_dimensions(lake, lake_con);

    // Final lakefraction
    if (lake->new_ice_area > lake->surface[0]) {
        lake->sarea = lake->new_ice_area;
    }
    else {
        lake->sarea = lake->surface[0];
    }
    newfraction = lake->sarea / lake_con.basin[0];

    /*******************************************************************/

    /* Adjust temperature distribution if number of nodes has changed.
       Note:  This approach (and the lake model in general) does not preserve
       the thermal energy of the water column. */

    index = 0.;
    if (lake->activenod != isave_n) {
        for (k = 0; k < lake->activenod; k++) {
            Tnew[k] = 0.0;
            for (i_dbl = 0; i_dbl < isave_n; i_dbl++) {
                index += (1. / lake->activenod);
                Tnew[k] += lake->temp[(int) floor(index)];
            }
        }
        for (k = 0; k < lake->activenod; k++) {
            if (isave_n > 0) {
                lake->temp[k] = Tnew[k] / isave_n;
            }
            else {
                lake->temp[k] = all_vars->energy[iveg][band].Tsurf;
            }
        }
    }

    if (lake->activenod == isave_n && isave_n == 0) {
        lake->temp[k] = all_vars->energy[iveg][band].Tsurf;
    }

    /**********************************************************************
        5. Rescale the fluxes in the lake and the wetland by the change in lake area;
           Advect the storages.
           Final units of storages and fluxes should be in mm/final lake or wetland area.
    **********************************************************************/
    // Wetland
    if (newfraction < 1.0) { // wetland exists at end of time step
        advect_soil_veg_storage(lakefrac, max_newfraction, newfraction,
                                delta_moist, &soil_con, &veg_con,
                                &(cell[iveg][band]), &(veg_var[iveg][band]),
                                lake_con);
        rescale_soil_veg_fluxes((1 - lakefrac), (1 - newfraction),
                                &(cell[iveg][band]), &(veg_var[iveg][band]));
        advect_snow_storage(lakefrac, max_newfraction, newfraction,
                            &(snow[iveg][band]));
        rescale_snow_energy_fluxes((1 - lakefrac), (1 - newfraction),
                                   &(snow[iveg][band]), &(energy[iveg][band]));
        for (j = 0; j < options.Nlayer; j++) {
            moist[j] = cell[iveg][band].layer[j].moist;
        }
        ErrorFlag = distribute_node_moisture_properties(
            energy[iveg][band].moist, energy[iveg][band].ice,
            energy[iveg][band].kappa_node,
            energy[iveg][band].Cs_node,
            soil_con.Zsum_node,
            energy[iveg][band].T,
            soil_con.max_moist_node,
            soil_con.expt_node,
            soil_con.bubble_node,
            moist, soil_con.depth,
            soil_con.soil_dens_min,
            soil_con.bulk_dens_min,
            soil_con.quartz,
            soil_con.soil_density,
            soil_con.bulk_density,
            soil_con.organic, options.Nnode,
            options.Nlayer,
            soil_con.FS_ACTIVE);
        if (ErrorFlag == ERROR) {
            return (ERROR);
        }
    }
    else if (lakefrac < 1.0) { // wetland is gone at end of time step, but existed at beginning of step
        if (lakefrac > 0.0) { // lake also existed at beginning of step
            for (j = 0; j < options.Nlayer; j++) {
                lake->evapw += cell[iveg][band].layer[j].evap / MM_PER_M *
                               (1. - lakefrac) * lake_con.basin[0];
            }
            lake->evapw += veg_var[iveg][band].canopyevap / MM_PER_M *
                           (1. - lakefrac) * lake_con.basin[0];
            lake->evapw += snow[iveg][band].canopy_vapor_flux *
                           (1. - lakefrac) * lake_con.basin[0];
            lake->evapw += snow[iveg][band].vapor_flux *
                           (1. - lakefrac) * lake_con.basin[0];
        }
    }

    // Lake
    if (newfraction > 0.0) { // lake exists at end of time step
        // Copy moisture fluxes into lake->soil structure, mm over end-of-step lake area
        lake->soil.runoff = lake->runoff_out * MM_PER_M /
                            (newfraction * lake_con.basin[0]);
        lake->soil.baseflow = lake->baseflow_out * MM_PER_M /
                              (newfraction * lake_con.basin[0]);
        lake->soil.inflow = lake->baseflow_out * MM_PER_M /
                            (newfraction * lake_con.basin[0]);
        for (lindex = 0; lindex < options.Nlayer; lindex++) {
            lake->soil.layer[lindex].evap = 0;
        }
        lake->soil.layer[0].evap += lake->evapw * MM_PER_M /
                                    (newfraction * lake_con.basin[0]);
        // Rescale other fluxes and storages to mm over end-of-step lake area
        if (lakefrac > 0.0) { // lake existed at beginning of time step
            rescale_snow_storage(lakefrac, newfraction, &(lake->snow));
            rescale_snow_energy_fluxes(lakefrac, newfraction, &(lake->snow),
                                       &(lake->energy));
            if (lake->snow.swq > 0) {
                lake->snow.coverage = lake->new_ice_area / lake->sarea;
            }
            else {
                lake->snow.coverage = 0;
            }
        }
        else { // lake didn't exist at beginning of time step; create new lake
               // Reset all non-essential lake variables
            initialize_lake(lake, lake_con, &soil_con, &(cell[iveg][band]),
                            true);

            // Compute lake dimensions from current lake depth
            compute_derived_lake_dimensions(lake, lake_con);

            // Assign initial temperatures
            for (k = 0; k < lake->activenod; k++) {
                lake->temp[k] = all_vars->energy[iveg][band].T[0];
            }
            lake->tempavg = lake->temp[0];
            lake->energy.Tsurf = all_vars->energy[iveg][band].Tsurf;
            for (k = 0; k < options.Nnode; k++) {
                lake->energy.T[k] = all_vars->energy[iveg][band].T[k];
            }
            for (k = 0; k < options.Nlayer; k++) {
                lake->soil.layer[k].T = all_vars->cell[iveg][band].layer[k].T;
            }
        }
    }
    else if (lakefrac > 0.0) { // lake is gone at end of time step, but existed at beginning of step
        if (lakefrac < 1.0) { // wetland also existed at beginning of step
            cell[iveg][band].layer[0].evap += MM_PER_M * lake->evapw /
                                              ((1. -
                                                newfraction) *
                                               lake_con.basin[0]);
            cell[iveg][band].runoff += MM_PER_M * lake->runoff_out /
                                       ((1. - newfraction) * lake_con.basin[0]);
            cell[iveg][band].baseflow += MM_PER_M * lake->baseflow_out /
                                         ((1. -
                                           newfraction) * lake_con.basin[0]);
            cell[iveg][band].inflow += MM_PER_M * lake->baseflow_out /
                                       ((1. - newfraction) * lake_con.basin[0]);
        }
    }

    if (options.CARBON) {
        advect_carbon_storage(lakefrac, newfraction, lake, &(cell[iveg][band]));
    }

    free((char*) delta_moist);
    free((char*) moist);

    return(0);
}

/******************************************************************************
 * @brief    Function to update moisture storage in the wetland soil and veg to
 *           account for changes in wetland area.
 *****************************************************************************/
void
advect_soil_veg_storage(double            lakefrac,
                        double            max_newfraction,
                        double            newfraction,
                        double           *delta_moist,
                        soil_con_struct  *soil_con,
                        veg_con_struct   *veg_con,
                        cell_data_struct *cell,
                        veg_var_struct   *veg_var,
                        lake_con_struct   lake_con)
{
    extern option_struct options;
    int                  ilidx;
    size_t               lidx;
    size_t               fidx;
    double               new_moist[MAX_LAYERS];
    double               tmp_moist[MAX_LAYERS];
    double               tmp_runoff;

    if (lakefrac < 1.0) { // wetland existed during this step
        // Add delta_moist to wetland, using wetland's initial area (1-lakefrac)
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            new_moist[lidx] = cell->layer[lidx].moist + delta_moist[lidx]; // mm over (1-lakefrac)
            delta_moist[lidx] = 0;
            if (new_moist[lidx] > soil_con->max_moist[lidx]) {
                if (lidx < options.Nlayer - 1) {
                    delta_moist[lidx +
                                1] += new_moist[lidx] -
                                      soil_con->max_moist[lidx];
                }
                else {
                    delta_moist[lidx] += new_moist[lidx] -
                                         soil_con->max_moist[lidx];
                }
                new_moist[lidx] = soil_con->max_moist[lidx];
            }
        }
        for (ilidx = (int) options.Nlayer - 1; ilidx >= 0; ilidx--) {
            new_moist[ilidx] += delta_moist[ilidx]; // mm over (1-lakefrac)
            delta_moist[ilidx] = 0;
            if (new_moist[ilidx] > soil_con->max_moist[ilidx]) {
                if (ilidx > 0) {
                    delta_moist[ilidx -
                                1] += new_moist[ilidx] -
                                      soil_con->max_moist[ilidx];
                }
                else {
                    delta_moist[ilidx] += new_moist[ilidx] -
                                          soil_con->max_moist[ilidx];
                }
                new_moist[ilidx] = soil_con->max_moist[ilidx];
            }
        }

        // Any recharge that cannot be accomodated by wetland goes to baseflow
        if (delta_moist[0] > 0) {
            cell->baseflow += delta_moist[0] / MM_PER_M *
                              (1 - lakefrac) * lake_con.basin[0];            // m^3
            delta_moist[0] = 0;
        }

        // Rescale wetland moisture to wetland's final area (= 1-newfraction)
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            new_moist[lidx] *= (1 - lakefrac); // mm over lake/wetland tile
            new_moist[lidx] += soil_con->max_moist[lidx] *
                               (lakefrac - newfraction);                   // Add the saturated portion between lakefrac and newfraction; this works whether newfraction is > or < or == lakefrac
            new_moist[lidx] /= (1 - newfraction); // mm over final wetland area (1-newfraction)
            cell->layer[lidx].moist = new_moist[lidx]; // mm over (1-newfraction)
        }

        // Recompute saturated areas
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            tmp_moist[lidx] = cell->layer[lidx].moist;
        }
        compute_runoff_and_asat(soil_con, tmp_moist, 0, &(cell->asat),
                                &tmp_runoff);

        // Recompute zwt's
        wrap_compute_zwt(soil_con, cell);

        // Update Wdew
        if (max_newfraction <= lakefrac) { // lake only receded
            if (veg_var != NULL) {
                veg_var->Wdew *= (1 - lakefrac) / (1 - newfraction);
            }
        }
        else {
            if (veg_var != NULL) {
                veg_var->Wdew *= (1 - max_newfraction) / (1 - newfraction);
            }
        }
    }
    else { // Wetland didn't exist until now; create new wetland
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            cell->layer[lidx].moist = soil_con->max_moist[lidx];
            for (fidx = 0; fidx < options.Nfrost; fidx++) {
                cell->layer[lidx].ice[fidx] = 0.0;
            }
        }
        cell->asat = 1.0;
        cell->zwt = 0;
        cell->zwt_lumped = 0;

        if (veg_var != NULL) {
            veg_var->Wdew = 0.0;
        }
    }

    // Compute rootmoist and wetness
    cell->rootmoist = 0;
    cell->wetness = 0;
    for (lidx = 0; lidx < options.Nlayer; lidx++) {
        if (veg_con->root[lidx] > 0) {
            cell->rootmoist += cell->layer[lidx].moist;
        }
        cell->wetness +=
            (cell->layer[lidx].moist -
             soil_con->Wpwp[lidx]) /
            (soil_con->porosity[lidx] * soil_con->depth[lidx] * MM_PER_M -
             soil_con->Wpwp[lidx]);
    }
    cell->wetness /= options.Nlayer;
}

/******************************************************************************
 * @brief    Function to update wetland soil and veg moisture fluxes to account
 *           for changes in wetland area.
 *****************************************************************************/
void
rescale_soil_veg_fluxes(double            oldfrac,
                        double            newfrac,
                        cell_data_struct *cell,
                        veg_var_struct   *veg_var)
{
    extern option_struct options;
    size_t               lidx;

    if (newfrac < DBL_EPSILON) {
        newfrac = DBL_EPSILON;
    }

    if (oldfrac > 0.0) { // existed at beginning of time step
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            cell->layer[lidx].evap *= oldfrac / newfrac;
        }
        cell->baseflow *= oldfrac / newfrac;
        cell->inflow *= oldfrac / newfrac;
        cell->runoff *= oldfrac / newfrac;
        if (veg_var != NULL) {
            veg_var->canopyevap *= oldfrac / newfrac;
            veg_var->throughfall *= oldfrac / newfrac;
        }
    }
    else { // didn't exist at beginning of time step; set fluxes to 0
        for (lidx = 0; lidx < options.Nlayer; lidx++) {
            cell->layer[lidx].evap = 0.0;
        }
        cell->baseflow = 0.0;
        cell->inflow = 0.0;
        cell->runoff = 0.0;
        if (veg_var != NULL) {
            veg_var->canopyevap = 0.0;
            veg_var->throughfall = 0.0;
        }
    }
}

/******************************************************************************
 * @brief    Function to update moisture storage in the wetland snow pack to
 *           account for changes in wetland area.
 *****************************************************************************/
void
advect_snow_storage(double            lakefrac,
                    double            max_newfraction,
                    double            newfraction,
                    snow_data_struct *snow)
{
    if ((1 - newfraction) < DBL_EPSILON) {
        newfraction = 1 - DBL_EPSILON;
    }

    if (lakefrac < 1.0) { // wetland existed during this step
        if (max_newfraction > lakefrac) {
            snow->depth *= (1 - max_newfraction) / (1 - newfraction);
            snow->pack_water *= (1 - max_newfraction) / (1 - newfraction);
            snow->snow_canopy *= (1 - max_newfraction) / (1 - newfraction);
            snow->surf_water *= (1 - max_newfraction) / (1 - newfraction);
            snow->swq *= (1 - max_newfraction) / (1 - newfraction);
        }
        else {
            snow->depth *= (1 - lakefrac) / (1 - newfraction);
            snow->pack_water *= (1 - lakefrac) / (1 - newfraction);
            snow->snow_canopy *= (1 - lakefrac) / (1 - newfraction);
            snow->surf_water *= (1 - lakefrac) / (1 - newfraction);
            snow->swq *= (1 - lakefrac) / (1 - newfraction);
        }
    }
    else { // Wetland didn't exist until now; create new wetland
        snow->depth = 0.0;
        snow->pack_water = 0.0;
        snow->snow_canopy = 0.0;
        snow->surf_water = 0.0;
        snow->swq = 0.0;
    }
}

/******************************************************************************
 * @brief    Function to update moisture storage in the lake ice snowpack to
 *           account for changes in lake area.
 *****************************************************************************/
void
rescale_snow_storage(double            oldfrac,
                     double            newfrac,
                     snow_data_struct *snow)
{
    if (newfrac < DBL_EPSILON) {
        newfrac = DBL_EPSILON;
    }

    snow->depth *= oldfrac / newfrac;
    snow->pack_water *= oldfrac / newfrac;
    snow->snow_canopy *= oldfrac / newfrac;
    snow->surf_water *= oldfrac / newfrac;
    snow->swq *= oldfrac / newfrac;
}

/******************************************************************************
 * @brief    Function to update the wetland snow and energy fluxes to account
 *           for changes in wetland area.
 *****************************************************************************/
void
rescale_snow_energy_fluxes(double             oldfrac,
                           double             newfrac,
                           snow_data_struct  *snow,
                           energy_bal_struct *energy)
{
    if (newfrac < DBL_EPSILON) {
        newfrac = DBL_EPSILON;
    }

    if (oldfrac > 0.0) { // existed at beginning of time step
        snow->blowing_flux *= oldfrac / newfrac;
        snow->melt *= oldfrac / newfrac;
        snow->surface_flux *= oldfrac / newfrac;
        snow->vapor_flux *= oldfrac / newfrac;
        energy->advected_sensible *= oldfrac / newfrac;
        energy->advection *= oldfrac / newfrac;
        energy->AtmosError *= oldfrac / newfrac;
        energy->AtmosLatent *= oldfrac / newfrac;
        energy->AtmosLatentSub *= oldfrac / newfrac;
        energy->AtmosSensible *= oldfrac / newfrac;
        energy->canopy_advection *= oldfrac / newfrac;
        energy->canopy_latent *= oldfrac / newfrac;
        energy->canopy_latent_sub *= oldfrac / newfrac;
        energy->canopy_refreeze *= oldfrac / newfrac;
        energy->canopy_sensible *= oldfrac / newfrac;
        energy->deltaCC *= oldfrac / newfrac;
        energy->deltaH *= oldfrac / newfrac;
        energy->error *= oldfrac / newfrac;
        energy->fusion *= oldfrac / newfrac;
        energy->grnd_flux *= oldfrac / newfrac;
        energy->latent *= oldfrac / newfrac;
        energy->latent_sub *= oldfrac / newfrac;
        energy->longwave *= oldfrac / newfrac;
        energy->LongOverIn *= oldfrac / newfrac;
        energy->LongUnderIn *= oldfrac / newfrac;
        energy->LongUnderOut *= oldfrac / newfrac;
        energy->melt_energy *= oldfrac / newfrac;
        energy->NetLongAtmos *= oldfrac / newfrac;
        energy->NetLongOver *= oldfrac / newfrac;
        energy->NetLongUnder *= oldfrac / newfrac;
        energy->NetShortAtmos *= oldfrac / newfrac;
        energy->NetShortGrnd *= oldfrac / newfrac;
        energy->NetShortOver *= oldfrac / newfrac;
        energy->NetShortUnder *= oldfrac / newfrac;
        energy->out_long_canopy *= oldfrac / newfrac;
        energy->out_long_surface *= oldfrac / newfrac;
        energy->refreeze_energy *= oldfrac / newfrac;
        energy->sensible *= oldfrac / newfrac;
        energy->shortwave *= oldfrac / newfrac;
        energy->ShortOverIn *= oldfrac / newfrac;
        energy->ShortUnderIn *= oldfrac / newfrac;
        energy->snow_flux *= oldfrac / newfrac;
    }
    else { // didn't exist at beginning of time step; set fluxes to 0
        snow->blowing_flux = 0.0;
        snow->melt = 0.0;
        snow->surface_flux = 0.0;
        snow->vapor_flux = 0.0;
        energy->advected_sensible = 0.0;
        energy->advection = 0.0;
        energy->AtmosError = 0.0;
        energy->AtmosLatent = 0.0;
        energy->AtmosLatentSub = 0.0;
        energy->AtmosSensible = 0.0;
        energy->canopy_advection = 0.0;
        energy->canopy_latent = 0.0;
        energy->canopy_latent_sub = 0.0;
        energy->canopy_refreeze = 0.0;
        energy->canopy_sensible = 0.0;
        energy->deltaCC = 0.0;
        energy->deltaH = 0.0;
        energy->error = 0.0;
        energy->fusion = 0.0;
        energy->grnd_flux = 0.0;
        energy->latent = 0.0;
        energy->latent_sub = 0.0;
        energy->longwave = 0.0;
        energy->LongOverIn = 0.0;
        energy->LongUnderIn = 0.0;
        energy->LongUnderOut = 0.0;
        energy->melt_energy = 0.0;
        energy->NetLongAtmos = 0.0;
        energy->NetLongOver = 0.0;
        energy->NetLongUnder = 0.0;
        energy->NetShortAtmos = 0.0;
        energy->NetShortGrnd = 0.0;
        energy->NetShortOver = 0.0;
        energy->NetShortUnder = 0.0;
        energy->out_long_canopy = 0.0;
        energy->out_long_surface = 0.0;
        energy->refreeze_energy = 0.0;
        energy->sensible = 0.0;
        energy->shortwave = 0.0;
        energy->ShortOverIn = 0.0;
        energy->ShortUnderIn = 0.0;
        energy->snow_flux = 0.0;
    }
}

/******************************************************************************
 * @brief    Function to update carbon storage in the lake and wetland soil
 *           columns to account for changes in wetland area.
 *****************************************************************************/
void
advect_carbon_storage(double            lakefrac,
                      double            newfraction,
                      lake_var_struct  *lake,
                      cell_data_struct *cell)
{
    if (newfraction > lakefrac) { // lake grew, wetland shrank
        if (newfraction < DBL_EPSILON) {
            newfraction = DBL_EPSILON;
        }
        lake->soil.CLitter =
            (lakefrac * lake->soil.CLitter +
             (newfraction - lakefrac) * cell->CLitter) / newfraction;
        lake->soil.CInter =
            (lakefrac * lake->soil.CInter +
             (newfraction - lakefrac) * cell->CInter) / newfraction;
        lake->soil.CSlow =
            (lakefrac * lake->soil.CSlow +
             (newfraction - lakefrac) * cell->CSlow) / newfraction;
    }
    else if (newfraction < lakefrac) { // lake shrank, wetland grew
        if ((1 - newfraction) < DBL_EPSILON) {
            newfraction = 1 - DBL_EPSILON;
        }
        cell->CLitter =
            ((lakefrac -
              newfraction) * lake->soil.CLitter +
             (1 - lakefrac) * cell->CLitter) / (1 - newfraction);
        cell->CInter =
            ((lakefrac -
              newfraction) * lake->soil.CInter +
             (1 - lakefrac) * cell->CInter) / (1 - newfraction);
        cell->CSlow =
            ((lakefrac -
              newfraction) * lake->soil.CSlow +
             (1 - lakefrac) * cell->CSlow) / (1 - newfraction);
    }
}
