/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine computes the surface energy balance for bare soil and
 * vegetation uncovered by snow.  It computes outgoing longwave, sensible heat
 * flux, ground heat flux, and storage of heat in the thin upper layer, based
 * on the given surface temperature.
 *
 * The Energy Balance Equation used comes from Xu Liang's Paper "Insights of
 * the Ground Heat Flux in Land Surface Parameterization Schemes."
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
 * @brief    Calculate the surface energy balance.
 *****************************************************************************/
double
func_surf_energy_bal(double  Ts,
                     va_list ap)
{
    extern parameters_struct param;
    extern option_struct     options;

    /* define routine input variables */

    /* general model terms */
    size_t             i;
    int                VEG;
    int                veg_class;
    int                Error;

    double             delta_t;

    /* soil layer terms */
    double             Cs1;
    double             Cs2;
    double             D1;
    double             D2;
    double             T1_old;
    double             T2;
    double             Ts_old;
    double            *Told_node;
    double             b_infilt;
    double             bubble;
    double             dp;
    double             expt;
    double             ice0;
    double             kappa1;
    double             kappa2;
    double             max_moist;
    double             moist;

    double            *Wmax;
    double            *Wcr;
    double            *Wpwp;
    double            *depth;
    double            *resid_moist;
    double            *bulk_dens_min;
    double            *soil_dens_min;
    double            *quartz;
    double            *bulk_density;
    double            *soil_density;
    double            *organic;

    double            *root;
    double            *CanopLayerBnd;

    /* meteorological forcing terms */
    int                UnderStory;
    int                overstory;

    double             NetShortBare; // net SW that reaches bare ground
    double             NetShortGrnd; // net SW that penetrates snowpack
    double             NetShortSnow; // net SW that reaches snow surface
    double             Tair; // temperature of canopy air or atmosphere
    double             atmos_density;
    double             atmos_pressure;
    double             elevation;
    double             emissivity;
    double             LongBareIn; // incoming LW to snow-free surface
    double             LongSnowIn; // incoming LW to snow surface - if INCLUDE_SNOW
    double             surf_atten;
    double             vp;
    double             vpd;
    double             shortwave;
    double             Catm;
    double            *dryFrac;

    double            *Wdew;
    double            *displacement;
    double            *ra;
    double            *Ra_veg;
    double            *Ra_used;
    double             rainfall;
    double            *ref_height;
    double            *roughness;
    double            *wind;

    /* latent heat terms */
    double             Le;

    /* snowpack terms */
    double             Advection;
    double             OldTSurf;
    double             Tsnow_surf;
    double             kappa_snow; // snow conductance / depth
    double             melt_energy; // energy consumed in reducing the snowpack coverage
    double             snow_coverage; // snowpack coverage fraction
    double             snow_density;
    double             snow_swq;
    double             snow_water;

    double            *deltaCC;
    double            *refreeze_energy;
    double            *vapor_flux;
    double            *blowing_flux;
    double            *surface_flux;

    /* soil node terms */
    int                Nnodes;

    double            *Cs_node;
    double            *T_node;
    double            *Tnew_node;
    char              *Tnew_fbflag;
    unsigned          *Tnew_fbcount;
    double            *alpha;
    double            *beta;
    double            *bubble_node;
    double            *Zsum_node;
    double            *expt_node;
    double            *gamma;
    double            *ice_node;
    double            *kappa_node;
    double            *max_moist_node;
    double            *moist_node;

    /* spatial frost terms */
    double            *frost_fract;

    /* model structures */
    soil_con_struct   *soil_con;
    layer_data_struct *layer;
    veg_var_struct    *veg_var;

    /* control flags */
    int                INCLUDE_SNOW;
    int                FS_ACTIVE;
    int                NOFLUX;
    int                EXP_TRANS;
    int                SNOWING;

    int               *FIRST_SOLN;

    /* returned energy balance terms */
    double            *NetLongBare; // net LW from snow-free ground
    double            *NetLongSnow; // net longwave from snow surface - if INCLUDE_SNOW
    double            *T1;
    double            *deltaH;
    double            *fusion;
    double            *grnd_flux;
    double            *latent_heat;
    double            *latent_heat_sub;
    double            *sensible_heat;
    double            *snow_flux;
    double            *store_error;

    /* Define internal routine variables */
    double             Evap; /** Total evap in m/s **/
    double             LongBareOut; // outgoing LW from snow-free ground
    double             NetBareRad;
    double             TMean;
    double             error;
    double             ice;
    double             temp_latent_heat;
    double             temp_latent_heat_sub;
    double             VaporMassFlux;
    double             BlowingMassFlux;
    double             SurfaceMassFlux;
    double             T11;
    double             T1_minus;
    double             T1_plus;
    double             D1_minus;
    double             D1_plus;
    double            *transp;
    double             Ra_bare[3];
    double             tmp_wind[3];
    double             tmp_height;
    double             tmp_displacement[3];
    double             tmp_roughness[3];
    double             tmp_ref_height[3];
    double             ga_veg;
    double             ga_bare;
    double             ga_average;

    /************************************
       Read variables from variable list
    ************************************/

    /* general model terms */
    VEG = (int) va_arg(ap, int);
    veg_class = (int) va_arg(ap, int);
    delta_t = (double) va_arg(ap, double);

    /* soil layer terms */
    Cs1 = (double) va_arg(ap, double);
    Cs2 = (double) va_arg(ap, double);
    D1 = (double) va_arg(ap, double);
    D2 = (double) va_arg(ap, double);
    T1_old = (double) va_arg(ap, double);
    T2 = (double) va_arg(ap, double);
    Ts_old = (double) va_arg(ap, double);
    Told_node = (double *) va_arg(ap, double *);
    bubble = (double) va_arg(ap, double);
    dp = (double) va_arg(ap, double);
    expt = (double) va_arg(ap, double);
    ice0 = (double) va_arg(ap, double);
    kappa1 = (double) va_arg(ap, double);
    kappa2 = (double) va_arg(ap, double);
    max_moist = (double) va_arg(ap, double);
    moist = (double) va_arg(ap, double);

    root = (double  *) va_arg(ap, double  *);
    CanopLayerBnd = (double *) va_arg(ap, double *);

    /* meteorological forcing terms */
    UnderStory = (int) va_arg(ap, int);
    overstory = (int) va_arg(ap, int);

    NetShortBare = (double) va_arg(ap, double);
    NetShortGrnd = (double) va_arg(ap, double);
    NetShortSnow = (double) va_arg(ap, double);
    Tair = (double) va_arg(ap, double);
    atmos_density = (double) va_arg(ap, double);
    atmos_pressure = (double) va_arg(ap, double);
    emissivity = (double) va_arg(ap, double);
    LongBareIn = (double) va_arg(ap, double);
    LongSnowIn = (double) va_arg(ap, double);
    surf_atten = (double) va_arg(ap, double);
    vp = (double) va_arg(ap, double);
    vpd = (double) va_arg(ap, double);
    shortwave = (double) va_arg(ap, double);
    Catm = (double) va_arg(ap, double);
    dryFrac = (double *) va_arg(ap, double *);

    Wdew = (double *) va_arg(ap, double *);
    displacement = (double *) va_arg(ap, double *);
    ra = (double *) va_arg(ap, double *);
    Ra_veg = (double *) va_arg(ap, double *);
    Ra_used = (double *) va_arg(ap, double *);
    rainfall = (double) va_arg(ap, double);
    ref_height = (double *) va_arg(ap, double *);
    roughness = (double *) va_arg(ap, double *);
    wind = (double *) va_arg(ap, double *);

    /* latent heat terms */
    Le = (double) va_arg(ap, double);

    /* snowpack terms */
    Advection = (double) va_arg(ap, double);
    OldTSurf = (double) va_arg(ap, double);
    Tsnow_surf = (double) va_arg(ap, double);
    kappa_snow = (double) va_arg(ap, double);
    melt_energy = (double) va_arg(ap, double);
    snow_coverage = (double) va_arg(ap, double);
    snow_density = (double) va_arg(ap, double);
    snow_swq = (double) va_arg(ap, double);
    snow_water = (double) va_arg(ap, double);

    deltaCC = (double *) va_arg(ap, double *);
    refreeze_energy = (double *) va_arg(ap, double *);
    vapor_flux = (double *) va_arg(ap, double *);
    blowing_flux = (double *) va_arg(ap, double *);
    surface_flux = (double *) va_arg(ap, double *);

    /* soil node terms */
    Nnodes = (int) va_arg(ap, int);

    Cs_node = (double *) va_arg(ap, double *);
    T_node = (double *) va_arg(ap, double *);
    Tnew_node = (double *) va_arg(ap, double *);
    Tnew_fbflag = (char *) va_arg(ap, char *);
    Tnew_fbcount = (unsigned *) va_arg(ap, unsigned *);
    alpha = (double *) va_arg(ap, double *);
    beta = (double *) va_arg(ap, double *);
    bubble_node = (double *) va_arg(ap, double *);
    Zsum_node = (double *) va_arg(ap, double *);
    expt_node = (double *) va_arg(ap, double *);
    gamma = (double *) va_arg(ap, double *);
    ice_node = (double *) va_arg(ap, double *);
    kappa_node = (double *) va_arg(ap, double *);
    max_moist_node = (double *) va_arg(ap, double *);
    moist_node = (double *) va_arg(ap, double *);

    /* model structures */
    soil_con = (soil_con_struct *) va_arg(ap, soil_con_struct *);
    layer = (layer_data_struct *) va_arg(ap, layer_data_struct *);
    veg_var = (veg_var_struct *) va_arg(ap, veg_var_struct *);

    /* control flags */
    INCLUDE_SNOW = (int) va_arg(ap, int);
    NOFLUX = (int) va_arg(ap, int);
    EXP_TRANS = (int) va_arg(ap, int);
    SNOWING = (int) va_arg(ap, int);

    FIRST_SOLN = (int *) va_arg(ap, int *);

    /* returned energy balance terms */
    NetLongBare = (double *) va_arg(ap, double *);
    NetLongSnow = (double *) va_arg(ap, double *);
    T1 = (double *) va_arg(ap, double *);
    deltaH = (double *) va_arg(ap, double *);
    fusion = (double *) va_arg(ap, double *);
    grnd_flux = (double *) va_arg(ap, double *);
    latent_heat = (double *) va_arg(ap, double *);
    latent_heat_sub = (double *) va_arg(ap, double *);
    sensible_heat = (double *) va_arg(ap, double *);
    snow_flux = (double *) va_arg(ap, double *);
    store_error = (double *) va_arg(ap, double *);

    /* take additional variables from soil_con structure */
    b_infilt = soil_con->b_infilt;
    Wmax = soil_con->max_moist;
    Wcr = soil_con->Wcr;
    Wpwp = soil_con->Wpwp;
    depth = soil_con->depth;
    resid_moist = soil_con->resid_moist;
    elevation = (double)soil_con->elevation;
    frost_fract = soil_con->frost_fract;
    FS_ACTIVE = soil_con->FS_ACTIVE;
    /* more soil layer terms for IMPLICIT option*/
    bulk_dens_min = soil_con->bulk_dens_min;
    soil_dens_min = soil_con->soil_dens_min;
    quartz = soil_con->quartz;
    bulk_density = soil_con->bulk_density;
    soil_density = soil_con->soil_density;
    organic = soil_con->organic;


    /***************
       MAIN ROUTINE
    ***************/

    Error = 0;

    TMean = Ts;

    transp = calloc(options.Nlayer, sizeof(*transp));
    for (i = 0; i < options.Nlayer; i++) {
        transp[i] = 0.;
    }

    /**********************************************
       Compute Surface Temperature at Half Time Step
    **********************************************/

    if (snow_coverage > 0 && !INCLUDE_SNOW) {
        /****************************************
           Compute energy flux through snow pack
        ****************************************/

        *snow_flux = (kappa_snow * (Tsnow_surf - TMean));
    }
    else if (INCLUDE_SNOW) {
        *snow_flux = 0.;
        Tsnow_surf = TMean;
    }
    else {
        *snow_flux = 0.;
    }

    /***************************************************************
       Estimate soil temperatures for ground heat flux calculations
    ***************************************************************/

    if (options.QUICK_FLUX) {
        /**************************************************************
           Use Liang et al. 1999 Equations to Calculate Ground Heat Flux
           NOTE: T2 is not the temperature of layer 2, nor of node 2, nor at depth dp;
           T2 is the constant temperature at depths much greater than dp to which
           the soil temperature profile is asymptotic.
        **************************************************************/
        *T1 = estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs2,
                          dp, delta_t);

        /*****************************************************
           Compute the Ground Heat Flux from the Top Soil Layer
        *****************************************************/
        if (options.GRND_FLUX_TYPE == GF_406) {
            *grnd_flux =
                (snow_coverage +
                 (1. -
                  snow_coverage) *
                 surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
        }
        else {
            *grnd_flux =
                (snow_coverage +
                 (1. -
                  snow_coverage) *
                 surf_atten) * (kappa1 / D1 * ((*T1) - TMean) +
                                (
                                    kappa2 / D2 *
                                    (1. -
                                     exp(-D1 /
                                         dp)) * (T2 - (*T1)))) / 2.;
        }
    }
    else {
        /*************************************************************
           Use Finite Difference Method to Solve Ground Heat
           Flux at Soil Thermal Nodes (Cherkauer and Lettenmaier, 1999)
        *************************************************************/
        T_node[0] = TMean;

        /* IMPLICIT Solution */
        if (options.IMPLICIT) {
            Error = solve_T_profile_implicit(Tnew_node, T_node, Tnew_fbflag,
                                             Tnew_fbcount, Zsum_node,
                                             kappa_node, Cs_node, moist_node,
                                             delta_t, max_moist_node,
                                             bubble_node, expt_node, ice_node,
                                             alpha, beta, gamma, dp, Nnodes,
                                             FIRST_SOLN, NOFLUX, EXP_TRANS,
                                             bulk_dens_min, soil_dens_min,
                                             quartz, bulk_density,
                                             soil_density, organic, depth);

            if (FIRST_SOLN[1]) {
                FIRST_SOLN[1] = false;
            }
        }

        /* EXPLICIT Solution, or if IMPLICIT Solution Failed */
        if (!options.IMPLICIT || Error == 1) {
            if (options.IMPLICIT) {
                FIRST_SOLN[0] = true;
            }
            Error = solve_T_profile(Tnew_node, T_node, Tnew_fbflag,
                                    Tnew_fbcount, Zsum_node, kappa_node,
                                    Cs_node, moist_node, delta_t,
                                    max_moist_node, bubble_node,
                                    expt_node, ice_node, alpha, beta, gamma, dp,
                                    Nnodes, FIRST_SOLN, FS_ACTIVE, NOFLUX,
                                    EXP_TRANS);
        }

        if ((int) Error == ERROR) {
            log_err("func_surf_energy_bal calling solve_T_profile");
        }
        /* Compute temperatures for calculations of ground heat flux, delta_H, and fusion */
        if (!options.EXP_TRANS) {
            *T1 = Tnew_node[1];
            T11 = Tnew_node[2];
        }
        else {
            i = 0;
            while (soil_con->Zsum_node[i] < D1) {
                i++;
            }
            D1_minus = soil_con->Zsum_node[i - 1];
            D1_plus = soil_con->Zsum_node[i];
            T1_minus = Tnew_node[i - 1];
            T1_plus = Tnew_node[i];
            *T1 = T1_minus +
                  (T1_plus - T1_minus) * (D1 - D1_minus) / (D1_plus - D1_minus);
        }


        /*****************************************************
           Compute the Ground Heat Flux between layers 0 and 1
        *****************************************************/
        if (!options.EXP_TRANS) {
            if (options.GRND_FLUX_TYPE == GF_406) {
                *grnd_flux =
                    (snow_coverage +
                     (1. -
                      snow_coverage) *
                     surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
            }
            else {
                *grnd_flux =
                    (snow_coverage +
                     (1. -
                      snow_coverage) *
                     surf_atten) * (kappa1 / D1 * ((*T1) - TMean) +
                                    (
                                        kappa2 / D2 *
                                        (T11 - (*T1)))) / 2.;
            }
        }
        else {
            *grnd_flux =
                (snow_coverage +
                 (1. -
                  snow_coverage) *
                 surf_atten) * (kappa1 / (D1 - D1_minus) * ((*T1) - T1_minus) +
                                (
                                    kappa2 /
                                    (D1_plus - D1) * (T1_plus - (*T1)))) / 2.;
        }
    }

    /******************************************************
       Compute the change in heat storage in the region between layers 0 and 1
    ******************************************************/
    if (!options.EXP_TRANS) {
        *deltaH =
            (Cs1 * ((Ts_old + T1_old) - (TMean + *T1)) * D1 / delta_t / 2.);
    }
    else {
        *deltaH = 0;
        i = 0;
        while (soil_con->Zsum_node[i + 1] < D1) {
            *deltaH +=
                (Cs1 *
                 ((Told_node[i] +
                   Told_node[i +
                             1]) -
                  (Tnew_node[i] +
                   Tnew_node[i +
                             1])) *
                 (soil_con->Zsum_node[i +
                                      1] -
                  soil_con->Zsum_node[i]) / delta_t / 2.);
            i++;
        }
        *deltaH +=
            (Cs1 *
             ((Told_node[i] +
               T1_old) -
              (Tnew_node[i] +
               *T1)) * (D1 - soil_con->Zsum_node[i]) / delta_t / 2.);
    }

    /******************************************************
       Compute the change in heat due to solid-liquid phase changes in the region between layers 0 and 1
    ******************************************************/
    if (FS_ACTIVE && options.FROZEN_SOIL) {
        if (!options.EXP_TRANS) {
            if ((TMean + *T1) / 2. < 0.) {
                ice = moist - maximum_unfrozen_water((TMean + *T1) / 2.,
                                                     max_moist, bubble, expt);
                if (ice < 0.) {
                    ice = 0.;
                }
            }
            else {
                ice = 0.;
            }
            *fusion =
                (-CONST_RHOICE * CONST_LATICE * (ice0 - ice) * D1 / delta_t);
        }
        else {
            *fusion = 0;
            i = 0;
            while (soil_con->Zsum_node[i + 1] < D1) {
                if ((Told_node[i] + Told_node[i + 1]) / 2. < 0.) {
                    ice0 = moist -
                           maximum_unfrozen_water(
                        (Told_node[i] + Told_node[i + 1]) / 2.,
                        max_moist, bubble,
                        expt);
                    if (ice0 < 0.) {
                        ice0 = 0.;
                    }
                }
                else {
                    ice0 = 0.;
                }
                if ((Tnew_node[i] + Tnew_node[i + 1]) / 2. < 0.) {
                    ice = moist -
                          maximum_unfrozen_water(
                        (Tnew_node[i] + Tnew_node[i + 1]) / 2.,
                        max_moist, bubble,
                        expt);
                    if (ice < 0.) {
                        ice = 0.;
                    }
                }
                else {
                    ice = 0.;
                }
                *fusion +=
                    (-CONST_RHOICE * CONST_LATICE *
                     (ice0 -
                      ice) *
                     (soil_con->Zsum_node[i +
                                          1] -
                      soil_con->Zsum_node[i]) / delta_t);
                i++;
            }
            if ((Told_node[i] + T1_old) / 2. < 0.) {
                ice0 = moist - maximum_unfrozen_water(
                    (Told_node[i] + T1_old) / 2.,
                    max_moist, bubble, expt);
                if (ice0 < 0.) {
                    ice0 = 0.;
                }
            }
            else {
                ice0 = 0.;
            }
            if ((Tnew_node[i] + *T1) / 2. < 0.) {
                ice = moist - maximum_unfrozen_water((Tnew_node[i] + *T1) / 2.,
                                                     max_moist, bubble, expt);
                if (ice < 0.) {
                    ice = 0.;
                }
            }
            else {
                ice = 0.;
            }
            *fusion +=
                (-CONST_RHOICE * CONST_LATICE *
                 (ice0 - ice) * (D1 - soil_con->Zsum_node[i]) / delta_t);
        }
    }

    /* if thin snowpack, compute the change in energy stored in the pack */
    if (INCLUDE_SNOW) {
        if (TMean > 0.) {
            *deltaCC = CONST_CPICE *
                       (snow_swq - snow_water) * (0. - OldTSurf) / delta_t;
        }
        else {
            *deltaCC = CONST_CPICE *
                       (snow_swq - snow_water) * (TMean - OldTSurf) / delta_t;
        }
        *refreeze_energy = (snow_water * CONST_LATICE * snow_density) / delta_t;
        *deltaCC *= snow_coverage; // adjust for snow cover fraction
        *refreeze_energy *= snow_coverage; // adjust for snow cover fraction
    }

    /** Compute net surface radiation of snow-free area for evaporation estimates **/
    LongBareOut = calc_outgoing_longwave(TMean + CONST_TKFRZ, param.EMISS_GRND);
    if (INCLUDE_SNOW) { // compute net LW at snow surface
        (*NetLongSnow) = (LongSnowIn - snow_coverage * LongBareOut);
    }
    (*NetLongBare) = (LongBareIn - (1. - snow_coverage) * LongBareOut); // net LW snow-free area
    NetBareRad =
        (NetShortBare + (*NetLongBare) + *grnd_flux + *deltaH + *fusion);

    /** Compute atmospheric stability correction **/
    if (wind[UnderStory] > 0.0 && overstory && SNOWING) {
        Ra_veg[0] = ra[UnderStory] /
                    StabilityCorrection(ref_height[UnderStory], 0., TMean,
                                        Tair,
                                        wind[UnderStory],
                                        roughness[UnderStory]);
    }
    else if (wind[UnderStory] > 0.0) {
        Ra_veg[0] = ra[UnderStory] /
                    StabilityCorrection(ref_height[UnderStory],
                                        displacement[UnderStory],
                                        TMean, Tair, wind[UnderStory],
                                        roughness[UnderStory]);
    }
    else {
        Ra_veg[0] = param.HUGE_RESIST;
    }
    Ra_used[0] = Ra_veg[0];

    /*************************************************
       Compute aerodynamic resistance for the case of exposed soil between
       plants (or gaps in canopy).  Assume plants and exposed soil are well-
       mixed, i.e., exposed soil is neither as disconnected from the
       atmosphere as soil under veg or canopy, nor as exposed as a large
       area of exposed soil.  Rather, it is subject to some wind attenuation
       from the surrounding plants.  Thus, compute as the area-weighted
       average of "pure" understory and "pure" exposed conditions.

       NOTE: can't average the resistances; must convert to conductances,
       average the conductances, and then convert the average conductance
       to a resistance.

       NOTE 2: since a resistance of 0 corresponds to an infinite conductance,
       if either Ra_veg (resistance under veg) or Ra_bare (resistance over
       exposed soil) are 0, then Ra_used must necessarily be 0 as well.
    *************************************************/
    if (veg_var->fcanopy < 1) {
        /** If Ra_veg is non-zero, use it to compute area-weighted average **/
        if (Ra_veg[0] > 0) {
            /** aerodynamic conductance under vegetation **/
            ga_veg = 1 / Ra_veg[0];
            /** compute aerodynamic resistance over exposed soil (Ra_bare) **/
            tmp_wind[0] = wind[0];
            tmp_wind[1] = MISSING; // unused
            tmp_wind[2] = MISSING; // unused
            tmp_height = soil_con->rough / param.VEG_RATIO_RL_HEIGHT;
            tmp_displacement[0] = calc_veg_displacement(tmp_height);
            tmp_roughness[0] = soil_con->rough;
            tmp_ref_height[0] = param.SOIL_WINDH; // wind height over bare soil
            Error = CalcAerodynamic(0, 0, 0, soil_con->snow_rough,
                                    soil_con->rough, 0, Ra_bare, tmp_wind,
                                    tmp_displacement, tmp_ref_height,
                                    tmp_roughness);
            Ra_bare[0] /= StabilityCorrection(tmp_ref_height[0],
                                              tmp_displacement[0], TMean,
                                              Tair, tmp_wind[0],
                                              tmp_roughness[0]);

            /** if Ra_bare is non-zero, compute area-weighted average
                aerodynamic conductance **/
            if (Ra_bare[0] > 0) {
                /** aerodynamic conductance over exposed soil **/
                ga_bare = 1 / Ra_bare[0];
                /** area-weighted average aerodynamic conductance **/
                ga_average = veg_var->fcanopy * ga_veg +
                             (1 - veg_var->fcanopy) * ga_bare;
                /** aerodynamic resistance is inverse of conductance **/
                Ra_used[0] = 1 / ga_average;
            }
            /** else aerodynamic resistance is zero **/
            else {
                Ra_used[0] = 0;
            }
        }
        /** else aerodynamic resistance is zero **/
        else {
            Ra_used[0] = 0;
        }
    }

    /*************************************************
       Compute Evapotranspiration if not snow covered

       Should evapotranspiration be active when the
       ground is only partially covered with snow????

       Use Arno Evap in the exposed soil portion, and/or
       if LAI is zero.
    *************************************************/
    if (VEG && !SNOWING && veg_var->fcanopy > 0) {
        Evap = canopy_evap(layer, veg_var, true,
                           veg_class, Wdew, delta_t, NetBareRad, vpd,
                           NetShortBare, Tair, Ra_veg[1], elevation, rainfall,
                           Wmax, Wcr, Wpwp, frost_fract, root, dryFrac,
                           shortwave, Catm, CanopLayerBnd);
        if (veg_var->fcanopy < 1) {
            for (i = 0; i < options.Nlayer; i++) {
                transp[i] = layer[i].evap;
                layer[i].evap = 0.;
            }
            Evap *= veg_var->fcanopy;
            Evap += (1 - veg_var->fcanopy) *
                    arno_evap(layer, surf_atten * NetBareRad, Tair, vpd,
                              depth[0], max_moist * depth[0] * MM_PER_M,
                              elevation, b_infilt, Ra_used[0], delta_t,
                              resid_moist[0], frost_fract);
            for (i = 0; i < options.Nlayer; i++) {
                layer[i].evap = veg_var->fcanopy * transp[i] +
                                (1 - veg_var->fcanopy) * layer[i].evap;
                if (layer[i].evap > 0.) {
                    layer[i].bare_evap_frac = 1 -
                                              (veg_var->fcanopy *
                                               transp[i]) / layer[i].evap;
                }
                else {
                    layer[i].bare_evap_frac = 0.;
                }
            }
            veg_var->throughfall =
                (1 -
                 veg_var->fcanopy) * rainfall + veg_var->fcanopy *
                veg_var->throughfall;
            veg_var->canopyevap *= veg_var->fcanopy;
            veg_var->Wdew *= veg_var->fcanopy;
        }
        else {
            for (i = 0; i < options.Nlayer; i++) {
                layer[i].bare_evap_frac = 0.;
            }
        }
    }
    else if (!SNOWING) {
        Evap = arno_evap(layer, NetBareRad, Tair, vpd,
                         depth[0], max_moist * depth[0] * MM_PER_M,
                         elevation, b_infilt, Ra_used[0], delta_t,
                         resid_moist[0], frost_fract);
        for (i = 0; i < options.Nlayer; i++) {
            layer[i].bare_evap_frac = 1;
        }
    }
    else {
        Evap = 0.;
    }

    free(transp);

    /**********************************************************************
       Compute the Latent Heat Flux from the Surface and Covering Vegetation
    **********************************************************************/
    *latent_heat = -CONST_RHOFW * Le * Evap;
    *latent_heat_sub = 0.;

    /** Compute the latent heat flux from a thin snowpack if present **/
    if (INCLUDE_SNOW) {
        /* Convert sublimation terms from m/timestep to kg/m2s */
        VaporMassFlux = *vapor_flux * CONST_RHOICE / delta_t;
        BlowingMassFlux = *blowing_flux * CONST_RHOICE / delta_t;
        SurfaceMassFlux = *surface_flux * CONST_RHOICE / delta_t;

        latent_heat_from_snow(atmos_density, vp, Le, atmos_pressure,
                              Ra_used[0], TMean, vpd, &temp_latent_heat,
                              &temp_latent_heat_sub, &VaporMassFlux,
                              &BlowingMassFlux, &SurfaceMassFlux);
        *latent_heat += temp_latent_heat * snow_coverage;
        *latent_heat_sub = temp_latent_heat_sub * snow_coverage;

        /* Convert sublimation terms from kg/m2s to m/timestep */
        *vapor_flux = VaporMassFlux * delta_t / CONST_RHOICE;
        *blowing_flux = BlowingMassFlux * delta_t / CONST_RHOICE;
        *surface_flux = SurfaceMassFlux * delta_t / CONST_RHOICE;
    }
    else {
        *latent_heat *= (1. - snow_coverage);
    }

    /************************************************
       Compute the Sensible Heat Flux from the Surface
    ************************************************/
    if (snow_coverage < 1 || INCLUDE_SNOW) {
        *sensible_heat =
            calc_sensible_heat(atmos_density, Tair, TMean, Ra_used[0]);
        if (!INCLUDE_SNOW) {
            (*sensible_heat) *= (1. - snow_coverage);
        }
    }
    else {
        *sensible_heat = 0.;
    }

    /*************************************
       Compute Surface Energy Balance Error
    *************************************/
    error = (NetBareRad + // net radiation on snow-free area
             NetShortGrnd + NetShortSnow + // net surface SW
             emissivity * (*NetLongSnow)) + // net surface LW
            *sensible_heat + // surface sensible heat
            (*latent_heat + *latent_heat_sub)  // surface latent heats
            /* heat flux through snowpack - for snow covered fraction */
            + *snow_flux * snow_coverage
            /* energy used in reducing snow coverage area */
            + melt_energy
            /* snow energy terms - values are 0 unless INCLUDE_SNOW */
            + Advection - *deltaCC;

    if (INCLUDE_SNOW) {
        if (Tsnow_surf == 0.0 && error > -(*refreeze_energy)) {
            *refreeze_energy = -error;
            error = 0.0;
        }
        else {
            error += *refreeze_energy;
        }
    }

    *store_error = error;

    return error;
}
