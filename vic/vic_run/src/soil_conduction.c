/******************************************************************************
* @section DESCRIPTION
*
* Calculate soil thermal conduction.
*
* @section LICENSE
*
* The Variable Infiltration Capacity (VIC) macroscale hydrological model
* Copyright (C) 2014  The Land Surface Hydrology Group, Department of Civil
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
******************************************************************************/

#include <vic_run.h>

/******************************************************************************
* @brief    Soil thermal conductivity calculated using Johansen's method.
*
* @note     Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
*               Chapter 7: Methods for Calculating the Thermal Conductivity
*               of Soils
******************************************************************************/
double
soil_conductivity(double moist,
                  double Wu,
                  double soil_dens_min,
                  double bulk_dens_min,
                  double quartz,
                  double soil_density,
                  double bulk_density,
                  double organic)
{
    double Ke;
    double Ki = 2.2;    /* thermal conductivity of ice (W/mK) */
    double Kw = 0.57;   /* thermal conductivity of water (W/mK) */
    double Ksat;
    double Kdry;        /* Dry thermal conductivity of soil (W/mK), including mineral and organic fractions */
    double Kdry_org = 0.05; /* Dry thermal conductivity of organic fraction (W/mK) (Farouki 1981) */
    double Kdry_min;    /* Dry thermal conductivity of mineral fraction (W/mK) */
    double Ks;          /* thermal conductivity of solid (W/mK), including mineral and organic fractions */
    double Ks_org = 0.25; /* thermal conductivity of organic fraction of solid (W/mK) (Farouki 1981) */
    double Ks_min;      /* thermal conductivity of mineral fraction of solid (W/mK) */
    double Sr;          /* fractional degree of saturation */
    double K;
    double porosity;

    /* Calculate dry conductivity as weighted average of mineral and organic fractions. */
    Kdry_min =
        (0.135 * bulk_dens_min +
         64.7) / (soil_dens_min - 0.947 * bulk_dens_min);
    Kdry = (1 - organic) * Kdry_min + organic * Kdry_org;

    if (moist > 0.) {
        porosity = 1.0 - bulk_density / soil_density; // NOTE: if excess_ice present,
                                                      // this is actually effective_porosity

        Sr = moist / porosity;

        // Compute Ks of mineral soil; here "quartz" is the fraction (quartz volume / mineral soil volume)
        if (quartz < .2) {
            Ks_min = pow(7.7, quartz) * pow(3.0, 1.0 - quartz); // when quartz is less than 0.2
        }
        else {
            Ks_min = pow(7.7, quartz) * pow(2.2, 1.0 - quartz); // when quartz is greater than 0.2
        }
        Ks = (1 - organic) * Ks_min + organic * Ks_org;

        if (Wu == moist) {
            /** Soil unfrozen **/
            Ksat = pow(Ks, 1.0 - porosity) * pow(Kw, porosity);
            Ke = 0.7 * log10(Sr) + 1.0;
        }
        else {
            /** Soil frozen **/
            Ksat =
                pow(Ks, 1.0 - porosity) * pow(Ki, porosity - Wu) * pow(Kw, Wu);
            Ke = Sr;
        }

        K = (Ksat - Kdry) * Ke + Kdry;
        if (K < Kdry) {
            K = Kdry;
        }
    }
    else {
        K = Kdry;
    }

    return (K);
}

/******************************************************************************
* @brief    This subroutine calculates the soil volumetric heat capacity
            based on the fractional volume of its component parts.
******************************************************************************/
double
volumetric_heat_capacity(double soil_fract,
                         double water_fract,
                         double ice_fract,
                         double organic_fract)
{
    double Cs;

    // Constant values are volumetric heat capacities in J/m^3/K
    Cs = 2.0e6 * soil_fract * (1 - organic_fract);
    Cs += 2.7e6 * soil_fract * organic_fract;
    Cs += 4.2e6 * water_fract;
    Cs += 1.9e6 * ice_fract;
    Cs += 1.3e3 * (1. - (soil_fract + water_fract + ice_fract)); // air

    return (Cs);
}

/******************************************************************************
* @brief    This subroutine sets the thermal node soil parameters to constant
*           values based on those defined for the current grid cells soil type.
*           Thermal node propertiers for the energy balance solution are also
*           set (these constants are used to reduce the solution time required
*           within each iteration).
******************************************************************************/
void
set_node_parameters(double *Zsum_node,
                    double *max_moist_node,
                    double *expt_node,
                    double *bubble_node,
                    double *alpha,
                    double *beta,
                    double *gamma,
                    double *depth,
                    double *max_moist,
                    double *expt,
                    double *bubble,
                    int     Nnodes,
                    int     Nlayers)
{
    extern option_struct options;

    char                 PAST_BOTTOM;
    int                  nidx, lidx;
    double               Lsum; /* cumulative depth of moisture layer */

    PAST_BOTTOM = false;
    lidx = 0;
    Lsum = 0.;

    /* set node parameters */
    for (nidx = 0; nidx < Nnodes; nidx++) {
        if (Zsum_node[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx !=
            Nlayers - 1) {
            /* node on layer boundary */
            max_moist_node[nidx] = (max_moist[lidx] / depth[lidx] +
                                    max_moist[lidx +
                                              1] /
                                    depth[lidx + 1]) / MM_PER_M / 2.;
            expt_node[nidx] = (expt[lidx] + expt[lidx + 1]) / 2.;
            bubble_node[nidx] = (bubble[lidx] + bubble[lidx + 1]) / 2.;
        }
        else {
            /* node completely in layer */
            max_moist_node[nidx] = max_moist[lidx] / depth[lidx] / MM_PER_M;
            expt_node[nidx] = expt[lidx];
            bubble_node[nidx] = bubble[lidx];
        }
        if (Zsum_node[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
            Lsum += depth[lidx];
            lidx++;
            if (lidx == Nlayers) {
                PAST_BOTTOM = true;
                lidx = Nlayers - 1;
            }
        }
    }

    /* set constant variables for thermal calculations */
    for (nidx = 0; nidx < Nnodes - 2; nidx++) {
        alpha[nidx] = Zsum_node[nidx + 2] - Zsum_node[nidx];
        beta[nidx] = Zsum_node[nidx + 1] - Zsum_node[nidx];
        gamma[nidx] = Zsum_node[nidx + 2] - Zsum_node[nidx + 1];
    }
    if (options.NOFLUX) {
        /* no flux bottom boundary activated */
        alpha[Nnodes -
              2] = 2. * (Zsum_node[Nnodes - 1] - Zsum_node[Nnodes - 2]);
        beta[nidx] = Zsum_node[Nnodes - 1] - Zsum_node[Nnodes - 2];
        gamma[nidx] = Zsum_node[Nnodes - 1] - Zsum_node[Nnodes - 2];
    }
}

/******************************************************************************
* @brief    This subroutine determines the moisture and ice contents of each
*           soil thermal node based on the current node temperature and layer
*           moisture content.  Thermal conductivity and volumetric heat
*           capacity are then estimated for each node based on the division of
*           moisture contents.
******************************************************************************/
int
distribute_node_moisture_properties(double *moist_node,
                                    double *ice_node,
                                    double *kappa_node,
                                    double *Cs_node,
                                    double *Zsum_node,
                                    double *T_node,
                                    double *max_moist_node,
                                    double *expt_node,
                                    double *bubble_node,
                                    double *moist,
                                    double *depth,
                                    double *soil_dens_min,
                                    double *bulk_dens_min,
                                    double *quartz,
                                    double *soil_density,
                                    double *bulk_density,
                                    double *organic,
                                    int     Nnodes,
                                    int     Nlayers,
                                    char    FS_ACTIVE)
{
    extern option_struct     options;
    extern parameters_struct param;

    char                     PAST_BOTTOM;
    int                      nidx, lidx;
    double                   Lsum; /* cumulative depth of moisture layer */

    lidx = 0;
    Lsum = 0.;
    PAST_BOTTOM = false;

    /* node estimates */
    for (nidx = 0; nidx < Nnodes; nidx++) {
        if (!PAST_BOTTOM || param.SOIL_SLAB_MOIST_FRACT < 0) {
            if (Zsum_node[nidx] == Lsum + depth[lidx] && nidx != 0 && lidx !=
                Nlayers - 1) {
                /* node on layer boundary */
                moist_node[nidx] = (moist[lidx] / depth[lidx] +
                                    moist[lidx +
                                          1] / depth[lidx + 1]) / MM_PER_M / 2.;
            }
            else {
                /* node completely in layer */
                moist_node[nidx] = moist[lidx] / depth[lidx] / MM_PER_M;
            }
        }
        else {
            // use constant soil moisture below bottom soil layer
            moist_node[nidx] = param.SOIL_SLAB_MOIST_FRACT *
                               max_moist_node[nidx];
        }


        // Check that node moisture does not exceed maximum node moisture
        if (moist_node[nidx] - max_moist_node[nidx] > 0) {
            // HACK!!!!!!!!!!!
            moist_node[nidx] = max_moist_node[nidx];
        }
        if (T_node[nidx] < 0 && (FS_ACTIVE && options.FROZEN_SOIL)) {
            /* compute moisture and ice contents */
            ice_node[nidx] =
                moist_node[nidx] - maximum_unfrozen_water(T_node[nidx],
                                                          max_moist_node[nidx],
                                                          bubble_node[nidx],
                                                          expt_node[nidx]);
            if (ice_node[nidx] < 0) {
                ice_node[nidx] = 0;
            }

            /* compute thermal conductivity */
            kappa_node[nidx] =
                soil_conductivity(moist_node[nidx], moist_node[nidx] -
                                  ice_node[nidx],
                                  soil_dens_min[lidx], bulk_dens_min[lidx],
                                  quartz[lidx],
                                  soil_density[lidx], bulk_density[lidx],
                                  organic[lidx]);
        }
        else {
            /* compute moisture and ice contents */
            ice_node[nidx] = 0;
            /* compute thermal conductivity */
            kappa_node[nidx] =
                soil_conductivity(moist_node[nidx], moist_node[nidx],
                                  soil_dens_min[lidx], bulk_dens_min[lidx],
                                  quartz[lidx],
                                  soil_density[lidx], bulk_density[lidx],
                                  organic[lidx]);
        }
        /* compute volumetric heat capacity */
        Cs_node[nidx] = volumetric_heat_capacity(
            bulk_density[lidx] / soil_density[lidx],
            moist_node[nidx] -
            ice_node[nidx], ice_node[nidx], organic[lidx]);

        if (Zsum_node[nidx] > Lsum + depth[lidx] && !PAST_BOTTOM) {
            Lsum += depth[lidx];
            lidx++;
            if (lidx == Nlayers) {
                PAST_BOTTOM = true;
                lidx = Nlayers - 1;
            }
        }
    }
    return (0);
}

/******************************************************************************
* @brief    This subroutine estimates the temperature and depth of all frost
*           sub-ares within each soil thermal node based on the distribution of
*           soil thermal node temperatures.
******************************************************************************/
void
estimate_frost_temperature_and_depth(double ***tmpT,
                                     double  **tmpZ,
                                     double   *Zsum_node,
                                     double   *T,
                                     double   *depth,
                                     double   *frost_fract,
                                     double    frost_slope,
                                     size_t    Nnodes,
                                     size_t    Nlayers)
{
    extern option_struct options;

    size_t               nidx, min_nidx;
    size_t               lidx, frost_area, max_nidx;
    double               Lsum[MAX_LAYERS + 1];
    double               min_temp, max_temp, tmp_fract;

    // compute cumulative layer depths
    Lsum[0] = 0;
    for (lidx = 1; lidx <= Nlayers; lidx++) {
        Lsum[lidx] = depth[lidx - 1] + Lsum[lidx - 1];
    }

    // estimate soil layer average variables
    for (lidx = 0; lidx < Nlayers; lidx++) {
        // Bracket current layer between nodes
        min_nidx = Nnodes - 2;
        while (Lsum[lidx] < Zsum_node[min_nidx] && min_nidx > 0) {
            min_nidx--;
        }
        max_nidx = 1;
        while (Lsum[lidx + 1] > Zsum_node[max_nidx] && max_nidx < Nnodes) {
            max_nidx++;
        }
        if (max_nidx >= Nnodes) {
            log_warn("Soil thermal nodes do not extend below bottom "
                     "soil layer; using deepest node temperature for "
                     "all deeper depths.");
            // If we get here, soil thermal nodes don't extend all the way
            // down to the bottom of the lowest layer.  In this case, just
            // use the deepest node to represent all deeper temperatures.
            max_nidx = Nnodes - 1;
        }

        // Get soil node temperatures for current layer
        if (Zsum_node[min_nidx] < Lsum[lidx]) {
            tmpT[lidx][min_nidx][options.Nfrost] = linear_interp(
                Lsum[lidx],
                Zsum_node[min_nidx],
                Zsum_node[min_nidx + 1],
                T[min_nidx],
                T[min_nidx + 1]);
        }
        else {
            tmpT[lidx][min_nidx][options.Nfrost] = T[min_nidx];
        }
        tmpZ[lidx][min_nidx] = Lsum[lidx];
        for (nidx = min_nidx + 1; nidx < max_nidx; nidx++) {
            tmpT[lidx][nidx][options.Nfrost] = T[nidx];
            tmpZ[lidx][nidx] = Zsum_node[nidx];
        }
        if (Zsum_node[max_nidx] > Lsum[lidx + 1]) {
            tmpT[lidx][max_nidx][options.Nfrost] = linear_interp(
                Lsum[lidx + 1],
                Zsum_node[max_nidx - 1],
                Zsum_node[max_nidx],
                T[max_nidx - 1],
                T[max_nidx]);
        }
        else {
            tmpT[lidx][max_nidx][options.Nfrost] = T[max_nidx];
        }
        tmpZ[lidx][max_nidx] = Lsum[lidx + 1];

        // distribute temperatures for sub-areas
        for (nidx = min_nidx; nidx <= max_nidx; nidx++) {
            min_temp = tmpT[lidx][nidx][options.Nfrost] - frost_slope / 2.;
            max_temp = min_temp + frost_slope;
            for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
                if (options.Nfrost > 1) {
                    if (frost_area == 0) {
                        tmp_fract = frost_fract[0] / 2.;
                    }
                    else {
                        tmp_fract += (frost_fract[frost_area - 1] / 2. +
                                      frost_fract[frost_area] / 2.);
                    }
                    tmpT[lidx][nidx][frost_area] = linear_interp(tmp_fract, 0,
                                                                 1,
                                                                 min_temp,
                                                                 max_temp);
                }
                else {
                    tmpT[lidx][nidx][frost_area] =
                        tmpT[lidx][nidx][options.Nfrost];
                }
            }
        }
    }
}

/******************************************************************************
* @brief    This subroutine estimates the temperature of all soil moisture
*           layers based on the distribution of soil thermal node temperatures.
******************************************************************************/
int
estimate_layer_temperature(layer_data_struct *layer,
                           double          ***tmpT,
                           double           **tmpZ,
                           double            *Zsum_node,
                           double            *depth,
                           size_t             Nnodes,
                           size_t             Nlayers)
{
    extern option_struct options;

    size_t               nidx, min_nidx;
    size_t               lidx, max_nidx;
    double               Lsum[MAX_LAYERS + 1];

    // compute cumulative layer depths
    Lsum[0] = 0;
    for (lidx = 1; lidx <= Nlayers; lidx++) {
        Lsum[lidx] = depth[lidx - 1] + Lsum[lidx - 1];
    }

    // estimate soil layer average temperature
    for (lidx = 0; lidx < Nlayers; lidx++) {
        // Initialize layer temperature
        layer[lidx].T = 0.;

        // Bracket current layer between nodes
        min_nidx = Nnodes - 2;
        while (Lsum[lidx] < Zsum_node[min_nidx] && min_nidx > 0) {
            min_nidx--;
        }
        max_nidx = 1;
        while (Lsum[lidx + 1] > Zsum_node[max_nidx] && max_nidx < Nnodes) {
            max_nidx++;
        }
        if (max_nidx >= Nnodes) {
            log_warn("Soil thermal nodes do not extend below bottom "
                     "soil layer; using deepest node temperature for "
                     "all deeper depths.");
            // If we get here, soil thermal nodes don't extend all the way
            // down to the bottom of the lowest layer.  In this case, just
            // use the deepest node to represent all deeper temperatures.
            max_nidx = Nnodes - 1;
        }

        // Compute average soil layer temperature
        for (nidx = min_nidx; nidx < max_nidx; nidx++) {
            layer[lidx].T +=
                (tmpZ[lidx][nidx + 1] - tmpZ[lidx][nidx]) *
                (tmpT[lidx][nidx + 1][options.Nfrost] +
                 tmpT[lidx][nidx][options.Nfrost]) / 2.;
        }
        layer[lidx].T /= depth[lidx];
    }

    return (0);
}

/******************************************************************************
* @brief    This subroutine estimates the ice content of all soil moisture
*           layers based on the distribution of soil thermal node temperatures.
******************************************************************************/
int
estimate_layer_ice_content(layer_data_struct *layer,
                           double          ***tmpT,
                           double           **tmpZ,
                           double            *Zsum_node,
                           double            *depth,
                           double            *max_moist,
                           double            *expt,
                           double            *bubble,
                           size_t             Nnodes,
                           size_t             Nlayers,
                           char               FS_ACTIVE)
{
    extern option_struct options;

    size_t               nidx, min_nidx;
    size_t               lidx, frost_area, max_nidx;
    double               Lsum[MAX_LAYERS + 1];
    double               tmp_ice[MAX_NODES][MAX_FROST_AREAS];

    // compute cumulative layer depths
    Lsum[0] = 0;
    for (lidx = 1; lidx <= Nlayers; lidx++) {
        Lsum[lidx] = depth[lidx - 1] + Lsum[lidx - 1];
    }

    // estimate soil layer average ice content
    for (lidx = 0; lidx < Nlayers; lidx++) {
        // Initialize layer ice content
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            layer[lidx].ice[frost_area] = 0.;
        }

        // Bracket current layer between nodes
        min_nidx = Nnodes - 2;
        while (Lsum[lidx] < Zsum_node[min_nidx] && min_nidx > 0) {
            min_nidx--;
        }
        max_nidx = 1;
        while (Lsum[lidx + 1] > Zsum_node[max_nidx] && max_nidx < Nnodes) {
            max_nidx++;
        }
        if (max_nidx >= Nnodes) {
            log_warn("Soil thermal nodes do not extend below bottom "
                     "soil layer; using deepest node temperature for "
                     "all deeper depths.");
            // If we get here, soil thermal nodes don't extend all the way
            // down to the bottom of the lowest layer.  In this case, just
            // use the deepest node to represent all deeper temperatures.
            max_nidx = Nnodes - 1;
        }

        // Get soil node ice content for current layer
        if (options.FROZEN_SOIL && FS_ACTIVE) {
            for (nidx = min_nidx; nidx <= max_nidx; nidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    tmp_ice[nidx][frost_area] = layer[lidx].moist -
                                                maximum_unfrozen_water(
                        tmpT[lidx][nidx][frost_area], max_moist[lidx],
                        bubble[lidx],
                        expt[lidx]);
                    if (tmp_ice[nidx][frost_area] < 0) {
                        tmp_ice[nidx][frost_area] = 0.;
                    }
                }
            }
        }
        else {
            for (nidx = min_nidx; nidx <= max_nidx; nidx++) {
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    tmp_ice[nidx][frost_area] = 0;
                }
            }
        }

        // Compute average soil layer ice content
        for (nidx = min_nidx; nidx < max_nidx; nidx++) {
            for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
                layer[lidx].ice[frost_area] +=
                    (tmpZ[lidx][nidx + 1] - tmpZ[lidx][nidx]) *
                    (tmp_ice[nidx + 1][frost_area] +
                     tmp_ice[nidx][frost_area]) / 2.;
            }
        }
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            layer[lidx].ice[frost_area] /= depth[lidx];
        }
    }

    return (0);
}

/******************************************************************************
* @brief    This subroutine estimates the temperature of all soil moisture
*           layers based on the simplified soil T profile described in Liang
*           et al. (1999), and used when QUICK_FLUX is TRUE.
*
* @note     These temperature estimates are much less accurate than those of
*           the finite element method (Cherkauer et al. (1999);
*           QUICK_FLUX FALSE). Since the Liang et al. (1999) approximation does
*           not account for the latent heat fluxes associated with freezing and
*           melting of ice, this function should not be called when FROZEN_SOIL
*           is TRUE.
******************************************************************************/
int
estimate_layer_temperature_quick_flux(layer_data_struct *layer,
                                      double            *depth,
                                      double             Dp,
                                      double             Tsurf,
                                      double             T1,
                                      double             Tp)
{
    extern option_struct options;

    size_t               lidx;
    double               Lsum[MAX_LAYERS + 1];

    // compute cumulative layer depths
    Lsum[0] = 0;
    for (lidx = 1; lidx <= options.Nlayer; lidx++) {
        Lsum[lidx] = depth[lidx - 1] + Lsum[lidx - 1];
    }

    // estimate soil layer average temperatures
    layer[0].T = 0.5 * (Tsurf + T1); // linear profile in topmost layer
    for (lidx = 1; lidx < options.Nlayer; lidx++) {
        layer[lidx].T = Tp - Dp / (depth[lidx]) *
                        (T1 -
                         Tp) *
                        (exp(-(Lsum[lidx + 1] -
                               Lsum[1]) /
                             Dp) - exp(-(Lsum[lidx] - Lsum[1]) / Dp));
    }

    return (0);
}

/******************************************************************************
* @brief    This subroutine estimates the ice content of all soil moisture
*           layers based on the simplified soil T profile described in Liang
*           et al. (1999), and used when QUICK_FLUX is TRUE.
*
* @note     This function should be called after estimate_layer_temperature_
*           quick_flux
******************************************************************************/
int
estimate_layer_ice_content_quick_flux(layer_data_struct *layer,
                                      double            *depth,
                                      double            *max_moist,
                                      double            *expt,
                                      double            *bubble,
                                      double            *frost_fract,
                                      double             frost_slope,
                                      char               FS_ACTIVE)
{
    extern option_struct options;

    size_t               lidx, frost_area;
    double               Lsum[MAX_LAYERS + 1];
    double               tmpT, tmp_fract, tmp_ice;
    double               min_temp, max_temp;

    // compute cumulative layer depths
    Lsum[0] = 0;
    for (lidx = 1; lidx <= options.Nlayer; lidx++) {
        Lsum[lidx] = depth[lidx - 1] + Lsum[lidx - 1];
    }

    // estimate soil layer ice contents
    for (lidx = 0; lidx < options.Nlayer; lidx++) {
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            layer[lidx].ice[frost_area] = 0;
        }

        if (options.FROZEN_SOIL && FS_ACTIVE) {
            min_temp = layer[lidx].T - frost_slope / 2.;
            max_temp = min_temp + frost_slope;
            for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
                if (frost_area == 0) {
                    tmp_fract = frost_fract[0] / 2.;
                }
                else {
                    tmp_fract +=
                        (frost_fract[frost_area - 1] /
                         2. + frost_fract[frost_area] / 2.);
                }
                tmpT = linear_interp(tmp_fract, 0, 1, min_temp, max_temp);
                tmp_ice = layer[lidx].moist -
                          maximum_unfrozen_water(tmpT, max_moist[lidx],
                                                 bubble[lidx], expt[lidx]);
                layer[lidx].ice[frost_area] = frost_fract[frost_area] * tmp_ice;
                if (layer[lidx].ice[frost_area] < 0) {
                    layer[lidx].ice[frost_area] = 0;
                }
                if (layer[lidx].ice[frost_area] > layer[lidx].moist) {
                    layer[lidx].ice[frost_area] = layer[lidx].moist;
                }
            }
        }
    }

    return (0);
}

/******************************************************************************
* @brief    This subroutine computes the thermal conductivity and volumetric
*           heat capacity of each soil layer based on its current moisture and
*           ice contents.  Ice is only present if the frozen soil algorithm is
*           activated.
******************************************************************************/
void
compute_soil_layer_thermal_properties(layer_data_struct *layer,
                                      double            *depth,
                                      double            *bulk_dens_min,
                                      double            *soil_dens_min,
                                      double            *quartz,
                                      double            *bulk_density,
                                      double            *soil_density,
                                      double            *organic,
                                      double            *frost_fract,
                                      size_t             Nlayers)
{
    extern option_struct options;
    size_t               lidx;
    size_t               frost_area;
    double               moist, ice;

    /* compute layer thermal properties */
    for (lidx = 0; lidx < Nlayers; lidx++) {
        moist = layer[lidx].moist / depth[lidx] / MM_PER_M;
        ice = 0;
        for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
            ice += layer[lidx].ice[frost_area] / depth[lidx] / MM_PER_M *
                   frost_fract[frost_area];
        }
        layer[lidx].kappa =
            soil_conductivity(moist, moist - ice,
                              soil_dens_min[lidx], bulk_dens_min[lidx],
                              quartz[lidx],
                              soil_density[lidx], bulk_density[lidx],
                              organic[lidx]);
        layer[lidx].Cs =
            volumetric_heat_capacity(bulk_density[lidx] / soil_density[lidx],
                                     moist - ice, ice, organic[lidx]);
    }
}

/******************************************************************************
* @brief    This subroutine reads through the soil thermal nodes and determines
*           the depths of all thawing and freezing fronts that are present.
******************************************************************************/
void
find_0_degree_fronts(energy_bal_struct *energy,
                     double            *Zsum_node,
                     double            *T,
                     int                Nnodes)
{
    int    nidx, fidx;
    int    Nthaw; /* number of thawing fronts found */
    int    Nfrost; /* number of frost fronts found */
    double tdepth[MAX_FRONTS]; /* thawing frost depths */
    double fdepth[MAX_FRONTS]; /* freezing front depths */

    /* Initialize parameters */
    Nthaw = Nfrost = 0;
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        fdepth[fidx] = MISSING;
        tdepth[fidx] = MISSING;
    }

    /* find 0 degree fronts */
    for (nidx = Nnodes - 2; nidx >= 0; nidx--) {
        if (T[nidx] > 0 && T[nidx + 1] <= 0 && Nthaw < MAX_FRONTS) {
            tdepth[Nthaw] =
                linear_interp(0, T[nidx], T[nidx + 1], Zsum_node[nidx],
                              Zsum_node[nidx + 1]);
            Nthaw++;
        }
        else if (T[nidx] < 0 && T[nidx + 1] >= 0 && Nfrost < MAX_FRONTS) {
            fdepth[Nfrost] =
                linear_interp(0, T[nidx], T[nidx + 1], Zsum_node[nidx],
                              Zsum_node[nidx + 1]);
            Nfrost++;
        }
    }

    /* store thaw depths */
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        energy->tdepth[fidx] = tdepth[fidx];
    }
    /* store frost depths */
    for (fidx = 0; fidx < MAX_FRONTS; fidx++) {
        energy->fdepth[fidx] = fdepth[fidx];
    }
    energy->Nthaw = Nthaw;
    energy->Nfrost = Nfrost;
}

/******************************************************************************
* @brief    This subroutine computes the maximum amount of unfrozen water that
*           can exist at the current temperature.
******************************************************************************/
double
maximum_unfrozen_water(double T,
                       double max_moist,
                       double bubble,
                       double expt)
{
    double unfrozen;

    if (T < 0.) {
        unfrozen = max_moist *
                   pow((-CONST_LATICE *
                        T) / (CONST_TKTRIP) / (CONST_G * bubble / (CM_PER_M)),
                       -(2.0 / (expt - 3.0)));
        if (unfrozen > max_moist) {
            unfrozen = max_moist;
        }
        if (unfrozen < 0.) {
            unfrozen = 0.;
        }
    }
    else {
        unfrozen = max_moist;
    }

    return (unfrozen);
}
