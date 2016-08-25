/******************************************************************************
 * @section DESCRIPTION
 *
 * Collect coupled fields and export to coupler
 *
 * Sign convention: positive value <=> downward flux
 * Units: see notes in code or seq_flds_mod.F90
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

#include <vic_driver_cesm.h>

void
vic_cesm_put_data()
{
    extern all_vars_struct    *all_vars;
    extern force_data_struct  *force;
    extern dmy_struct          dmy_current;
    extern domain_struct       local_domain;
    extern soil_con_struct    *soil_con;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;
    extern l2x_data_struct    *l2x_vic;
    extern x2l_data_struct    *x2l_vic;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern parameters_struct   param;

    bool                       IsWet = false; // TODO: add lake fraction
    bool                       overstory;
    bool                       HasVeg;
    size_t                     i;
    size_t                     veg;
    size_t                     band;
    size_t                     index;
    double                     AreaFactor;
    double                     AreaFactorSum;
    double                     TreeAdjustFactor = 1.;
    double                     lakefactor = 1.;
    double                     rad_temp;
    double                     albedo;
    double                     aero_resist;
    double                     roughness;
    double                     wind_stress;
    double                     wind_stress_x;
    double                     wind_stress_y;
    double                     evap;
    cell_data_struct           cell;
    energy_bal_struct          energy;
    snow_data_struct           snow;
    veg_var_struct             veg_var;

    for (i = 0; i < local_domain.ncells_active; i++) {
        // Zero l2x vars (leave unused fields as MISSING values)
        l2x_vic[i].l2x_Sl_t = 0;
        l2x_vic[i].l2x_Sl_tref = 0;
        l2x_vic[i].l2x_Sl_qref = 0;
        l2x_vic[i].l2x_Sl_avsdr = 0;
        l2x_vic[i].l2x_Sl_anidr = 0;
        l2x_vic[i].l2x_Sl_avsdf = 0;
        l2x_vic[i].l2x_Sl_anidf = 0;
        l2x_vic[i].l2x_Sl_snowh = 0;
        l2x_vic[i].l2x_Sl_u10 = 0;
        // l2x_vic[i].l2x_Sl_ddvel = 0;
        l2x_vic[i].l2x_Sl_fv = 0;
        l2x_vic[i].l2x_Sl_ram1 = 0;
        l2x_vic[i].l2x_Sl_logz0 = 0;
        l2x_vic[i].l2x_Fall_taux = 0;
        l2x_vic[i].l2x_Fall_tauy = 0;
        l2x_vic[i].l2x_Fall_lat = 0;
        l2x_vic[i].l2x_Fall_sen = 0;
        l2x_vic[i].l2x_Fall_lwup = 0;
        l2x_vic[i].l2x_Fall_evap = 0;
        l2x_vic[i].l2x_Fall_swnet = 0;
        // l2x_vic[i].l2x_Fall_fco2_lnd = 0;
        // l2x_vic[i].l2x_Fall_flxdst1 = 0;
        // l2x_vic[i].l2x_Fall_flxdst2 = 0;
        // l2x_vic[i].l2x_Fall_flxdst3 = 0;
        // l2x_vic[i].l2x_Fall_flxdst4 = 0;
        // l2x_vic[i].l2x_Fall_flxvoc = 0;
        l2x_vic[i].l2x_Flrl_rofliq = 0;
        // l2x_vic[i].l2x_Flrl_rofice = 0;

        // running sum to make sure we get the full grid cell
        AreaFactorSum = 0;

        for (veg = 0; veg <= local_domain.locations[i].nveg; veg++) {
            overstory = veg_lib[i][veg_con[i][veg].veg_class].overstory;
            if (veg <= local_domain.locations[i].nveg - 1) {
                HasVeg = true;
            }
            else {
                HasVeg = false;
            }

            for (band = 0; band < options.SNOW_BAND; band++) {
                cell = all_vars[i].cell[veg][band];
                energy = all_vars[i].energy[veg][band];
                snow = all_vars[i].snow[veg][band];
                veg_var = all_vars[i].veg_var[veg][band];

                // TODO: Consider treeline and lake factors
                AreaFactor = (veg_con[i][veg].Cv *
                              soil_con[i].AreaFract[band] *
                              TreeAdjustFactor * lakefactor);
                if (AreaFactor < DBL_EPSILON) {
                    // Skip this patch since the area factor is zero
                    continue;
                }
                AreaFactorSum += AreaFactor;

                // temperature
                // CESM units: K
                if (overstory && snow.snow && !(options.LAKES && IsWet)) {
                    rad_temp = energy.Tfoliage + CONST_TKFRZ;
                }
                else {
                    rad_temp = energy.Tsurf + CONST_TKFRZ;
                }
                l2x_vic[i].l2x_Sl_t += AreaFactor * rad_temp;

                // 2m reference temperature
                // CESM units: K
                l2x_vic[i].l2x_Sl_tref += AreaFactor * force->air_temp[NR];

                // 2m reference specific humidity
                // CESM units: g/g
                l2x_vic[i].l2x_Sl_qref += AreaFactor * CONST_EPS *
                                          force->vp[NR] / force->pressure[NR];

                // Albedo Note: VIC does not partition its albedo, all returned
                // values will be the same

                // albedo: direct, visible
                // CESM units: unitless
                // force->shortwave is the incoming shortwave (+ down)
                // force->NetShortAtmos net shortwave flux (+ down)
                // SWup = force->shortwave[NR] - energy.NetShortAtmos
                // Set the albedo to zero for the case where there is no shortwave down
                if (force->shortwave[NR] > 0.) {
                    albedo = AreaFactor *
                             (force->shortwave[NR] - energy.NetShortAtmos) /
                             force->shortwave[NR];
                }
                else {
                    albedo = 0.;
                }
                l2x_vic[i].l2x_Sl_avsdr += albedo;

                // albedo: direct , near-ir
                // CESM units: unitless
                l2x_vic[i].l2x_Sl_anidr += albedo;

                // albedo: diffuse, visible
                // CESM units: unitless
                l2x_vic[i].l2x_Sl_avsdf += albedo;

                // albedo: diffuse, near-ir
                // CESM units: unitless
                l2x_vic[i].l2x_Sl_anidf += albedo;

                // snow height
                // CESM units: m
                l2x_vic[i].l2x_Sl_snowh += AreaFactor * snow.depth;

                // 10m wind
                // CESM units: m/s
                l2x_vic[i].l2x_Sl_u10 += AreaFactor * force->wind[NR];

                // dry deposition velocities (optional)
                // CESM units: ?
                // l2x_vic[i].l2x_Sl_ddvel;

                // aerodynamical resistance
                // CESM units: s/m
                if (overstory) {
                    aero_resist = cell.aero_resist[1];
                }
                else {
                    aero_resist = cell.aero_resist[0];
                }

                if (aero_resist < DBL_EPSILON) {
                    log_warn("aero_resist (%f) is < %f", aero_resist,
                             DBL_EPSILON);
                    aero_resist = param.HUGE_RESIST;
                }

                l2x_vic[i].l2x_Sl_ram1 += AreaFactor * aero_resist;

                // log z0
                // CESM units: m
                if (snow.snow) {
                    // snow roughness
                    roughness = soil_con[i].snow_rough;
                }
                else if (HasVeg) {
                    // bare soil roughness
                    roughness =
                        veg_lib[i][veg_con[i][veg].veg_class].roughness[
                            dmy_current.month - 1];
                }
                else {
                    roughness = soil_con[i].rough;
                }
                if (roughness < DBL_EPSILON) {
                    log_warn("roughness (%f) is < %f", roughness, DBL_EPSILON);
                    roughness = DBL_EPSILON;
                }
                l2x_vic[i].l2x_Sl_logz0 += AreaFactor * log(roughness);

                // wind stress, zonal
                // CESM units: N m-2
                wind_stress_x = -1 * force[i].density[NR] *
                                x2l_vic[i].x2l_Sa_u / aero_resist;
                l2x_vic[i].l2x_Fall_taux += AreaFactor * wind_stress_x;

                // wind stress, meridional
                // CESM units: N m-2
                wind_stress_y = -1 * force[i].density[NR] *
                                x2l_vic[i].x2l_Sa_v / aero_resist;
                l2x_vic[i].l2x_Fall_tauy += AreaFactor * wind_stress_y;

                // friction velocity
                // CESM units: m s-1
                wind_stress =
                    sqrt(pow(wind_stress_x, 2) + pow(wind_stress_y, 2));
                l2x_vic[i].l2x_Sl_fv += AreaFactor *
                                        (wind_stress / force[i].density[NR]);

                // latent heat flux
                // CESM units: W m-2
                l2x_vic[i].l2x_Fall_lat += -1 * AreaFactor * energy.AtmosLatent;

                // sensible heat flux
                // CESM units: W m-2
                l2x_vic[i].l2x_Fall_sen += -1 * AreaFactor *
                                           energy.AtmosSensible;

                // upward longwave heat flux
                // CESM units: W m-2
                l2x_vic[i].l2x_Fall_lwup += AreaFactor *
                                            (force->longwave[NR] -
                                             energy.NetLongAtmos);

                // evaporation water flux
                // CESM units: kg m-2 s-1
                evap = 0.0;
                for (index = 0; index < options.Nlayer; index++) {
                    evap += cell.layer[index].evap;
                }
                evap += snow.vapor_flux * MM_PER_M;
                if (HasVeg) {
                    evap += snow.canopy_vapor_flux * MM_PER_M;
                    evap += veg_var.canopyevap;
                }
                l2x_vic[i].l2x_Fall_evap += -1 * AreaFactor * evap /
                                            global_param.dt;

                // heat flux shortwave net
                l2x_vic[i].l2x_Fall_swnet += AreaFactor *
                                             (force->shortwave[NR] -
                                              energy.NetShortAtmos);

                // co2 flux **For testing set to 0
                // l2x_vic[i].l2x_Fall_fco2_lnd;

                // dust flux size bin 1
                // l2x_vic[i].l2x_Fall_flxdst1;

                // dust flux size bin 2
                // l2x_vic[i].l2x_Fall_flxdst2;

                // dust flux size bin 3
                // l2x_vic[i].l2x_Fall_flxdst3;

                // dust flux size bin 4
                // l2x_vic[i].l2x_Fall_flxdst4;

                // MEGAN fluxes
                // l2x_vic[i].l2x_Fall_flxvoc;

                // lnd->rtm input fluxes
                l2x_vic[i].l2x_Flrl_rofliq += AreaFactor *
                                              (cell.runoff +
                                               cell.baseflow) / global_param.dt;

                // lnd->rtm input fluxes
                // l2x_vic[i].l2x_Flrl_rofice;

                // vars set flag
                l2x_vic[i].l2x_vars_set = true;
            }
        }

        if (!assert_close_double(AreaFactorSum, 1., 0., 1e-3)) {
            log_warn("AreaFactorSum (%f) is not 1",
                     AreaFactorSum);
        }
    }
}
