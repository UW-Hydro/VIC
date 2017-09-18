/******************************************************************************
 * @section DESCRIPTION
 *
 * Unpack forcing data.
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

/******************************************************************************
 * @brief    Read atmospheric forcing data.
 *****************************************************************************/
void
vic_force(void)
{
    extern size_t              NF;
    extern size_t              NR;
    extern size_t              current;
    extern force_data_struct  *force;
    extern x2l_data_struct    *x2l_vic;
    extern dmy_struct          dmy_current;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern veg_lib_struct    **veg_lib;
    extern parameters_struct   param;

    double                     t_offset;
    size_t                     i;
    size_t                     j;
    size_t                     v;
    int                        vidx;

    // Check to make sure variables have been set by coupler
    for (i = 0; i < local_domain.ncells_active; i++) {
        if (!x2l_vic[i].x2l_vars_set) {
            if (current == 0) {
                make_dummy_forcings(&x2l_vic[i]);
            }
            else {
                log_err("x2l_vars_set is false");
            }
        }
    }

    // Air temperature
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: K
            // VIC units: C
            force[i].air_temp[j] = x2l_vic[i].x2l_Sa_tbot - CONST_TKFRZ;
        }
    }

    // Precipitation
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: km m-2 s-1
            // VIC units: mm / timestep
            // Note: VIC does not use liquid/solid precip partitioning
            force[i].prec[j] = (x2l_vic[i].x2l_Faxa_rainc +
                                x2l_vic[i].x2l_Faxa_rainl +
                                x2l_vic[i].x2l_Faxa_snowc +
                                x2l_vic[i].x2l_Faxa_snowl) *
                               global_param.snow_dt;
        }
    }

    // Downward solar radiation
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: W m-2
            // VIC units: W m-2
            // Note: VIC does not use partitioned shortwave fluxes.
            force[i].shortwave[j] = (x2l_vic[i].x2l_Faxa_swndr +
                                     x2l_vic[i].x2l_Faxa_swvdr +
                                     x2l_vic[i].x2l_Faxa_swndf +
                                     x2l_vic[i].x2l_Faxa_swvdf);
        }
    }

    // Downward longwave radiation
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: W m-2
            // VIC units: W m-2
            force[i].longwave[j] = x2l_vic[i].x2l_Faxa_lwdn;
        }
    }

    // Wind speed
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: m s-1
            // VIC units: m s-1
            // Note: VIC does not use partitioned wind speeds
            force[i].wind[j] = sqrt(pow(x2l_vic[i].x2l_Sa_u, 2) +
                                    pow(x2l_vic[i].x2l_Sa_v, 2));
        }
    }

    // Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: Pa
            // VIC units: Pa
            // Note: Image Driver uses kPa inputs and
            // converts to Pa
            force[i].pressure[j] = x2l_vic[i].x2l_Sa_pbot;
        }
    }

    // Vapor Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: shum is specific humidity (g/g)
            // VIC units: Pa
            force[i].vp[j] = q_to_vp(x2l_vic[i].x2l_Sa_shum,
                                     force[i].pressure[j]);
        }
    }

    if (options.CARBON) {
        // Fraction of incoming shortwave that is direct
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: n/a (calculated from SW fluxes)
                // VIC units: fraction
                if (force[i].shortwave[j] != 0.) {
                    force[i].fdir[j] = (x2l_vic[i].x2l_Faxa_swndr +
                                        x2l_vic[i].x2l_Faxa_swvdr) /
                                       (x2l_vic[i].x2l_Faxa_swndf +
                                        x2l_vic[i].x2l_Faxa_swvdf);
                }
                else {
                    force[i].fdir[j] = 0.;
                }
            }
        }

        // Concentration of CO2
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: 1e-6 mol/mol
                // VIC units: mol CO2/ mol air
                force[i].Catm[j] = 1e6 * x2l_vic[i].x2l_Sa_co2prog;
            }
        }

        // Cosine of solar zenith angle
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                force[i].coszen[j] = compute_coszen(
                    local_domain.locations[i].latitude,
                    local_domain.locations[i].longitude,
                    soil_con[i].time_zone_lng, dmy_current.day_in_year,
                    dmy_current.dayseconds);
            }
        }
    }

    if (options.LAKES) {
        // incoming channel inflow
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: kg m-2 s-1
                // VIC units: mm
                force[i].channel_in[j] = x2l_vic[i].x2l_Flrr_flood *
                                         global_param.snow_dt;
            }
        }
    }

    if (options.SNOW_BAND > 1) {
        log_err("SNOW_BAND not implemented");
    }
    else {
        t_offset = 0;
    }

    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < NF; j++) {
            // vapor pressure deficit
            force[i].vpd[j] = svp(force[i].air_temp[j]) - force[i].vp[j];
            if (force[i].vpd[j] < 0) {
                log_warn("Vapor pressure deficit is %f which is < 0, "
                         "setting vapor pressure deficit to 0 and calculating "
                         "saturated vapor pressure using air temperature %f.",
                         force[i].vpd[j], force[i].air_temp[j]);
                force[i].vpd[j] = 0;
                force[i].vp[j] = svp(force[i].air_temp[j]);
            }
            // photosynthetically active radiation
            // TODO: Add CARBON_SW2PAR back to the parameters structure
            // force[i].par[j] = param.CARBON_SW2PAR * force[i].shortwave[j];
            // air density
            force[i].density[j] = air_density(force[i].air_temp[j],
                                              force[i].pressure[j]);
            // snow flag
            force[i].snowflag[j] = will_it_snow(&(force[i].air_temp[j]),
                                                t_offset,
                                                param.SNOW_MAX_SNOW_TEMP,
                                                &(force[i].prec[j]), 1);
        }
    }


    // Put average value in NR field
    for (i = 0; i < local_domain.ncells_active; i++) {
        force[i].air_temp[NR] = average(force[i].air_temp, NF);
        // For precipitation put total
        force[i].prec[NR] = average(force[i].prec, NF) * NF;
        force[i].shortwave[NR] = average(force[i].shortwave, NF);
        force[i].longwave[NR] = average(force[i].longwave, NF);
        force[i].pressure[NR] = average(force[i].pressure, NF);
        force[i].wind[NR] = average(force[i].wind, NF);
        force[i].vp[NR] = average(force[i].vp, NF);
        force[i].vpd[NR] = (svp(force[i].air_temp[NR]) - force[i].vp[NR]);
        force[i].density[NR] = air_density(force[i].air_temp[NR],
                                           force[i].pressure[NR]);
        force[i].snowflag[NR] = will_it_snow(force[i].air_temp, t_offset,
                                             param.SNOW_MAX_SNOW_TEMP,
                                             force[i].prec, NF);

        // Optional inputs
        if (options.LAKES) {
            force[i].channel_in[NR] = average(force[i].channel_in, NF) * NF;
        }
        if (options.CARBON) {
            force[i].Catm[NR] = average(force[i].Catm, NF);
            force[i].fdir[NR] = average(force[i].fdir, NF);
            force[i].par[NR] = average(force[i].par, NF);
            // for coszen, use value at noon
            force[i].coszen[NR] = compute_coszen(
                local_domain.locations[i].latitude,
                local_domain.locations[i].longitude, soil_con[i].time_zone_lng,
                dmy_current.day_in_year, SEC_PER_DAY / 2);
        }
    }

    // Update the veg_hist structure with the current vegetation parameters.
    // Currently only implemented for climatological values in image mode
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    veg_hist[i][vidx].albedo[j] =
                        veg_con[i][vidx].albedo[dmy_current.month - 1];
                    veg_hist[i][vidx].displacement[j] =
                        veg_con[i][vidx].displacement[dmy_current.month - 1];
                    veg_hist[i][vidx].fcanopy[j] =
                        veg_con[i][vidx].fcanopy[dmy_current.month - 1];
                    veg_hist[i][vidx].LAI[j] =
                        veg_con[i][vidx].LAI[dmy_current.month - 1];
                    veg_hist[i][vidx].roughness[j] =
                        veg_con[i][vidx].roughness[dmy_current.month - 1];
                }
                // not the correct way to calculate average albedo, but leave
                // for now
                veg_hist[i][vidx].albedo[NR] = average(veg_hist[i][vidx].albedo,
                                                       NF);
                veg_hist[i][vidx].displacement[NR] = average(
                    veg_hist[i][vidx].displacement, NF);
                veg_hist[i][vidx].fcanopy[NR] = average(
                    veg_hist[i][vidx].fcanopy, NF);
                veg_hist[i][vidx].LAI[NR] = average(veg_hist[i][vidx].LAI, NF);
                veg_hist[i][vidx].roughness[NR] = average(
                    veg_hist[i][vidx].roughness, NF);
            }
        }
    }
}

/******************************************************************************
 * @brief   dummy forcings for initialization (should be removed or never used)
 *****************************************************************************/
void
make_dummy_forcings(x2l_data_struct *x2l)
{
    extern x2l_data_struct *x2l_vic;
    extern domain_struct    local_domain;

    x2l->x2l_Sa_z = 10;  /** bottom atm level height */
    x2l->x2l_Sa_u = 1.;  /** bottom atm level zon wind */
    x2l->x2l_Sa_v = 1.;  /** bottom atm level mer wind */
    x2l->x2l_Sa_ptem = 1.;  /** bottom atm level pot temp */
    x2l->x2l_Sa_shum = 0.02;  /** bottom atm level spec hum */
    x2l->x2l_Sa_pbot = 101325.;  /** bottom atm level pressure */
    x2l->x2l_Sa_tbot = 275.0;  /** bottom atm level temp */
    x2l->x2l_Faxa_lwdn = 50.;  /** downward lw heat flux */
    x2l->x2l_Faxa_rainc = 0.;  /** prec: liquid "convective" */
    x2l->x2l_Faxa_rainl = 0.;  /** prec: liquid "large scale" */
    x2l->x2l_Faxa_snowc = 0.;  /** prec: frozen "convective" */
    x2l->x2l_Faxa_snowl = 0.;  /** prec: frozen "large scale" */
    x2l->x2l_Faxa_swndr = 1.;  /** sw: nir direct  downward */
    x2l->x2l_Faxa_swvdr = 1.;  /** sw: vis direct  downward */
    x2l->x2l_Faxa_swndf = 1.;  /** sw: nir diffuse downward */
    x2l->x2l_Faxa_swvdf = 1.;  /** sw: vis diffuse downward */
    x2l->x2l_Sa_co2prog = 0.;  /** bottom atm level prognostic co2 */
    x2l->x2l_Sa_co2diag = 0.;  /** bottom atm level diagnostic co2 */
    x2l->x2l_Faxa_bcphidry = 0.;  /** flux: Black Carbon hydrophilic dry deposition */
    x2l->x2l_Faxa_bcphodry = 0.;  /** flux: Black Carbon hydrophobic dry deposition */
    x2l->x2l_Faxa_bcphiwet = 0.;  /** flux: Black Carbon hydrophilic wet deposition */
    x2l->x2l_Faxa_ocphidry = 0.;  /** flux: Organic Carbon hydrophilic dry deposition */
    x2l->x2l_Faxa_ocphodry = 0.;  /** flux: Organic Carbon hydrophobic dry deposition */
    x2l->x2l_Faxa_ocphiwet = 0.;  /** flux: Organic Carbon hydrophilic dry deposition */
    x2l->x2l_Faxa_dstwet1 = 0.;  /** flux: Size 1 dust -- wet deposition */
    x2l->x2l_Faxa_dstwet2 = 0.;  /** flux: Size 2 dust -- wet deposition */
    x2l->x2l_Faxa_dstwet3 = 0.;  /** flux: Size 3 dust -- wet deposition */
    x2l->x2l_Faxa_dstwet4 = 0.;  /** flux: Size 4 dust -- wet deposition */
    x2l->x2l_Faxa_dstdry1 = 0.;  /** flux: Size 1 dust -- dry deposition */
    x2l->x2l_Faxa_dstdry2 = 0.;  /** flux: Size 2 dust -- dry deposition */
    x2l->x2l_Faxa_dstdry3 = 0.;  /** flux: Size 3 dust -- dry deposition */
    x2l->x2l_Faxa_dstdry4 = 0.;  /** flux: Size 4 dust -- dry deposition */
    x2l->x2l_Flrr_flood = 0.;  /** rtm->lnd rof (flood) flux */
    x2l->x2l_vars_set = true; /** x2l set flag */
}
