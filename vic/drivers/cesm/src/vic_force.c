/******************************************************************************
 * @section DESCRIPTION
 *
 * Read atmospheric forcing data.
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
    extern atmos_data_struct  *atmos;
    extern x2l_data_struct    *x2l_vic;
    extern dmy_struct          dmy;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern filenames_struct    filenames;
    extern global_param_struct global_param;
    extern option_struct       options;
    extern veg_con_map_struct *veg_con_map;
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
            atmos[i].air_temp[j] = x2l_vic[i].x2l_Sa_tbot - CONST_TKFRZ;
        }
    }

    // Precipitation
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: km m-2 s-1
            // VIC units: mm / timestep
            // Note: VIC does not use liquid/solid precip partitioning
            atmos[i].prec[j] = (x2l_vic[i].x2l_Faxa_rainc +
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
            atmos[i].shortwave[j] = (x2l_vic[i].x2l_Faxa_swndr +
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
            atmos[i].longwave[j] = x2l_vic[i].x2l_Faxa_lwdn;
        }
    }

    // Wind speed
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: m s-1
            // VIC units: m s-1
            // Note: VIC does not use partitioned wind speeds
            atmos[i].wind[j] = sqrt(pow(x2l_vic[i].x2l_Sa_u, 2) +
                                    pow(x2l_vic[i].x2l_Sa_v, 2));
        }
    }

    // Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: Pa
            // VIC units: kPa
            atmos[i].pressure[j] = x2l_vic[i].x2l_Sa_pbot / PA_PER_KPA;
        }
    }

    // Vapor Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            // CESM units: shum is specific humidity (g/g)
            // VIC units: kPa
            atmos[i].vp[j] = q_to_vp(x2l_vic[i].x2l_Sa_shum,
                                     atmos[i].pressure[j]);
        }
    }

    if (options.CARBON) {
        // Fraction of incoming shortwave that is direct
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: n/a (calculated from SW fluxes)
                // VIC units: fraction
                if (atmos[i].shortwave[j] != 0.) {
                    atmos[i].fdir[j] = (x2l_vic[i].x2l_Faxa_swndr +
                                        x2l_vic[i].x2l_Faxa_swvdr) /
                                       (x2l_vic[i].x2l_Faxa_swndf +
                                        x2l_vic[i].x2l_Faxa_swvdf);
                }
                else {
                    atmos[i].fdir[j] = 0.;
                }
            }
        }

        // Concentration of CO2
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: 1e-6 mol/mol
                // VIC units: mol CO2/ mol air
                atmos[i].Catm[j] = 1e6 * x2l_vic[i].x2l_Sa_co2prog;
            }
        }
    }

    if (options.LAKES) {
        // incoming channel inflow
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                // CESM units: kg m-2 s-1
                // VIC units: mm
                atmos[i].channel_in[j] = x2l_vic[i].x2l_Flrr_flood *
                                         global_param.snow_dt;
            }
        }
    }

    if (options.SNOW_BAND > 1) {
        log_err("SNOW_BAND not implemented in vic_force()");
    }
    else {
        t_offset = 0;
    }

    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < NF; j++) {
            // vapor pressure deficit
            atmos[i].vpd[j] = svp(atmos[i].air_temp[j]) - atmos[i].vp[j];
            // photosynthetically active radiation
            // TODO: Add CARBON_SW2PAR back to the parameters structure
            // atmos[i].par[j] = param.CARBON_SW2PAR * atmos[i].shortwave[j];
            // air density
            atmos[i].density[j] = air_density(atmos[i].air_temp[j],
                                              atmos[i].pressure[j]);
            // snow flag
            atmos[i].snowflag[j] = will_it_snow(&(atmos[i].air_temp[j]),
                                                t_offset,
                                                param.SNOW_MAX_SNOW_TEMP,
                                                &(atmos[i].prec[j]), 1);
        }
    }


    // Put average value in NR field
    for (i = 0; i < local_domain.ncells_active; i++) {
        atmos[i].air_temp[NR] = average(atmos[i].air_temp, NF);
        // For precipitation put total
        atmos[i].prec[NR] = average(atmos[i].prec, NF) * NF;
        atmos[i].shortwave[NR] = average(atmos[i].shortwave, NF);
        atmos[i].longwave[NR] = average(atmos[i].longwave, NF);
        atmos[i].pressure[NR] = average(atmos[i].pressure, NF);
        atmos[i].wind[NR] = average(atmos[i].wind, NF);
        atmos[i].vp[NR] = average(atmos[i].vp, NF);
        atmos[i].vpd[NR] = (svp(atmos[i].air_temp[NR]) - atmos[i].vp[NR]);
        atmos[i].density[NR] = air_density(atmos[i].air_temp[NR],
                                           atmos[i].pressure[NR]);
        atmos[i].snowflag[NR] = will_it_snow(atmos[i].air_temp, t_offset,
                                             param.SNOW_MAX_SNOW_TEMP,
                                             atmos[i].prec, NF);

        // Optional inputs
        if (options.LAKES) {
            atmos[i].channel_in[NR] = average(atmos[i].channel_in, NF) * NF;
        }
        if (options.CARBON) {
            atmos[i].Catm[NR] = average(atmos[i].Catm, NF);
            atmos[i].fdir[NR] = average(atmos[i].fdir, NF);
            atmos[i].par[NR] = average(atmos[i].par, NF);
        }
    }

    // TBD: coszen (used for some of the carbon functions)

    // Update the veg_hist structure with the current vegetation parameters.
    // Currently only implemented for climatological values in image mode
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    veg_hist[i][vidx].albedo[j] =
                        veg_lib[i][v].albedo[dmy.month - 1];
                    veg_hist[i][vidx].LAI[j] =
                        veg_lib[i][v].LAI[dmy.month - 1];
                    veg_hist[i][vidx].fcanopy[j] =
                        veg_lib[i][v].fcanopy[dmy.month - 1];
                }
                // not the correct way to calculate average albedo, but leave
                // for now
                veg_hist[i][vidx].albedo[NR] = average(veg_hist[i][vidx].albedo,
                                                       NF);
                veg_hist[i][vidx].LAI[NR] = average(veg_hist[i][vidx].LAI, NF);
                veg_hist[i][vidx].fcanopy[NR] = average(
                    veg_hist[i][vidx].fcanopy, NF);
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
