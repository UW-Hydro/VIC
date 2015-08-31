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

#include <vic_def.h>
#include <vic_run.h>
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
    extern dmy_struct         *dmy;
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

    // Air temperature
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells; i++) {
            // CESM units: K
            // VIC units: C
            atmos[i].air_temp[j] = x2l_vic[i].x2l_Sa_tbot - CONST_TKFRZ;
        }
    }

    // Precipitation
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells; i++) {
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
        for (i = 0; i < local_domain.ncells; i++) {
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
        for (i = 0; i < local_domain.ncells; i++) {
            // CESM units: W m-2
            // VIC units: W m-2
            atmos[i].longwave[j] = x2l_vic[i].x2l_Faxa_lwdn;
        }
    }

    // Wind speed
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells; i++) {
            // CESM units: m s-1
            // VIC units: m s-1
            // Note: VIC does not use partitioned wind speeds
            atmos[i].wind[j] = sqrt(pow(x2l_vic[i].x2l_Sa_u, 2) +
                                    pow(x2l_vic[i].x2l_Sa_v, 2));
        }
    }

    // Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells; i++) {
            // CESM units: Pa
            // VIC units: kPa
            atmos[i].pressure[j] = x2l_vic[i].x2l_Sa_pbot / PA_PER_KPA;
        }
    }

    // Vapor Pressure
    for (j = 0; j < NF; j++) {
        for (i = 0; i < local_domain.ncells; i++) {
            // CESM units: shum is specific humidity (g/g)
            // VIC units: kPa
            atmos[i].vp[j] = q_to_vp(x2l_vic[i].x2l_Sa_shum,
                                     atmos[i].pressure[j]);
        }
    }

    if (options.SNOW_BAND > 1) {
        log_err("SNOW_BAND not implemented in vic_force()");
    }
    else {
        t_offset = 0;
    }

    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells; i++) {
        for (j = 0; j < NF; j++) {
            // vapor pressure deficit
            atmos[i].vpd[j] = svp(atmos[i].air_temp[j]) - atmos[i].vp[j];
            // photosynthetically active radiation
            atmos[i].par[j] = param.CARBON_SW2PAR * atmos[i].shortwave[j];
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
    for (i = 0; i < local_domain.ncells; i++) {
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
    }

    // TBD: coszen (used for some of the carbon functions), fdir (if needed)
    // Catm, fdir (not used as far as I can tell)

    // Update the veg_hist structure with the current vegetation parameters.
    // Currently only implemented for climatological values in image mode
    for (i = 0; i < local_domain.ncells; i++) {
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != -1) {
                for (j = 0; j < NF; j++) {
                    veg_hist[i][vidx].albedo[j] =
                        veg_lib[i][v].albedo[dmy[current].month - 1];
                    veg_hist[i][vidx].LAI[j] =
                        veg_lib[i][v].LAI[dmy[current].month - 1];
                    veg_hist[i][vidx].vegcover[j] =
                        veg_lib[i][v].vegcover[dmy[current].month - 1];
                }
                // not the correct way to calculate average albedo, but leave
                // for now
                veg_hist[i][vidx].albedo[NR] = average(veg_hist[i][vidx].albedo,
                                                       NF);
                veg_hist[i][vidx].LAI[NR] = average(veg_hist[i][vidx].LAI, NF);
                veg_hist[i][vidx].vegcover[NR] = average(
                    veg_hist[i][vidx].vegcover, NF);
            }
        }
    }
}

/******************************************************************************
 * @brief    calculate 1d average
 *****************************************************************************/
double
average(double *ar,
        size_t  n)
{
    size_t i;
    double sum = 0.;

    if (n <= 0) {
        log_err("Error in calc_average: divide by zero or negative");
    }
    else if (n == 1) {
        return ar[0];
    }
    else {
        for (i = 0; i < n; i++) {
            sum += ar[i];
        }
    }

    return sum / n;
}

/******************************************************************************
 * @brief   convert specific humidity (q) to vapor pressure (vp) based on
 *          pressure (p)
 *
 * @param q specific humidity
 * @param p pressure
 *
 * @return vp vapor pressure (units are the same as p)
 *****************************************************************************/
double
q_to_vp(double q,
        double p)
{
    double vp;

    // full equation
    // vp = q/(q+CONST_EPS*(1-q))*p;

    // approximation used in VIC
    vp = q * p / CONST_EPS;

    return vp;
}

/******************************************************************************
 * @brief   convert surface pressure (kPa) to density (kg/m3) based on
 *          pressure (p), vapor pressure (vp), and temperature
 *
 * @param t temperature
 * @param p pressure
 *
 * @return rho surface pressure
 *****************************************************************************/
double
air_density(double t,
            double p)
{
    double rho;

    // full equation
    // rho = (p*1000)/(Rd * *t+CONST_TKFRZ) + (pv*1000)/(Rv * *t+CONST_TKFRZ);

    // approximation used in VIC
    rho = 0.003486 * p / (275.0 + t);

    return rho;
}

/******************************************************************************
 * @brief   return 1 if it will snow, otherwise return 0
 *****************************************************************************/
char
will_it_snow(double *t,
             double  t_offset,
             double  max_snow_temp,
             double *prcp,
             size_t  n)
{
    size_t i;

    for (i = 0; i < n; i++) {
        if ((t[i] + t_offset) < max_snow_temp && prcp[i] > 0.) {
            return 1;
        }
    }

    return 0;
}
