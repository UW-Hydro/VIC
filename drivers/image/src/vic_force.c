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
#include <vic_driver_image.h>

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
    float                     *fvar = NULL;
    size_t                     i;
    size_t                     j;
    size_t                     v;
    int                        vidx;
    size_t                     d3count[3];
    size_t                     d3start[3];

    // allocate memory for variables to be read
    fvar = (float *) malloc(local_domain.ncells * sizeof(float));
    if (fvar == NULL) {
        log_err("Memory allocation error in vic_force().");
    }

    // for now forcing file is determined by the year
    sprintf(filenames.forcing[0], "%s%4d.nc", filenames.f_path_pfx[0],
            dmy[current].year);

    // global_param.forceoffset[0] resets every year since the met file restarts
    // every year
    if (current > 1 && (dmy[current].year != dmy[current - 1].year)) {
        global_param.forceoffset[0] = 0;
    }

    // only the time slice changes for the met file reads. The rest is constant
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    // Air temperature: tas
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "tas",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].air_temp[j] = (double) fvar[i];
        }
    }

    // Precipitation: prcp
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "prcp",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].prec[j] = (double) fvar[i];
        }
    }

    // Downward solar radiation: dswrf
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "dswrf",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].shortwave[j] = (double) fvar[i];
        }
    }

    // Downward longwave radiation: dlwrf
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "dlwrf",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].longwave[j] = (double) fvar[i];
        }
    }

    // Wind speed: wind
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "wind",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].wind[j] = (double) fvar[i];
        }
    }

    // Specific humidity: shum
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "shum",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].vp[j] = (double) fvar[i];
        }
    }

    // Pressure: pressure
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceoffset[0] + j;
        get_scatter_nc_field_float(filenames.forcing[0], "pres",
                                   d3start, d3count, fvar);
        for (i = 0; i < local_domain.ncells; i++) {
            atmos[i].pressure[j] = (double) fvar[i];
        }
    }

    // Update the offset counter
    global_param.forceoffset[0] += NF;

    if (options.SNOW_BAND > 1) {
        log_err("SNOW_BAND not implemented in vic_force()");
    }
    else {
        t_offset = 0;
    }
    // Convert forcings into what we need and calculate missing ones
    for (i = 0; i < local_domain.ncells; i++) {
        for (j = 0; j < NF; j++) {
            // temperature in CONST_TKFRZ
            atmos[i].air_temp[j] -= CONST_TKFRZ;
            // precipitation in mm/period
            atmos[i].prec[j] *= global_param.snow_dt;
            // pressure in kPa
            atmos[i].pressure[j] /= PA_PER_KPA;
            // vapor pressure in kPa (we read specific humidity in kg/kg)
            atmos[i].vp[j] = q_to_vp(atmos[i].vp[j], atmos[i].pressure[j]);
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
            }
            // not the correct way to calculate average albedo, but leave
            // for now
            veg_hist[i][v].albedo[NR] = average(veg_hist[i][v].albedo, NF);
            veg_hist[i][v].LAI[NR] = average(veg_hist[i][v].LAI, NF);
            veg_hist[i][v].vegcover[NR] = average(veg_hist[i][v].vegcover, NF);
        }
    }


    // cleanup
    free(fvar);
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
