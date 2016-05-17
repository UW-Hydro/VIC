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
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_hist_struct   **veg_hist;
    extern parameters_struct   param;
    extern param_set_struct    param_set;

    double                     t_offset;
    double                    *dvar = NULL;
    size_t                     i;
    size_t                     j;
    size_t                     v;
    int                        vidx;
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];

    // allocate memory for variables to be read
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
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
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[AIR_TEMP].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].air_temp[j] = (double) dvar[i];
        }
    }

    // Precipitation: prcp
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[PREC].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].prec[j] = (double) dvar[i];
        }
    }

    // Downward solar radiation: dswrf
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[SWDOWN].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].shortwave[j] = (double) dvar[i];
        }
    }

    // Downward longwave radiation: dlwrf
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[LWDOWN].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].longwave[j] = (double) dvar[i];
        }
    }

    // Wind speed: wind
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[WIND].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].wind[j] = (double) dvar[i];
        }
    }

    // Specific humidity: shum
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[VP].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].vp[j] = (double) dvar[i];
        }
    }

    // Pressure: pressure
    for (j = 0; j < NF; j++) {
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[PRESSURE].varname,
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            atmos[i].pressure[j] = (double) dvar[i];
        }
    }
    // Optional inputs
    if (options.LAKES) {
        // Channel inflow to lake
        d3start[0] = global_param.forceskip[0] + global_param.forceoffset[0] +
                     j;
        get_scatter_nc_field_double(filenames.forcing[0],
                                    param_set.TYPE[CHANNEL_IN].varname,
                                    d3start, d3count, dvar);
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                atmos[i].channel_in[j] = (double) dvar[i];
            }
        }
    }
    if (options.CARBON) {
        // Atmospheric CO2 mixing ratio
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(filenames.forcing[0],
                                        param_set.TYPE[CATM].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                atmos[i].Catm[j] = (double) dvar[i];
            }
        }
        // Cosine of solar zenith angle
        for (j = 0; j < NF; j++) {
            for (i = 0; i < local_domain.ncells_active; i++) {
                atmos[i].coszen[j] = compute_coszen(
                    local_domain.locations[i].latitude,
                    local_domain.locations[i].longitude,
                    soil_con[i].time_zone_lng, dmy[current].day_in_year,
                    dmy[current].dayseconds);
            }
        }
        // Fraction of shortwave that is direct
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(filenames.forcing[0],
                                        param_set.TYPE[FDIR].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                atmos[i].fdir[j] = (double) dvar[i];
            }
        }
        // Photosynthetically active radiation
        for (j = 0; j < NF; j++) {
            d3start[0] = global_param.forceskip[0] +
                         global_param.forceoffset[0] + j;
            get_scatter_nc_field_double(filenames.forcing[0],
                                        param_set.TYPE[PAR].varname,
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                atmos[i].par[j] = (double) dvar[i];
            }
        }
    }

    // Update the offset counter
    global_param.forceoffset[0] += NF;

    // Initialize the veg_hist structure with the current climatological
    // vegetation parameters.  This may be overwritten with the historical
    // forcing time series.
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    veg_hist[i][vidx].albedo[j] =
                        veg_con[i][vidx].albedo[dmy[current].month - 1];
                    veg_hist[i][vidx].displacement[j] =
                        veg_con[i][vidx].displacement[dmy[current].month - 1];
                    veg_hist[i][vidx].fcanopy[j] =
                        veg_con[i][vidx].fcanopy[dmy[current].month - 1];
                    veg_hist[i][vidx].LAI[j] =
                        veg_con[i][vidx].LAI[dmy[current].month - 1];
                    veg_hist[i][vidx].roughness[j] =
                        veg_con[i][vidx].roughness[dmy[current].month - 1];
                }
            }
        }
    }

    // Read veg_hist file
    if (options.LAI_SRC == FROM_VEGHIST ||
        options.FCAN_SRC == FROM_VEGHIST ||
        options.ALB_SRC == FROM_VEGHIST) {
        // for now forcing file is determined by the year
        sprintf(filenames.forcing[1], "%s%4d.nc", filenames.f_path_pfx[1],
                dmy[current].year);

        // global_param.forceoffset[1] resets every year since the met file restarts
        // every year
        if (current > 1 && (dmy[current].year != dmy[current - 1].year)) {
            global_param.forceoffset[1] = 0;
        }

        // only the time slice changes for the met file reads. The rest is constant
        d4start[2] = 0;
        d4start[3] = 0;
        d4count[0] = 1;
        d4count[1] = 1;
        d4count[2] = global_domain.n_ny;
        d4count[3] = global_domain.n_nx;

        // Leaf Area Index: lai
        if (options.LAI_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(filenames.forcing[1], "lai",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].LAI[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Partial veg cover fraction: fcov
        if (options.FCAN_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(filenames.forcing[1], "fcov",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].fcanopy[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Albedo: alb
        if (options.ALB_SRC == FROM_VEGHIST) {
            for (j = 0; j < NF; j++) {
                d4start[0] = global_param.forceskip[1] +
                             global_param.forceoffset[1] + j;
                for (v = 0; v < options.NVEGTYPES; v++) {
                    d4start[1] = v;
                    get_scatter_nc_field_double(filenames.forcing[1], "alb",
                                                d4start, d4count, dvar);
                    for (i = 0; i < local_domain.ncells_active; i++) {
                        vidx = veg_con_map[i].vidx[v];
                        if (vidx != NODATA_VEG) {
                            veg_hist[i][vidx].albedo[j] = (double) dvar[i];
                        }
                    }
                }
            }
        }

        // Update the offset counter
        global_param.forceoffset[1] += NF;
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
            // temperature in Celsius
            atmos[i].air_temp[j] -= CONST_TKFRZ;
            // precipitation in mm/period
            atmos[i].prec[j] *= global_param.snow_dt;
            // pressure in Pa
            // vapor pressure in Pa (we read specific humidity in kg/kg)
            atmos[i].vp[j] = q_to_vp(atmos[i].vp[j], atmos[i].pressure[j]);
            // vapor pressure deficit in Pa
            atmos[i].vpd[j] = svp(atmos[i].air_temp[j]) - atmos[i].vp[j];
            // air density in kg/m3
            atmos[i].density[j] = air_density(atmos[i].air_temp[j],
                                              atmos[i].pressure[j]);
            // snow flag
            atmos[i].snowflag[j] = will_it_snow(&(atmos[i].air_temp[j]),
                                                t_offset,
                                                param.SNOW_MAX_SNOW_TEMP,
                                                &(atmos[i].prec[j]), 1);
        }
        // Check on fcanopy
        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                for (j = 0; j < NF; j++) {
                    if ((veg_hist[i][vidx].fcanopy[j] < MIN_FCANOPY) &&
                        ((current == 0) ||
                         (options.FCAN_SRC == FROM_VEGHIST))) {
                        // Only issue this warning once if not using veg hist fractions
                        log_warn(
                            "cell %zu, veg` %d substep %zu fcanopy %f < minimum of %f; setting = %f\n", i, vidx, j,
                            veg_hist[i][vidx].fcanopy[j], MIN_FCANOPY,
                            MIN_FCANOPY);
                        veg_hist[i][vidx].fcanopy[j] = MIN_FCANOPY;
                    }
                }
            }
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

        for (v = 0; v < options.NVEGTYPES; v++) {
            vidx = veg_con_map[i].vidx[v];
            if (vidx != NODATA_VEG) {
                // not the correct way to calculate average albedo in general,
                // but leave for now (it's correct if albedo is constant over
                // the model step)
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

        // Optional inputs
        if (options.LAKES) {
            atmos[i].channel_in[NR] = average(atmos[i].channel_in, NF) * NF;
        }
        if (options.CARBON) {
            atmos[i].Catm[NR] = average(atmos[i].Catm, NF);
            atmos[i].fdir[NR] = average(atmos[i].fdir, NF);
            atmos[i].par[NR] = average(atmos[i].par, NF);
            // for coszen, use value at noon
            atmos[i].coszen[NR] = compute_coszen(
                local_domain.locations[i].latitude,
                local_domain.locations[i].longitude, soil_con[i].time_zone_lng,
                dmy[current].day_in_year, SEC_PER_DAY / 2);
        }
    }


    // cleanup
    free(dvar);
}

/******************************************************************************
 * @brief    Determine timestep and start year, month, day, and seconds of forcing files
 *****************************************************************************/
void
get_forcing_file_info(param_set_struct *param_set,
                      size_t            file_num)
{
    extern global_param_struct global_param;
    extern filenames_struct    filenames;

    double                     nc_times[2];
    double                     nc_time_origin;
    size_t                     start = 0;
    size_t                     count = 2;
    char                      *nc_unit_chars = NULL;
    char                      *calendar_char = NULL;
    unsigned short int         time_units;
    unsigned short int         calendar;
    dmy_struct                 nc_origin_dmy;
    dmy_struct                 nc_start_dmy;

    // read time info from netcdf file
    get_nc_field_double(filenames.forcing[0], "time", &start, &count, nc_times);
    get_nc_var_attr(filenames.forcing[0], "time", "units", &nc_unit_chars);
    get_nc_var_attr(filenames.forcing[0], "time", "calendar", &calendar_char);

    // parse the calendar string and check to make sure it matches the global clock
    calendar = calendar_from_chars(calendar_char);

    // parse the time units
    parse_nc_time_units(nc_unit_chars, &time_units, &nc_origin_dmy);

    // Get date/time of the first entry in the forcing file.
    nc_time_origin =
        date2num(0., &nc_origin_dmy, 0., calendar, TIME_UNITS_DAYS);
    num2date(nc_time_origin, nc_times[0], 0., calendar, time_units,
             &nc_start_dmy);

    // Assign file start date/time
    global_param.forceyear[file_num] = nc_start_dmy.year;
    global_param.forcemonth[file_num] = nc_start_dmy.month;
    global_param.forceday[file_num] = nc_start_dmy.day;
    global_param.forcesec[file_num] = nc_start_dmy.dayseconds;

    // calculate timestep in forcing file
    if (time_units == TIME_UNITS_DAYS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(1. / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_HOURS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(HOURS_PER_DAY / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_MINUTES) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(MIN_PER_DAY / (nc_times[1] - nc_times[0]));
    }
    else if (time_units == TIME_UNITS_SECONDS) {
        param_set->force_steps_per_day[file_num] =
            (size_t) nearbyint(SEC_PER_DAY / (nc_times[1] - nc_times[0]));
    }

    // check that this forcing file will work
    if (param_set->force_steps_per_day[file_num] !=
        global_param.model_steps_per_day) {
        log_err("Forcing file timestep must match the model timestep.  "
                "Model timesteps per day is set to %zu and the forcing file "
                "timestep is set to %zu", global_param.model_steps_per_day,
                param_set->force_steps_per_day[file_num])
    }
    if (calendar != global_param.calendar) {
        log_err("Calendar in forcing file does not match the calendar of "
                "VIC's clock: %s", calendar_char);
    }

    // Free attribute character arrays
    free(nc_unit_chars);
    free(calendar_char);
}
