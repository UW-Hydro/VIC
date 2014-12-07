/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes atmospheric variables for both the model time step,
 * and the time step used by the snow algorithm (if different). Air temperature
 * is estimated using MTCLIM (see routine for reference), atmospheric moisture
 * is estimated using Kimball's algorithm (see routine for reference), and
 * radiation is estimated using Bras's algorithms (see routines for reference).
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
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Initialize atmospheric variables for both the model time step.
 *****************************************************************************/
void
initialize_atmos(atmos_data_struct    *atmos,
                 dmy_struct           *dmy,
                 FILE                **infile,
                 veg_lib_struct       *veg_lib,
                 veg_con_struct       *veg_con,
                 veg_hist_struct     **veg_hist,
                 soil_con_struct      *soil_con,
                 out_data_file_struct *out_data_files,
                 out_data_struct      *out_data)
{
    extern option_struct       options;
    extern param_set_struct    param_set;
    extern global_param_struct global_param;
    extern parameters_struct   param;
    extern size_t              NR, NF;


    size_t                     i;
    size_t                     j;
    size_t                     v;
    size_t                     band;
    unsigned short             day;
    unsigned short             hour;
    size_t                     rec;
    int                        idx;
    unsigned                   uidx;
    int                        k;
    int                       *tmaxhour;
    int                       *tminhour;
    double                     cell_area;
    double                     theta_l;
    double                     theta_s;
    double                     hour_offset;
    double                     phi;
    double                     elevation;
    double                     slope;
    double                     aspect;
    double                     ehoriz;
    double                     whoriz;
    double                     annual_prec;
    double                     avgJulyAirTemp;
    double                    *Tfactor;
    bool                      *AboveTreeLine;
    double                     min_Tfactor;
    double                    *hourlyrad;
    double                    *fdir;
    double                    *prec;
    double                    *tmax;
    double                    *tmin;
    double                    *tair;
    double                    *tskc;
    double                    *daily_vp;
    double                    *dailyrad;
    size_t                     Ndays;
    size_t                     stepspday;
    double                     sum, sum2;
    double                  ***veg_hist_data;
    double                  ***local_veg_hist_data;
    double                   **forcing_data;
    double                   **local_forcing_data;
    int                        type;
    double                     delta_t_minus;
    double                     delta_t_plus;
    int                        have_dewpt;
    int                        have_shortwave;
    int                        hour_offset_int;
    int                        tmp_starthour, tmp_endhour;
    short                      local_startyear;
    short                      local_startmonth;
    short                      local_startday;
    short                      local_starthour;
    unsigned short             day_in_year, year, month, days_in_month;
    size_t                     tmp_nrecs;
    size_t                     Ndays_local;
    dmy_struct                *dmy_local;
    dmy_struct                 dmy_tmp;
    int                        month_days[MONTHS_PER_YEAR] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
    int                        fstepspday;
    int                        tmp_int;
    double                     tmp_double;
    int                        save_prec_supplied;
    int                        save_wind_supplied;
    int                        save_vp_supplied;

    theta_l = (double)soil_con->time_zone_lng;
    theta_s = (double)soil_con->lng;
    hour_offset = (theta_l - theta_s) * HOURS_PER_DAY / 360;
    if (hour_offset < 0) {
        hour_offset_int = (int)(hour_offset - 0.5);
    }
    else {
        hour_offset_int = (int)(hour_offset + 0.5);
    }

    // hour_offset is now the distance from the center of local time zone
    hour_offset -= hour_offset_int;
    phi = soil_con->lat;
    elevation = soil_con->elevation;
    slope = soil_con->slope;
    aspect = soil_con->aspect;
    ehoriz = soil_con->ehoriz;
    whoriz = soil_con->whoriz;
    annual_prec = soil_con->annual_prec;
    cell_area = soil_con->cell_area;
    avgJulyAirTemp = soil_con->avgJulyAirTemp;
    Tfactor = soil_con->Tfactor;
    AboveTreeLine = soil_con->AboveTreeLine;
    save_prec_supplied = param_set.TYPE[PREC].SUPPLIED;
    save_wind_supplied = param_set.TYPE[WIND].SUPPLIED;
    save_vp_supplied = param_set.TYPE[VP].SUPPLIED;

    /* Check on minimum forcing requirements */
    if (!param_set.TYPE[PREC].SUPPLIED &&
        ((!param_set.TYPE[RAINF].SUPPLIED &&
          (!param_set.TYPE[LSRAINF].SUPPLIED ||
           !param_set.TYPE[CRAINF].SUPPLIED)) ||
         ((!param_set.TYPE[SNOWF].SUPPLIED &&
           (!param_set.TYPE[LSSNOWF].SUPPLIED ||
            !param_set.TYPE[CSNOWF].SUPPLIED))))) {
        log_err("Input meteorological forcing files must contain some form of "
                "precipitation (PREC, or { {RAINF or {LSRAINF and CRAINF}} "
                "and {SNOWF or {LSSNOWF and CSNOWF}} }); check input files.");
    }

    if (!(param_set.TYPE[TMAX].SUPPLIED &&
          param_set.FORCE_DT[param_set.TYPE[TMAX].SUPPLIED - 1] ==
          HOURS_PER_DAY &&
          param_set.TYPE[TMIN].SUPPLIED &&
          param_set.FORCE_DT[param_set.TYPE[TMIN].SUPPLIED - 1] ==
          HOURS_PER_DAY) &&
        !(param_set.TYPE[AIR_TEMP].SUPPLIED &&
          param_set.FORCE_DT[param_set.TYPE[AIR_TEMP].SUPPLIED - 1] <
          HOURS_PER_DAY)) {
        log_err("Input meteorological forcing files must contain either: a. "
                "Daily TMAX and TMIN (maximum and minimum air temperature) or "
                "b. sub-daily AIR_TEMP (air temperature); check input files.");
    }

    /* Assign N_ELEM for veg-dependent forcings */
    param_set.TYPE[LAI_IN].N_ELEM = veg_con[0].vegetat_type_num;
    param_set.TYPE[VEGCOVER].N_ELEM = veg_con[0].vegetat_type_num;
    param_set.TYPE[ALBEDO].N_ELEM = veg_con[0].vegetat_type_num;

    /* compute number of simulation days */
    tmp_starthour = 0;
    tmp_endhour = HOURS_PER_DAY - global_param.dt;
    tmp_nrecs = global_param.nrecs + global_param.starthour - tmp_starthour +
                tmp_endhour - dmy[global_param.nrecs - 1].hour;
    Ndays = (tmp_nrecs * global_param.dt) / HOURS_PER_DAY;

    /* compute number of full model time steps per day */
    stepspday = HOURS_PER_DAY / global_param.dt;

    /* Compute number of days for MTCLIM (in local time); for sub-daily, we
       must pad start and end with dummy records */
    Ndays_local = Ndays;
    if (hour_offset_int != 0) {
        Ndays_local = Ndays + 1;
    }

    local_starthour = global_param.starthour - hour_offset_int;
    local_startday = global_param.startday;
    local_startmonth = global_param.startmonth;
    local_startyear = global_param.startyear;
    if (local_starthour < 0) {
        local_starthour += HOURS_PER_DAY;
        local_startday--;
        if (local_startday < 1) {
            local_startmonth--;
            if (local_startmonth < 1) {
                local_startmonth = MONTHS_PER_YEAR;
                local_startyear--;
            }
            local_startday = month_days[local_startmonth - 1];
            if (local_startyear % 4 == 0 && local_startmonth == 2) {
                local_startday++;
            }
        }
    }

    /* compute local version of dmy array */
    dmy_local =
        (dmy_struct *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(dmy_struct));
    if (dmy_local == NULL) {
        log_err("Memory allocation failure in initialize_atmos()");
    }
    day_in_year = local_startday;
    for (month = 1; month < local_startmonth; month++) {
        days_in_month = month_days[month - 1];
        if (local_startyear % 4 == 0 && month == 2) {
            days_in_month++;
        }
        day_in_year += days_in_month;
    }
    year = local_startyear;
    month = local_startmonth;
    day = local_startday;
    rec = 0;
    hour = 0;
    while (rec < Ndays_local * HOURS_PER_DAY) {
        dmy_local[rec].day_in_year = day_in_year;
        dmy_local[rec].year = year;
        dmy_local[rec].month = month;
        dmy_local[rec].day = day;
        dmy_local[rec].hour = hour;
        rec++;
        hour++;
        if (hour == HOURS_PER_DAY) {
            hour = 0;
            day_in_year++;
            day++;
            days_in_month = month_days[month - 1];
            if (year % 4 == 0 && month == 2) {
                days_in_month++;
            }
            if (day > days_in_month) {
                day = 1;
                month++;
                if (month > MONTHS_PER_YEAR) {
                    day_in_year = 1;
                    month = 1;
                    year++;
                }
            }
        }
    }

    /* mtclim routine memory allocations */

    hourlyrad = (double *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(double));
    prec = (double *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(double));
    tair = (double *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(double));
    tmax = (double *) calloc(Ndays_local, sizeof(double));
    tmaxhour = (int *)    calloc(Ndays_local, sizeof(int));
    tmin = (double *) calloc(Ndays_local, sizeof(double));
    tminhour = (int *)    calloc(Ndays_local, sizeof(int));
    tskc = (double *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(double));
    daily_vp = (double *) calloc(Ndays_local, sizeof(double));
    dailyrad = (double *) calloc(Ndays_local, sizeof(double));
    fdir = (double *) calloc(Ndays_local * HOURS_PER_DAY, sizeof(double));

    if (hourlyrad == NULL || prec == NULL || tair == NULL || tmax == NULL ||
        tmaxhour == NULL || tmin == NULL || tminhour == NULL || tskc == NULL ||
        daily_vp == NULL || dailyrad == NULL || fdir == NULL) {
        log_err("Memory allocation failure in initialize_atmos()");
    }

    /*******************************
       read in meteorological data
    *******************************/

    forcing_data = read_forcing_data(infile, global_param, &veg_hist_data);

    log_info("Read meteorological forcing file");

    /*************************************************
       Pre-processing
    *************************************************/

    /*************************************************
       Convert units from ALMA to VIC standard, if necessary
    *************************************************/
    if (options.ALMA_INPUT) {
        for (type = 0; type < N_FORCING_TYPES; type++) {
            if (param_set.TYPE[type].SUPPLIED) {
                /* Convert moisture flux rates to accumulated moisture flux
                   per time step */
                if (type == PREC ||
                    type == RAINF ||
                    type == CRAINF ||
                    type == LSRAINF ||
                    type == SNOWF ||
                    type == CSNOWF ||
                    type == LSSNOWF ||
                    type == CHANNEL_IN
                    ) {
                    for (i = 0; i < (global_param.nrecs * NF); i++) {
                        forcing_data[type][i] *=
                            param_set.FORCE_DT[param_set.TYPE[type].SUPPLIED -
                                               1] * SEC_PER_HOUR;
                    }
                }
                /* Convert temperatures from K to C */
                else if (type == AIR_TEMP || type == TMIN || type == TMAX) {
                    for (i = 0; i < (global_param.nrecs * NF); i++) {
                        forcing_data[type][i] -= CONST_TKFRZ;
                    }
                }
            }
        }
    }
    else {
        for (type = 0; type < N_FORCING_TYPES; type++) {
            if (param_set.TYPE[type].SUPPLIED) {
                /* Convert pressures from kPa to Pa */
                if (type == PRESSURE || type == VP) {
                    for (i = 0; i < (global_param.nrecs * NF); i++) {
                        forcing_data[type][i] *= PA_PER_KPA;
                    }
                }
            }
        }
    }

    /*************************************************
       If provided, translate rainfall and snowfall
       into total precipitation
       NOTE: this overwrites any PREC data that was supplied
    *************************************************/

    if (param_set.TYPE[RAINF].SUPPLIED && param_set.TYPE[SNOWF].SUPPLIED) {
        /* rainfall and snowfall supplied */
        if (forcing_data[PREC] == NULL) {
            forcing_data[PREC] = (double *)calloc((global_param.nrecs * NF),
                                                  sizeof(double));
        }
        for (i = 0; i < (global_param.nrecs * NF); i++) {
            forcing_data[PREC][i] = forcing_data[RAINF][i] +
                                    forcing_data[SNOWF][i];
        }
        param_set.TYPE[PREC].SUPPLIED = param_set.TYPE[RAINF].SUPPLIED;
    }
    else if (param_set.TYPE[CRAINF].SUPPLIED &&
             param_set.TYPE[LSRAINF].SUPPLIED &&
             param_set.TYPE[CSNOWF].SUPPLIED &&
             param_set.TYPE[LSSNOWF].SUPPLIED) {
        /* convective and large-scale rainfall and snowfall supplied */
        if (forcing_data[PREC] == NULL) {
            forcing_data[PREC] = (double *)calloc((global_param.nrecs * NF),
                                                  sizeof(double));
        }
        for (i = 0; i < (global_param.nrecs * NF); i++) {
            forcing_data[PREC][i] = forcing_data[CRAINF][i] +
                                    forcing_data[LSRAINF][i] +
                                    forcing_data[CSNOWF][i] +
                                    forcing_data[LSSNOWF][i];
        }
        param_set.TYPE[PREC].SUPPLIED = param_set.TYPE[LSRAINF].SUPPLIED;
    }

    /*************************************************
       If provided, translate WIND_E and WIND_N into WIND
       NOTE: this overwrites any WIND data that was supplied
    *************************************************/

    if (param_set.TYPE[WIND_E].SUPPLIED && param_set.TYPE[WIND_N].SUPPLIED) {
        /* specific wind_e and wind_n supplied */
        if (forcing_data[WIND] == NULL) {
            forcing_data[WIND] = (double *)calloc((global_param.nrecs * NF),
                                                  sizeof(double));
        }
        for (i = 0; i < (global_param.nrecs * NF); i++) {
            forcing_data[WIND][i] = sqrt(forcing_data[WIND_E][i] *
                                         forcing_data[WIND_E][i] +
                                         forcing_data[WIND_N][i] *
                                         forcing_data[WIND_N][i]);
        }
        param_set.TYPE[WIND].SUPPLIED = param_set.TYPE[WIND_E].SUPPLIED;
    }

    /*************************************************
       Create new forcing arrays referenced to local time
       This will simplify subsequent data processing
    *************************************************/

    local_forcing_data = (double **) calloc(N_FORCING_TYPES, sizeof(double*));
    if (local_forcing_data == NULL) {
        log_err("Memory allocation failure in initialize_atmos() for "
                "local_forcing_data");
    }
    local_veg_hist_data =
        (double ***) calloc(N_FORCING_TYPES, sizeof(double**));
    if (local_veg_hist_data == NULL) {
        log_err("Memory allocation failure in initialize_atmos() for "
                "local_fveg_hist_data");
    }
    for (type = 0; type < N_FORCING_TYPES; type++) {
        // Allocate enough space for hourly data
        if (type != ALBEDO && type != LAI_IN && type != VEGCOVER) {
            if ((local_forcing_data[type] =
                     (double *)calloc(Ndays_local * HOURS_PER_DAY,
                                      sizeof(double))) == NULL) {
                log_err("Memory allocation failure in initialize_atmos()");
            }
        }
        else {
            if ((local_veg_hist_data[type] =
                     (double **)calloc(param_set.TYPE[type].N_ELEM,
                                       sizeof(double*))) == NULL) {
                log_err("Memory allocation failure in initialize_atmos()");
            }
            for (v = 0; v < param_set.TYPE[type].N_ELEM; v++) {
                if ((local_veg_hist_data[type][v] =
                         (double *)calloc(Ndays_local * HOURS_PER_DAY,
                                          sizeof(double))) == NULL) {
                    log_err("Memory allocation failure in initialize_atmos()");
                }
            }
        }
        if (param_set.TYPE[type].SUPPLIED) {
            if (param_set.FORCE_DT[param_set.TYPE[type].SUPPLIED - 1] ==
                HOURS_PER_DAY) {
                /* Daily forcings in non-local time will straddle local day
                   boundaries and need to be padded with an extra day at
                   start or end */
                for (uidx = 0; uidx < Ndays_local; uidx++) {
                    k = uidx;
                    if (hour_offset_int > 0) {
                        k--;            // W. Hemisphere, in GMT time
                    }
                    if (k < 0) {
                        k = 0; // W. Hemisphere, in GMT time; pad extra day in front
                    }
                    if (k >= (int)Ndays) {
                        k = (int)Ndays - 1; // E. Hemisphere, in GMT time; pad extra day at end
                    }
                    if (type != ALBEDO && type != LAI_IN && type != VEGCOVER) {
                        local_forcing_data[type][uidx] = forcing_data[type][k];
                    }
                    else {
                        for (v = 0; v < param_set.TYPE[type].N_ELEM; v++) {
                            local_veg_hist_data[type][v][uidx] =
                                veg_hist_data[type][v][k];
                        }
                    }
                }
            }
            else {
                /* Local sub-daily forcings will be hourly for coding convenience
                   Sub-daily forcings need to a) start at hour 0, local time
                   and b) draw from the correct element of the supplied
                   forcings (if the supplied forcings are not in local time) */
                fstepspday = HOURS_PER_DAY /
                             param_set.FORCE_DT[param_set.TYPE[type].SUPPLIED -
                                                1];
                for (uidx = 0; uidx < (Ndays_local * HOURS_PER_DAY); uidx++) {
                    k = (uidx - global_param.starthour + hour_offset_int) /
                        param_set.FORCE_DT[param_set.TYPE[type].SUPPLIED - 1];
                    if (k < 0) {
                        k += fstepspday;
                    }
                    if (k >= (int)(Ndays * fstepspday)) {
                        k -= fstepspday;
                    }
                    if (type == PREC ||
                        type == RAINF ||
                        type == CRAINF ||
                        type == LSRAINF ||
                        type == SNOWF ||
                        type == CSNOWF ||
                        type == LSSNOWF ||
                        type == CHANNEL_IN
                        ) {
                        /* Amounts per step need to be scaled to new step length */
                        local_forcing_data[type][uidx] = forcing_data[type][k] /
                                                         param_set.FORCE_DT[
                            param_set.TYPE[type].SUPPLIED - 1];
                    }
                    else {
                        /* All other forcings are assumed constant over hourly
                           substeps */
                        if (type != ALBEDO && type != LAI_IN && type !=
                            VEGCOVER) {
                            local_forcing_data[type][uidx] =
                                forcing_data[type][k];
                        }
                        else {
                            for (v = 0; v < param_set.TYPE[type].N_ELEM; v++) {
                                local_veg_hist_data[type][v][uidx] =
                                    veg_hist_data[type][v][k];
                            }
                        }
                    }
                }
            }
        }
    }

    /*************************************************
       Incoming Channel Flow
    *************************************************/

    if (param_set.TYPE[CHANNEL_IN].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[CHANNEL_IN].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily channel_in provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].channel_in[j] =
                        local_forcing_data[CHANNEL_IN][uidx] /
                        (double)(NF * stepspday);                                                  // divide evenly over the day
                    atmos[rec].channel_in[j] *= MM_PER_M / cell_area; // convert to mm over grid cell
                    sum += atmos[rec].channel_in[j];
                }
                if (NF > 1) {
                    atmos[rec].channel_in[NR] = sum;
                }
            }
        }
        else {
            /* sub-daily channel_in provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    atmos[rec].channel_in[i] = 0;
                    while (hour < rec * global_param.dt +
                           (i +
                            1) * options.SNOW_STEP + global_param.starthour -
                           hour_offset_int) {
                        idx = hour;
                        if (idx < 0) {
                            idx += HOURS_PER_DAY;
                        }
                        atmos[rec].channel_in[i] +=
                            local_forcing_data[CHANNEL_IN][idx];
                        hour++;
                    }
                    atmos[rec].channel_in[i] *= MM_PER_M / cell_area; // convert to mm over grid cell
                    sum += atmos[rec].channel_in[i];
                }
                if (NF > 1) {
                    atmos[rec].channel_in[NR] = sum;
                }
            }
        }
    }
    else {
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                atmos[rec].channel_in[i] = 0;
                sum += atmos[rec].channel_in[i];
            }
            if (NF > 1) {
                atmos[rec].channel_in[NR] = sum;
            }
        }
    }

    /*************************************************
       Precipitation
    *************************************************/

    if (param_set.FORCE_DT[param_set.TYPE[PREC].SUPPLIED - 1] ==
        HOURS_PER_DAY) {
        /* daily precipitation provided */
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (j = 0; j < NF; j++) {
                hour = rec * global_param.dt + j * options.SNOW_STEP +
                       global_param.starthour - hour_offset_int;
                if (global_param.starthour - hour_offset_int < 0) {
                    hour += HOURS_PER_DAY;
                }
                uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                atmos[rec].prec[j] = local_forcing_data[PREC][uidx] /
                                     (double)(NF * stepspday);  // divide evenly over the day
                sum += atmos[rec].prec[j];
            }
            if (NF > 1) {
                atmos[rec].prec[NR] = sum;
            }
        }
        for (day = 0; day < Ndays_local; day++) {
            prec[day] = local_forcing_data[PREC][day];
        }
    }
    else {
        /* sub-daily precipitation provided */
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                hour = rec * global_param.dt + i * options.SNOW_STEP +
                       global_param.starthour - hour_offset_int;
                if (global_param.starthour - hour_offset_int < 0) {
                    hour += HOURS_PER_DAY;
                }
                atmos[rec].prec[i] = 0;
                for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                    atmos[rec].prec[i] += local_forcing_data[PREC][uidx];
                }
                sum += atmos[rec].prec[i];
            }
            if (NF > 1) {
                atmos[rec].prec[NR] = sum;
            }
        }
        for (day = 0; day < Ndays_local; day++) {
            prec[day] = 0;
            for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                prec[day] +=
                    local_forcing_data[PREC][day * HOURS_PER_DAY + hour];
            }
        }
    }

    /*************************************************
       Wind Speed
    *************************************************/

    if (param_set.TYPE[WIND].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[WIND].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily wind provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].wind[j] = local_forcing_data[WIND][uidx]; // assume constant over the day
                    sum += atmos[rec].wind[j];
                }
                if (NF > 1) {
                    atmos[rec].wind[NR] = sum / (double)NF;
                }
                if (global_param.dt == HOURS_PER_DAY) {
                    if (atmos[rec].wind[j] < param.WIND_SPEED_MIN) {
                        atmos[rec].wind[j] = param.WIND_SPEED_MIN;
                    }
                }
            }
        }
        else {
            /* sub-daily wind provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].wind[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        if (local_forcing_data[WIND][uidx] <
                            param.WIND_SPEED_MIN) {
                            atmos[rec].wind[i] += param.WIND_SPEED_MIN;
                        }
                        else {
                            atmos[rec].wind[i] +=
                                local_forcing_data[WIND][uidx];
                        }
                    }
                    atmos[rec].wind[i] /= options.SNOW_STEP;
                    sum += atmos[rec].wind[i];
                }
                if (NF > 1) {
                    atmos[rec].wind[NR] = sum / (double)NF;
                }
            }
        }
    }
    else {
        /* no wind data provided, use default constant */
        for (rec = 0; rec < global_param.nrecs; rec++) {
            for (i = 0; i < NF; i++) {
                atmos[rec].wind[i] = param.WIND_SPEED_DEFAULT;
            }
            atmos[rec].wind[NR] = param.WIND_SPEED_DEFAULT;
        }
    }

    /*************************************************
       Air Temperature, part 1.
    *************************************************/

    /************************************************
       Set maximum daily air temperature if provided
    ************************************************/

    if (param_set.TYPE[TMAX].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[TMAX].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily tmax provided */
            for (day = 0; day < Ndays_local; day++) {
                tmax[day] = local_forcing_data[TMAX][day];
            }
        }
        else {
            /* sub-daily tmax provided */
            for (day = 0; day < Ndays_local; day++) {
                tmax[day] = local_forcing_data[TMAX][day * HOURS_PER_DAY];
            }
        }
    }

    /************************************************
       Set minimum daily air temperature if provided
    ************************************************/

    if (param_set.TYPE[TMIN].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[TMIN].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily tmin provided */
            for (day = 0; day < Ndays_local; day++) {
                tmin[day] = local_forcing_data[TMIN][day];
            }
        }
        else {
            /* sub-daily tmin provided */
            for (day = 0; day < Ndays_local; day++) {
                tmin[day] = local_forcing_data[TMIN][day * HOURS_PER_DAY];
            }
        }
    }

    /*************************************************
       Store sub-daily air temperature if provided
    *************************************************/

    if (param_set.TYPE[AIR_TEMP].SUPPLIED) {
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                hour = rec * global_param.dt + i * options.SNOW_STEP +
                       global_param.starthour - hour_offset_int;
                if (global_param.starthour - hour_offset_int < 0) {
                    hour += HOURS_PER_DAY;
                }
                atmos[rec].air_temp[i] = 0;
                for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                    atmos[rec].air_temp[i] +=
                        local_forcing_data[AIR_TEMP][uidx];
                }
                atmos[rec].air_temp[i] /= options.SNOW_STEP;
                sum += atmos[rec].air_temp[i];
            }
            if (NF > 1) {
                atmos[rec].air_temp[NR] = sum / (double)NF;
            }
        }
    }

    /******************************************************
       Determine Tmax and Tmin from sub-daily temperatures
    ******************************************************/

    if (!(param_set.TYPE[TMAX].SUPPLIED && param_set.TYPE[TMIN].SUPPLIED)) {
        for (day = 0; day < Ndays_local; day++) {
            tmax[day] = tmin[day] = MISSING;
            for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                if (hour >= 9 &&
                    (tmax[day] == MISSING ||
                     local_forcing_data[AIR_TEMP][hour] >
                     tmax[day])) {
                    tmax[day] = local_forcing_data[AIR_TEMP][hour];
                }
                if (hour < 12 &&
                    (tmin[day] == MISSING ||
                     local_forcing_data[AIR_TEMP][hour] <
                     tmin[day])) {
                    tmin[day] = local_forcing_data[AIR_TEMP][hour];
                }
            }
        }
    }


    /*************************************************
       Vapor Pressure, part 1.
    *************************************************/

    if (!param_set.TYPE[VP].SUPPLIED) {
        /*************************************************
           If provided, translate specific humidity and atm. pressure
           into vapor pressure
           NOTE: if atm. pressure wasn't supplied, we must handle
           specific humidity after call to MTCLIM
        *************************************************/

        if (param_set.TYPE[QAIR].SUPPLIED &&
            param_set.TYPE[PRESSURE].SUPPLIED) {
            /* specific humidity and atm. pressure supplied */
            if (param_set.FORCE_DT[param_set.TYPE[QAIR].SUPPLIED - 1] ==
                HOURS_PER_DAY) {
                for (day = 0; day < Ndays_local; day++) {
                    if (param_set.FORCE_DT[param_set.TYPE[PRESSURE].SUPPLIED -
                                           1] < HOURS_PER_DAY) {
                        tmp_double = 0;
                        for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                            tmp_double +=
                                local_forcing_data[PRESSURE][day *
                                                             HOURS_PER_DAY +
                                                             hour];
                        }
                        tmp_double /= HOURS_PER_DAY;
                    }
                    else {
                        tmp_double = local_forcing_data[PRESSURE][day];
                    }
                    local_forcing_data[VP][day] =
                        local_forcing_data[QAIR][day] * tmp_double / CONST_EPS;
                    daily_vp[day] = local_forcing_data[VP][day];
                }
            }
            else {
                for (day = 0; day < Ndays_local; day++) {
                    daily_vp[day] = 0;
                    for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                        if (param_set.FORCE_DT[param_set.TYPE[PRESSURE].SUPPLIED
                                               -
                                               1] == HOURS_PER_DAY) {
                            tmp_double = local_forcing_data[PRESSURE][day];
                        }
                        else {
                            tmp_double =
                                local_forcing_data[PRESSURE][day *
                                                             HOURS_PER_DAY +
                                                             hour];
                        }
                        local_forcing_data[VP][day * HOURS_PER_DAY +
                                               hour] =
                            local_forcing_data[QAIR][day * HOURS_PER_DAY +
                                                     hour] * tmp_double /
                            CONST_EPS;
                        daily_vp[day] +=
                            local_forcing_data[VP][day * HOURS_PER_DAY + hour];
                    }
                    daily_vp[day] /= HOURS_PER_DAY;
                }
            }
            param_set.TYPE[VP].SUPPLIED = param_set.TYPE[QAIR].SUPPLIED;
        }

        /*************************************************
           If provided, translate relative humidity and air temperature
           into vapor pressure
           NOTE: if air temperature wasn't supplied, we must handle
           relative humidity after call to MTCLIM
        *************************************************/

        else if (param_set.TYPE[REL_HUMID].SUPPLIED &&
                 param_set.TYPE[AIR_TEMP].SUPPLIED) {
            /* relative humidity and air temperature supplied */
            if (param_set.FORCE_DT[param_set.TYPE[REL_HUMID].SUPPLIED - 1] ==
                HOURS_PER_DAY) {
                for (day = 0; day < Ndays_local; day++) {
                    if (param_set.FORCE_DT[param_set.TYPE[AIR_TEMP].SUPPLIED -
                                           1] < HOURS_PER_DAY) {
                        tmp_double = 0;
                        for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                            tmp_double +=
                                svp(
                                    local_forcing_data[AIR_TEMP][day *
                                                                 HOURS_PER_DAY +
                                                                 hour]);
                        }
                        tmp_double /= HOURS_PER_DAY;
                    }
                    else {
                        tmp_double = svp(local_forcing_data[AIR_TEMP][day]);
                    }
                    local_forcing_data[VP][day] =
                        local_forcing_data[REL_HUMID][day] * tmp_double /
                        FRACT_TO_PERCENT;
                    daily_vp[day] = local_forcing_data[VP][day];
                }
            }
            else {
                for (day = 0; day < Ndays_local; day++) {
                    daily_vp[day] = 0;
                    for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                        if (param_set.FORCE_DT[param_set.TYPE[AIR_TEMP].SUPPLIED
                                               -
                                               1] == HOURS_PER_DAY) {
                            tmp_double = svp(local_forcing_data[AIR_TEMP][day]);
                        }
                        else {
                            tmp_double =
                                svp(
                                    local_forcing_data[AIR_TEMP][day *
                                                                 HOURS_PER_DAY +
                                                                 hour]);
                        }
                        local_forcing_data[VP][day * HOURS_PER_DAY +
                                               hour] =
                            local_forcing_data[REL_HUMID][day * HOURS_PER_DAY +
                                                          hour] * tmp_double /
                            FRACT_TO_PERCENT;
                        daily_vp[day] +=
                            local_forcing_data[VP][day * HOURS_PER_DAY + hour];
                    }
                    daily_vp[day] /= HOURS_PER_DAY;
                }
            }
            param_set.TYPE[VP].SUPPLIED = param_set.TYPE[REL_HUMID].SUPPLIED;
        }
    } // end if VP not supplied

    /*************************************************
       If vapor pressure supplied, transfer to appropriate arrays
    *************************************************/

    if (param_set.TYPE[VP].SUPPLIED) {
        have_dewpt = 2; // flag for MTCLIM

        if (param_set.FORCE_DT[param_set.TYPE[VP].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily vp provided */
            for (day = 0; day < Ndays_local; day++) {
                daily_vp[day] = local_forcing_data[VP][day];
            }
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].vp[j] = local_forcing_data[VP][uidx]; // assume constant over the day
                    sum += atmos[rec].vp[j];
                }
                if (NF > 1) {
                    atmos[rec].vp[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily vp provided */
            for (day = 0; day < Ndays_local; day++) {
                daily_vp[day] = 0;
                for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                    daily_vp[day] +=
                        local_forcing_data[VP][day * HOURS_PER_DAY + hour];
                }
                daily_vp[day] /= HOURS_PER_DAY;
            }
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].vp[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].vp[i] += local_forcing_data[VP][uidx];
                    }
                    atmos[rec].vp[i] /= options.SNOW_STEP;
                    sum += atmos[rec].vp[i];
                }
                if (NF > 1) {
                    atmos[rec].vp[NR] = sum / (double)NF;
                }
            }
        }
    }
    else {
        have_dewpt = 0;
    } // end if VP supplied


    /*************************************************
       Shortwave, part 1.
    *************************************************/

    if (param_set.TYPE[SHORTWAVE].SUPPLIED) {
        have_shortwave = 1; // flag for MTCLIM
        for (day = 0; day < Ndays_local; day++) {
            for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                if (param_set.FORCE_DT[param_set.TYPE[SHORTWAVE].SUPPLIED -
                                       1] == HOURS_PER_DAY) {
                    hourlyrad[day * HOURS_PER_DAY +
                              hour] = local_forcing_data[SHORTWAVE][day];
                }
                else {
                    hourlyrad[day * HOURS_PER_DAY +
                              hour] =
                        local_forcing_data[SHORTWAVE][day * HOURS_PER_DAY +
                                                      hour];
                }
            }
        }
    }
    else {
        have_shortwave = 0;
    }

    /**************************************************
       Use MTCLIM algorithms to estimate hourly shortwave,
       daily vapor pressure, and cloud radiation attenuation.

       Requires prec, tmax, and tmin.

       If we already have observations of shortwave and/or
       vp, MTCLIM will use them to compute the other variables
       more accurately.
    **************************************************/
    mtclim_wrapper(have_dewpt, have_shortwave, hour_offset, elevation, slope,
                   aspect, ehoriz, whoriz, annual_prec, phi, Ndays_local,
                   dmy_local, prec, tmax, tmin, tskc, daily_vp, hourlyrad,
                   fdir);

    /***********************************************************
       Shortwave, part 2.
       Transfer the hourly shortwave from MTCLIM to atmos array.
       This hourly shortwave is one of the following:
       a) exactly equal to the supplied shortwave, if supplied shortwave was hourly
       b) equal to the supplied shortwave when aggregated up to the DT of the supplied shortwave (with hourly variability estimated by MTCLIM)
       c) completely estimated by MTCLIM, if no shortwave was supplied as a forcing
    ***********************************************************/

    // Ignore MTCLIM estimates if sub-daily SW was supplied
    if (param_set.TYPE[SHORTWAVE].SUPPLIED &&
        param_set.FORCE_DT[param_set.TYPE[SHORTWAVE].SUPPLIED - 1] <
        HOURS_PER_DAY) {
        for (day = 0; day < Ndays_local; day++) {
            for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                hourlyrad[day * HOURS_PER_DAY +
                          hour] =
                    local_forcing_data[SHORTWAVE][day * HOURS_PER_DAY + hour];
            }
        }
    }
    // Transfer hourlyrad to atmos structure
    for (rec = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        for (i = 0; i < NF; i++) {
            hour = rec * global_param.dt + i * options.SNOW_STEP +
                   global_param.starthour - hour_offset_int;
            if (global_param.starthour - hour_offset_int < 0) {
                hour += HOURS_PER_DAY;
            }
            atmos[rec].shortwave[i] = 0;
            for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                atmos[rec].shortwave[i] += hourlyrad[uidx];
            }
            atmos[rec].shortwave[i] /= options.SNOW_STEP;
            sum += atmos[rec].shortwave[i];
        }
        if (NF > 1) {
            atmos[rec].shortwave[NR] = sum / (double)NF;
        }
    }

    /**************************************************************************
       Air Temperature, part 2.
    **************************************************************************/

    /**************************************************************************
       Calculate the hours at which the minimum and maximum temperatures occur
       (if sub-daily air_temp will be estimated) and/or at which daily vapor
       pressure will occur (if daily vapor pressure is estimated)
    **************************************************************************/
    set_max_min_hour(hourlyrad, Ndays_local, tmaxhour, tminhour);

    if (!param_set.TYPE[AIR_TEMP].SUPPLIED) {
        /**********************************************************************
           Calculate the subdaily and daily temperature based on tmax and tmin
        **********************************************************************/
        HourlyT(1, Ndays_local, tmaxhour, tmax, tminhour, tmin, tair);
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                hour = rec * global_param.dt + i * options.SNOW_STEP +
                       global_param.starthour - hour_offset_int;
                if (global_param.starthour - hour_offset_int < 0) {
                    hour += HOURS_PER_DAY;
                }
                atmos[rec].air_temp[i] = 0;
                for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                    atmos[rec].air_temp[i] += tair[uidx];
                }
                atmos[rec].air_temp[i] /= options.SNOW_STEP;
                sum += atmos[rec].air_temp[i];
            }
            if (NF > 1) {
                atmos[rec].air_temp[NR] = sum / (double)NF;
            }
        }
    }


    /**************************************************************************
       Atmospheric Pressure and Density
    **************************************************************************/

    /*************************************************
       Store atmospheric density if provided (kg/m^3)
    *************************************************/

    if (param_set.TYPE[DENSITY].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[DENSITY].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily density provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].density[j] = local_forcing_data[DENSITY][uidx]; // assume constant over the day
                    sum += atmos[rec].density[j];
                }
                if (NF > 1) {
                    atmos[rec].density[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily density provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].density[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].density[i] +=
                            local_forcing_data[DENSITY][uidx];
                    }
                    atmos[rec].density[i] /= options.SNOW_STEP;
                    sum += atmos[rec].density[i];
                }
                if (NF > 1) {
                    atmos[rec].density[NR] = sum / (double)NF;
                }
            }
        }
    }

    /**************************************
       Estimate Atmospheric Pressure (Pa)
    **************************************/

    if (!param_set.TYPE[PRESSURE].SUPPLIED) {
        if (!param_set.TYPE[DENSITY].SUPPLIED) {
            /* Estimate pressure */
            if (options.PLAPSE) {
                /* Assume average virtual temperature in air column
                   between ground and sea level = CONST_TKFRZ+atmos[rec].air_temp[NR] + 0.5*elevation*param.LAPSE_RATE */
                for (rec = 0; rec < global_param.nrecs; rec++) {
                    atmos[rec].pressure[NR] = CONST_PSTD *
                                              exp(-elevation * CONST_G /
                                                  (CONST_RDAIR *
                                                   (CONST_TKFRZ +
                                                    atmos[rec].air_temp[NR] +
                                                    0.5 * elevation *
                                                    param.LAPSE_RATE)));
                    for (i = 0; i < NF; i++) {
                        atmos[rec].pressure[i] = CONST_PSTD *
                                                 exp(-elevation * CONST_G /
                                                     (CONST_RDAIR *
                                                      (CONST_TKFRZ +
                                                       atmos[rec].air_temp[i] +
                                                       0.5 * elevation *
                                                       param.LAPSE_RATE)));
                    }
                }
            }
            else {
                /* set pressure to constant value */
                for (rec = 0; rec < global_param.nrecs; rec++) {
                    atmos[rec].pressure[NR] = 95500.;
                    for (i = 0; i < NF; i++) {
                        atmos[rec].pressure[i] = atmos[rec].pressure[NR];
                    }
                }
            }
        }
        else {
            /* use observed densities to estimate pressure */
            if (options.PLAPSE) {
                for (rec = 0; rec < global_param.nrecs; rec++) {
                    atmos[rec].pressure[NR] =
                        (CONST_TKFRZ +
                         atmos[rec].air_temp[NR]) * atmos[rec].density[NR] *
                        CONST_RDAIR;
                    for (i = 0; i < NF; i++) {
                        atmos[rec].pressure[i] =
                            (CONST_TKFRZ +
                             atmos[rec].air_temp[i]) * atmos[rec].density[i] *
                            CONST_RDAIR;
                    }
                }
            }
            else {
                for (rec = 0; rec < global_param.nrecs; rec++) {
                    atmos[rec].pressure[NR] =
                        (275.0 +
                         atmos[rec].air_temp[NR]) * atmos[rec].density[NR] /
                        0.003486;
                    for (i = 0; i < NF; i++) {
                        atmos[rec].pressure[i] =
                            (275.0 +
                             atmos[rec].air_temp[i]) * atmos[rec].density[i] /
                            0.003486;
                    }
                }
            }
        }
    }
    else {
        /* observed atmospheric pressure supplied */
        if (param_set.FORCE_DT[param_set.TYPE[PRESSURE].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily pressure provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].pressure[j] = local_forcing_data[PRESSURE][uidx]; // assume constant over the day
                    sum += atmos[rec].pressure[j];
                }
                if (NF > 1) {
                    atmos[rec].pressure[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily pressure provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].pressure[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].pressure[i] +=
                            local_forcing_data[PRESSURE][uidx];
                    }
                    atmos[rec].pressure[i] /= options.SNOW_STEP;
                    sum += atmos[rec].pressure[i];
                }
                if (NF > 1) {
                    atmos[rec].pressure[NR] = sum / (double)NF;
                }
            }
        }
    }

    /********************************************************
       Estimate Atmospheric Density if not provided (kg/m^3)
    ********************************************************/

    if (!param_set.TYPE[DENSITY].SUPPLIED) {
        /* use pressure to estimate density */
        if (options.PLAPSE) {
            for (rec = 0; rec < global_param.nrecs; rec++) {
                atmos[rec].density[NR] = atmos[rec].pressure[NR] /
                                         (CONST_RDAIR *
                                          (CONST_TKFRZ +
                                           atmos[rec].air_temp[NR]));
                for (i = 0; i < NF; i++) {
                    atmos[rec].density[i] = atmos[rec].pressure[i] /
                                            (CONST_RDAIR *
                                             (CONST_TKFRZ +
                                              atmos[rec].air_temp[i]));
                }
            }
        }
        else {
            for (rec = 0; rec < global_param.nrecs; rec++) {
                atmos[rec].density[NR] = 0.003486 * atmos[rec].pressure[NR] /
                                         (275.0 + atmos[rec].air_temp[NR]);
                for (i = 0; i < NF; i++) {
                    atmos[rec].density[i] = 0.003486 * atmos[rec].pressure[i] /
                                            (275.0 + atmos[rec].air_temp[i]);
                }
            }
        }
    }

    /**************************************************************************
       Vapor Pressure, part 2.
    **************************************************************************/

    if (!param_set.TYPE[VP].SUPPLIED) {
        /* handle cases of daily QAIR or RH supplied without pressure or temperature */

        if (param_set.TYPE[QAIR].SUPPLIED &&
            param_set.FORCE_DT[param_set.TYPE[QAIR].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /**************************************************************************
               If we arrive here, it means we couldn't use Qair earlier because
               atmospheric pressure wasn't available at that time.  Now it is
               available, so use Qair and pressure to estimate vp.
            **************************************************************************/
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    daily_vp[uidx] = local_forcing_data[QAIR][uidx] *
                                     atmos[rec].pressure[j] / CONST_EPS;
                }
            }
        } // end if QAIR supplied
        else if (param_set.TYPE[REL_HUMID].SUPPLIED &&
                 param_set.FORCE_DT[param_set.TYPE[REL_HUMID].SUPPLIED - 1] ==
                 HOURS_PER_DAY) {
            /**************************************************************************
               If we arrive here, it means we couldn't use RH earlier because
               air temperature wasn't available at that time.  Now it is
               available, so use RH and temperature to estimate vp.
            **************************************************************************/
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    daily_vp[uidx] = local_forcing_data[REL_HUMID][uidx] * svp(
                        atmos[rec].air_temp[j]) / FRACT_TO_PERCENT;
                }
            }
        } // end if REL_HUMID supplied
    } // end if VP not supplied

    if (!param_set.TYPE[VP].SUPPLIED ||
        param_set.FORCE_DT[param_set.TYPE[VP].SUPPLIED - 1] == HOURS_PER_DAY) {
        /**************************************************
           Either no observations of VP, QAIR, or REL_HUMID were supplied,
           in which case we will use MTCLIM's estimates of daily vapor pressure,
           or daily VP was supplied.
           Now, calculate subdaily vapor pressure
        **************************************************/

        if (options.VP_INTERP) {
            /* Linearly interpolate between daily VP values, assuming they
               occurred at time of tmin */

            for (day = 0; day < Ndays_local; day++) {
                if (day == 0 && Ndays_local == 1) {
                    delta_t_minus = HOURS_PER_DAY;
                    delta_t_plus = HOURS_PER_DAY;
                }
                else if (day == 0) {
                    delta_t_minus = HOURS_PER_DAY;
                    delta_t_plus =
                        tminhour[day + 1] + HOURS_PER_DAY - tminhour[day];
                }
                else if (day == Ndays_local - 1) {
                    delta_t_minus = tminhour[day] + HOURS_PER_DAY -
                                    tminhour[day - 1];
                    delta_t_plus = HOURS_PER_DAY;
                }
                else {
                    delta_t_minus = tminhour[day] + HOURS_PER_DAY -
                                    tminhour[day - 1];
                    delta_t_plus =
                        tminhour[day + 1] + HOURS_PER_DAY - tminhour[day];
                }
                for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                    if (hour < tminhour[day]) {
                        if (day > 0) {
                            local_forcing_data[VP][day * HOURS_PER_DAY +
                                                   hour] =
                                daily_vp[day -
                                         1] +
                                (daily_vp[day] -
                                 daily_vp[day -
                                          1]) *
                                (hour + HOURS_PER_DAY -
                                 tminhour[day - 1]) / delta_t_minus;
                        }
                        else {
                            local_forcing_data[VP][day * HOURS_PER_DAY +
                                                   hour] = daily_vp[day];
                        }
                    }
                    else {
                        if (day < Ndays_local - 1) {
                            local_forcing_data[VP][day * HOURS_PER_DAY +
                                                   hour] = daily_vp[day] +
                                                           (daily_vp[day +
                                                                     1] -
                                                            daily_vp[day]) *
                                                           (hour -
                                                            tminhour[day]) /
                                                           delta_t_plus;
                        }
                        else {
                            local_forcing_data[VP][day * HOURS_PER_DAY +
                                                   hour] = daily_vp[day];
                        }
                    }
                }
            }
        }
        else {
            /* Hold VP constant throughout day */

            for (day = 0; day < Ndays_local; day++) {
                for (hour = 0; hour < HOURS_PER_DAY; hour++) {
                    local_forcing_data[VP][day * HOURS_PER_DAY +
                                           hour] = daily_vp[day];
                }
            }
        }

        /* Transfer sub-daily VP to atmos array */
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                hour = rec * global_param.dt + i * options.SNOW_STEP +
                       global_param.starthour - hour_offset_int;
                if (global_param.starthour - hour_offset_int < 0) {
                    hour += HOURS_PER_DAY;
                }
                atmos[rec].vp[i] = 0;
                for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                    atmos[rec].vp[i] += local_forcing_data[VP][uidx];
                }
                atmos[rec].vp[i] /= options.SNOW_STEP;
                sum += atmos[rec].vp[i];
            }
            if (NF > 1) {
                atmos[rec].vp[NR] = sum / (double)NF;
            }
        }

        /**************************************************
           If sub-daily specific or relative humidity were supplied without
           pressure or temperature,
           overwrite the sub-daily VP from MTCLIM here.
        **************************************************/
        if (param_set.TYPE[QAIR].SUPPLIED &&
            param_set.FORCE_DT[param_set.TYPE[QAIR].SUPPLIED - 1] <
            HOURS_PER_DAY) {
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].vp[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].vp[i] += local_forcing_data[QAIR][uidx] *
                                            atmos[rec].pressure[j] / CONST_EPS;
                    }
                    atmos[rec].vp[i] /= options.SNOW_STEP;
                    sum += atmos[rec].vp[i];
                }
                if (NF > 1) {
                    atmos[rec].vp[NR] = sum / (double)NF;
                }
            }
        }
        else if (param_set.TYPE[REL_HUMID].SUPPLIED &&
                 param_set.FORCE_DT[param_set.TYPE[REL_HUMID].SUPPLIED - 1] <
                 HOURS_PER_DAY) {
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].vp[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].vp[i] +=
                            local_forcing_data[REL_HUMID][uidx] *
                            svp(atmos[rec].air_temp[j]) /
                            FRACT_TO_PERCENT;
                    }
                    atmos[rec].vp[i] /= options.SNOW_STEP;
                    sum += atmos[rec].vp[i];
                }
                if (NF > 1) {
                    atmos[rec].vp[NR] = sum / (double)NF;
                }
            }
        }
    } // end computation of sub-daily VP

    /*************************************************
       Vapor Pressure Deficit
    *************************************************/

    for (rec = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        sum2 = 0;
        for (i = 0; i < NF; i++) {
            atmos[rec].vpd[i] = svp(atmos[rec].air_temp[i]) - atmos[rec].vp[i];
            if (atmos[rec].vpd[i] < 0) {
                atmos[rec].vpd[i] = 0;
                atmos[rec].vp[i] = svp(atmos[rec].air_temp[i]);
            }
            sum += atmos[rec].vpd[i];
            sum2 += atmos[rec].vp[i];
        }
        if (param_set.TYPE[VP].SUPPLIED || options.VP_INTERP) { // ensure that vp[NR] and vpd[NR] are accurate averages of vp[i] and vpd[i]
            if (NF > 1) {
                atmos[rec].vpd[NR] = sum / (double)NF;
            }
            if (NF > 1) {
                atmos[rec].vp[NR] = sum2 / (double)NF;
            }
        }
        else { // do not recompute vp[NR]; vpd[NR] is computed relative to vp[NR] and air_temp[NR]
            atmos[rec].vpd[NR] =
                (svp(atmos[rec].air_temp[NR]) - atmos[rec].vp[NR]);
        }
    }

    /*************************************************
       Cloud Transmissivity (from MTCLIM)
    *************************************************/

    for (rec = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        for (j = 0; j < NF; j++) {
            hour = rec * global_param.dt + j * options.SNOW_STEP +
                   global_param.starthour - hour_offset_int;
            if (global_param.starthour - hour_offset_int < 0) {
                hour += HOURS_PER_DAY;
            }
            uidx = (unsigned)((double)hour / HOURS_PER_DAY);
            atmos[rec].tskc[j] = tskc[uidx]; // assume constant over the day
            sum += atmos[rec].tskc[j];
        }
        if (NF > 1) {
            atmos[rec].tskc[NR] = sum / (double)NF;
        }
    }

    /*************************************************
       Longwave
    *************************************************/

    if (!param_set.TYPE[LONGWAVE].SUPPLIED) {
        /** Incoming longwave radiation not supplied **/
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                calc_longwave(&(atmos[rec].longwave[i]), atmos[rec].tskc[i],
                              atmos[rec].air_temp[i], atmos[rec].vp[i]);
                sum += atmos[rec].longwave[i];
            }
            if (NF > 1) {
                atmos[rec].longwave[NR] = sum / (double)NF;
            }
        }
    }
    else {
        if (param_set.FORCE_DT[param_set.TYPE[LONGWAVE].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily incoming longwave radiation provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].longwave[j] = local_forcing_data[LONGWAVE][uidx]; // assume constant over the day
                    sum += atmos[rec].longwave[j];
                }
                if (NF > 1) {
                    atmos[rec].longwave[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily incoming longwave radiation provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].longwave[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].longwave[i] +=
                            local_forcing_data[LONGWAVE][uidx];
                    }
                    atmos[rec].longwave[i] /= options.SNOW_STEP;
                    sum += atmos[rec].longwave[i];
                }
                if (NF > 1) {
                    atmos[rec].longwave[NR] = sum / (double)NF;
                }
            }
        }
    }

    /****************************************************
       Albedo
    ****************************************************/

    /* First, assign default climatology */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
            for (j = 0; j < NF; j++) {
                veg_hist[rec][v].albedo[j] =
                    veg_lib[veg_con[v].veg_class].albedo[dmy[rec].month - 1];
            }
        }
    }

    if (param_set.TYPE[ALBEDO].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[ALBEDO].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily albedo provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (j = 0; j < NF; j++) {
                        hour = rec * global_param.dt + j * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        if (global_param.starthour - hour_offset_int < 0) {
                            hour += HOURS_PER_DAY;
                        }
                        uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                        if (local_veg_hist_data[ALBEDO][v][uidx] != NODATA_VH) {
                            veg_hist[rec][v].albedo[j] =
                                local_veg_hist_data[ALBEDO][v][uidx];            // assume constant over the day
                        }
                        sum += veg_hist[rec][v].albedo[j];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].albedo[NR] = sum / (double)NF;
                    }
                }
            }
        }
        else {
            /* sub-daily albedo provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (i = 0; i < NF; i++) {
                        hour = rec * global_param.dt + i * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        veg_hist[rec][v].albedo[i] = 0;
                        while (hour < rec * global_param.dt +
                               (i +
                                1) * options.SNOW_STEP +
                               global_param.starthour -
                               hour_offset_int) {
                            idx = hour;
                            if (idx < 0) {
                                idx += HOURS_PER_DAY;
                            }
                            if (local_veg_hist_data[ALBEDO][v][idx] !=
                                NODATA_VH) {
                                veg_hist[rec][v].albedo[i] =
                                    local_veg_hist_data[ALBEDO][v][idx];
                            }
                            hour++;
                        }
                        sum += veg_hist[rec][v].albedo[i];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].albedo[NR] = sum / (double)NF;
                    }
                }
            }
        }
    }

    /****************************************************
       Leaf Area Index (LAI)
    ****************************************************/

    /* First, assign default climatology */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
            for (j = 0; j < NF; j++) {
                veg_hist[rec][v].LAI[j] =
                    veg_lib[veg_con[v].veg_class].LAI[dmy[rec].month - 1];
            }
        }
    }

    if (param_set.TYPE[LAI_IN].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[LAI_IN].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily LAI provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (j = 0; j < NF; j++) {
                        hour = rec * global_param.dt + j * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        if (global_param.starthour - hour_offset_int < 0) {
                            hour += HOURS_PER_DAY;
                        }
                        uidx = (double)((double)hour / HOURS_PER_DAY);
                        if (local_veg_hist_data[LAI_IN][v][uidx] != NODATA_VH) {
                            veg_hist[rec][v].LAI[j] =
                                local_veg_hist_data[LAI_IN][v][uidx];         // assume constant over the day
                        }
                        sum += veg_hist[rec][v].LAI[j];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].LAI[NR] = sum / (double)NF;
                    }
                }
            }
        }
        else {
            /* sub-daily LAI provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (i = 0; i < NF; i++) {
                        hour = rec * global_param.dt + i * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        veg_hist[rec][v].LAI[i] = 0;
                        while (hour < rec * global_param.dt +
                               (i +
                                1) * options.SNOW_STEP +
                               global_param.starthour -
                               hour_offset_int) {
                            idx = hour;
                            if (idx < 0) {
                                idx += HOURS_PER_DAY;
                            }
                            if (local_veg_hist_data[LAI_IN][v][idx] !=
                                NODATA_VH) {
                                veg_hist[rec][v].LAI[i] =
                                    local_veg_hist_data[LAI_IN][v][idx];
                            }
                            hour++;
                        }
                        sum += veg_hist[rec][v].LAI[i];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].LAI[NR] = sum / (double)NF;
                    }
                }
            }
        }
    }

    /****************************************************
       Fractional Vegetation Cover
    ****************************************************/

    /* First, assign default climatology */
    for (rec = 0; rec < global_param.nrecs; rec++) {
        for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
            for (j = 0; j < NF; j++) {
                veg_hist[rec][v].vegcover[j] =
                    veg_lib[veg_con[v].veg_class].vegcover[dmy[rec].month - 1];
            }
        }
    }

    if (param_set.TYPE[VEGCOVER].SUPPLIED) {
        if (param_set.FORCE_DT[param_set.TYPE[VEGCOVER].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily vegcover provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (j = 0; j < NF; j++) {
                        hour = rec * global_param.dt + j * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        if (global_param.starthour - hour_offset_int < 0) {
                            hour += HOURS_PER_DAY;
                        }
                        uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                        if (local_veg_hist_data[VEGCOVER][v][uidx] !=
                            NODATA_VH) {
                            veg_hist[rec][v].vegcover[j] =
                                local_veg_hist_data[VEGCOVER][v][uidx];              // assume constant over the day
                            if (veg_hist[rec][v].vegcover[j] < MIN_VEGCOVER) {
                                veg_hist[rec][v].vegcover[j] = MIN_VEGCOVER;
                            }
                        }
                        sum += veg_hist[rec][v].vegcover[j];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].vegcover[NR] = sum / (double)NF;
                    }
                }
            }
        }
        else {
            /* sub-daily vegcover provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                for (v = 0; v < veg_con[0].vegetat_type_num; v++) {
                    sum = 0;
                    for (i = 0; i < NF; i++) {
                        hour = rec * global_param.dt + i * options.SNOW_STEP +
                               global_param.starthour - hour_offset_int;
                        veg_hist[rec][v].vegcover[i] = 0;
                        while (hour < rec * global_param.dt +
                               (i +
                                1) * options.SNOW_STEP +
                               global_param.starthour -
                               hour_offset_int) {
                            idx = hour;
                            if (idx < 0) {
                                idx += HOURS_PER_DAY;
                            }
                            if (local_veg_hist_data[VEGCOVER][v][idx] !=
                                NODATA_VH) {
                                veg_hist[rec][v].vegcover[i] =
                                    local_veg_hist_data[VEGCOVER][v][idx];
                                if (veg_hist[rec][v].vegcover[i] <
                                    MIN_VEGCOVER) {
                                    veg_hist[rec][v].vegcover[i] = MIN_VEGCOVER;
                                }
                            }
                            hour++;
                        }
                        sum += veg_hist[rec][v].vegcover[i];
                    }
                    if (NF > 1) {
                        veg_hist[rec][v].vegcover[NR] = sum / (double)NF;
                    }
                }
            }
        }
    }

    /*************************************************
       Cosine of Solar Zenith Angle
    *************************************************/

    for (rec = 0; rec < global_param.nrecs; rec++) {
        dmy_tmp.year = dmy[rec].year;
        dmy_tmp.month = dmy[rec].month;
        dmy_tmp.day = dmy[rec].day;
        dmy_tmp.day_in_year = dmy[rec].day_in_year;
        for (j = 0; j < NF; j++) {
            hour = rec * global_param.dt + j * options.SNOW_STEP +
                   global_param.starthour - hour_offset_int;
            if (global_param.starthour - hour_offset_int < 0) {
                hour += HOURS_PER_DAY;
            }
            dmy_tmp.hour = hour + 0.5 * options.SNOW_STEP;
            atmos[rec].coszen[j] =
                compute_coszen(phi, theta_s, theta_l, dmy_tmp);
        }
        if (NF > 1) {
            dmy_tmp.hour = dmy[rec].hour + 0.5 * global_param.dt;
            atmos[rec].coszen[NR] = compute_coszen(phi, theta_s, theta_l,
                                                   dmy_tmp);
        }
    }

    /*************************************************
       Direct Shortwave Fraction (from MTCLIM)
    *************************************************/

    for (rec = 0; rec < global_param.nrecs; rec++) {
        sum = 0;
        for (j = 0; j < NF; j++) {
            hour = rec * global_param.dt + j * options.SNOW_STEP +
                   global_param.starthour - hour_offset_int;
            if (global_param.starthour - hour_offset_int < 0) {
                hour += HOURS_PER_DAY;
            }
            uidx = (unsigned)((double)hour / HOURS_PER_DAY);
            atmos[rec].fdir[j] = fdir[uidx]; // assume constant over the day
            sum += atmos[rec].fdir[j];
        }
        if (NF > 1) {
            atmos[rec].fdir[NR] = sum / (double)NF;
        }
    }

    /*************************************************
       Photosynthetically Active Radiation
    *************************************************/

    if (!param_set.TYPE[PAR].SUPPLIED) {
        /** par not supplied **/
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                atmos[rec].par[i] = param.CARBON_SW2PAR *
                                    atmos[rec].shortwave[i];
                sum += atmos[rec].par[i];
            }
            if (NF > 1) {
                atmos[rec].par[NR] = sum / (double)NF;
            }
        }
    }
    else {
        if (param_set.FORCE_DT[param_set.TYPE[PAR].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily par provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    tmp_int = (int)(rec / stepspday) * stepspday;
                    tmp_double = 0;
                    for (i = 0; i < stepspday; i++) {
                        tmp_double += atmos[tmp_int + i].shortwave[NR];
                    }
                    tmp_double /= stepspday;
                    if (tmp_double > 0) {
                        atmos[rec].par[j] = local_forcing_data[PAR][uidx] *
                                            atmos[rec].shortwave[j] /
                                            tmp_double;
                    }
                    else {
                        atmos[rec].par[j] = 0;
                    }
                    sum += atmos[rec].par[j];
                }
                if (NF > 1) {
                    atmos[rec].par[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily par provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].par[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].par[i] += local_forcing_data[PAR][uidx];
                    }
                    atmos[rec].par[i] /= options.SNOW_STEP;
                    sum += atmos[rec].par[i];
                }
                if (NF > 1) {
                    atmos[rec].par[NR] = sum / (double)NF;
                }
            }
        }
    }

    /*************************************************
       Atmospheric Carbon Dioxide Mixing Ratio
    *************************************************/

    if (!param_set.TYPE[CATM].SUPPLIED) {
        /** Atmospheric carbon dioxide concentration not supplied **/
        for (rec = 0; rec < global_param.nrecs; rec++) {
            sum = 0;
            for (i = 0; i < NF; i++) {
                atmos[rec].Catm[i] = param.CARBON_CATMCURRENT * PPM_to_MIXRATIO; // convert ppm to mixing ratio
                sum += atmos[rec].Catm[i];
            }
            if (NF > 1) {
                atmos[rec].Catm[NR] = sum / (double)NF;
            }
        }
    }
    else {
        if (param_set.FORCE_DT[param_set.TYPE[CATM].SUPPLIED - 1] ==
            HOURS_PER_DAY) {
            /* daily atmospheric carbon dioxide concentration provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (j = 0; j < NF; j++) {
                    hour = rec * global_param.dt + j * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    uidx = (unsigned)((double)hour / HOURS_PER_DAY);
                    atmos[rec].Catm[j] = local_forcing_data[CATM][uidx] *
                                         PPM_to_MIXRATIO;                                  // convert ppm to mixing ratio
                    sum += atmos[rec].Catm[j];
                }
                if (NF > 1) {
                    atmos[rec].Catm[NR] = sum / (double)NF;
                }
            }
        }
        else {
            /* sub-daily atmospheric carbon dioxide concentration provided */
            for (rec = 0; rec < global_param.nrecs; rec++) {
                sum = 0;
                for (i = 0; i < NF; i++) {
                    hour = rec * global_param.dt + i * options.SNOW_STEP +
                           global_param.starthour - hour_offset_int;
                    if (global_param.starthour - hour_offset_int < 0) {
                        hour += HOURS_PER_DAY;
                    }
                    atmos[rec].Catm[i] = 0;
                    for (uidx = hour; uidx < hour + options.SNOW_STEP; uidx++) {
                        atmos[rec].Catm[i] += local_forcing_data[CATM][uidx] *
                                              PPM_to_MIXRATIO;                   // convert ppm to mixing ratio
                    }
                    atmos[rec].Catm[i] /= options.SNOW_STEP;
                    sum += atmos[rec].Catm[i];
                }
                if (NF > 1) {
                    atmos[rec].Catm[NR] = sum / (double)NF;
                }
            }
        }
    }

    /****************************************************
       Determine if Snow will Fall During Each Time Step
    ****************************************************/

    if (!options.OUTPUT_FORCE) {
        min_Tfactor = Tfactor[0];
        for (band = 1; band < options.SNOW_BAND; band++) {
            if (Tfactor[band] < min_Tfactor) {
                min_Tfactor = Tfactor[band];
            }
        }
        for (rec = 0; rec < global_param.nrecs; rec++) {
            atmos[rec].snowflag[NR] = false;
            for (i = 0; i < NF; i++) {
                if ((atmos[rec].air_temp[i] + min_Tfactor) <
                    param.SNOW_MAX_SNOW_TEMP &&
                    atmos[rec].prec[i] > 0) {
                    atmos[rec].snowflag[i] = true;
                    atmos[rec].snowflag[NR] = true;
                }
                else {
                    atmos[rec].snowflag[i] = false;
                }
            }
        }
    }

    param_set.TYPE[PREC].SUPPLIED = save_prec_supplied;
    param_set.TYPE[WIND].SUPPLIED = save_wind_supplied;
    param_set.TYPE[VP].SUPPLIED = save_vp_supplied;

    // Free temporary parameters
    free(hourlyrad);
    free(prec);
    free(tair);
    free(tmax);
    free(tmaxhour);
    free(tmin);
    free(tminhour);
    free(tskc);
    free(daily_vp);
    free(dailyrad);
    free(fdir);

    for (i = 0; i < N_FORCING_TYPES; i++) {
        if (param_set.TYPE[i].SUPPLIED) {
            if (i != ALBEDO && i != LAI_IN && i != VEGCOVER) {
                free(forcing_data[i]);
            }
            else {
                for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                    free(veg_hist_data[i][j]);
                }
                free(veg_hist_data[i]);
            }
        }
        if (i != ALBEDO && i != LAI_IN && i != VEGCOVER) {
            free(local_forcing_data[i]);
        }
        else {
            for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                free(local_veg_hist_data[i][j]);
            }
            free(local_veg_hist_data[i]);
        }
    }
    free(forcing_data);
    free(local_forcing_data);
    free(veg_hist_data);
    free(local_veg_hist_data);
    free((char *)dmy_local);

    if (!options.OUTPUT_FORCE) {
        // If COMPUTE_TREELINE is TRUE and the treeline computation hasn't
        // specifically been turned off for this cell (by supplying avgJulyAirTemp
        // and setting it to -999), calculate which snowbands are above the
        // treeline, based on average July air temperature.
        if (options.COMPUTE_TREELINE) {
            if (!(options.JULY_TAVG_SUPPLIED && avgJulyAirTemp == -999)) {
                compute_treeline(atmos, dmy, avgJulyAirTemp, Tfactor,
                                 AboveTreeLine);
            }
        }
    }
    else {
        // If OUTPUT_FORCE is TRUE then the full
        // forcing data array is dumped into a new set of files.
        write_forcing_file(atmos, global_param.nrecs, out_data_files, out_data);
    }
}
