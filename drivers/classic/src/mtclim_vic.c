/******************************************************************************
 * @section DESCRIPTION
 *
 * VIC version of MTCLIM 4.3
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

#include <time.h>
#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>
#include <vic_physical_constants.h>
#include <mtclim.h>

/******************************************************************************
 * @brief    allocates space for input and output data arrays
 *****************************************************************************/
int
data_alloc(const control_struct *ctrl,
           data_struct          *data)
{
    int ok = 1;
    int ndays;

    ndays = ctrl->ndays;

    if (ok && ctrl->inyear &&
        !(data->year = (int*) malloc(ndays * sizeof(int)))) {
        log_err("Error allocating for year array");
        ok = 0;
    }
    if (ok && !(data->yday = (int*) malloc(ndays * sizeof(int)))) {
        log_err("Error allocating for yearday array");
        ok = 0;
    }
    if (ok && !(data->tmax = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for tmax array");
        ok = 0;
    }
    if (ok && !(data->tmin = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for tmin array");
        ok = 0;
    }
    if (ok && !(data->prcp = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for prcp array");
        ok = 0;
    }
    if (ok && ctrl->indewpt &&
        !(data->tdew = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for input humidity array");
        ok = 0;
    }
    if (ok && !(data->s_tmax = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site Tmax array");
        ok = 0;
    }
    if (ok && !(data->s_tmin = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site Tmin array");
        ok = 0;
    }
    if (ok && !(data->s_tday = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site Tday array");
        ok = 0;
    }
    if (ok && !(data->s_prcp = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site prcp array");
        ok = 0;
    }
    if (ok && !(data->s_hum = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site VPD array");
        ok = 0;
    }
    if (ok && !(data->s_srad = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site radiation array");
        ok = 0;
    }
    if (ok && !(data->s_dayl = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site daylength array");
        ok = 0;
    }
    if (ok && !(data->s_swe = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site snowpack array");
        ok = 0;
    }
    /* start vic_change */
    if (ok && !(data->s_fdir = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for direct fraction array");
        ok = 0;
    }
    if (ok && !(data->s_tskc = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site cloudiness array");
        ok = 0;
    }
    if (ok && !(data->s_ppratio = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site pet/prcp ratio array");
        ok = 0;
    }
    if (ok && !(data->s_tfmax = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site pet/prcp ratio array");
        ok = 0;
    }
    if (ok && !(data->s_ttmax = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for site pet/prcp ratio array");
        ok = 0;
    }
    /* end vic_change */
    return (!ok);
} /* end of data_alloc() */

/******************************************************************************
 * @brief    calculates daily air temperatures
 *****************************************************************************/
int
calc_tair(const control_struct   *ctrl,
          const parameter_struct *p,
          data_struct            *data)
{
    int                      ok = 1;
    int                      i, ndays;
    double                   dz;
    double                   tmean, tmax, tmin;

    extern parameters_struct param;

    ndays = ctrl->ndays;
    /* calculate elevation difference in kilometers */
    dz = (p->site_elev - p->base_elev) / M_PER_KM;

    /* apply lapse rate corrections to tmax and tmin */

    /* Since tmax lapse rate usually has a larger absolute value than tmin
       lapse rate, it is possible at high elevation sites for these corrections
       to result in tmin > tmax. Check for that occurrence and force
       tmin = corrected tmax - 0.5 deg C. */
    for (i = 0; i < ndays; i++) {
        /* lapse rate corrections */
        data->s_tmax[i] = tmax = data->tmax[i] + (dz * p->tmax_lr);
        data->s_tmin[i] = tmin = data->tmin[i] + (dz * p->tmin_lr);

        /* derived temperatures */
        tmean = (tmax + tmin) / 2.0;
        data->s_tday[i] = ((tmax - tmean) * param.MTCLIM_TDAYCOEF) + tmean;
    }

    return (!ok);
}

/******************************************************************************
 * @brief    calculates daily total precipitation
 *****************************************************************************/
int
calc_prcp(const control_struct   *ctrl,
          const parameter_struct *p,
          data_struct            *data)
{
    int    ok = 1;
    int    i, ndays;
    double ratio;

    ndays = ctrl->ndays;

    /* start vic_change */
    ratio = -1.;
    if (p->site_isoh < 1e-10 && p->base_isoh < 1e-10) {
        /* If base_isoh and site_isoh are both small, set the ratio to 1.
           This handles the case in which annual precip is 0, resulting in
           base_isoh and site_isoh being 0 and their ratio being undefined. */
        ratio = 1.;
    }
    else if (p->base_isoh == 0) {
        log_err("Error in calc_prcp(): base_isoh == 0 and site_isoh/base_isoh "
                "== NaN.");
    }
    else {
        ratio = p->site_isoh / p->base_isoh;
    }
    /* end vic_change */

    if (ok) {
        for (i = 0; i < ndays; i++) {
            data->s_prcp[i] = data->prcp[i] * ratio;
        }
    }

    return (!ok);
}

/******************************************************************************
 * @brief    estimates the accumulation and melt of snow for radiation
 *           algorithm corrections
 *****************************************************************************/
int
snowpack(const control_struct *ctrl,
         data_struct          *data)
{
    int                      ok = 1;
    int                      i, ndays, count;
    int                      start_yday, prev_yday;
    double                   snowpack, newsnow, snowmelt, sum;

    extern parameters_struct param;

    ndays = ctrl->ndays;

    /* first pass to initialize SWE array */
    snowpack = 0.0;
    for (i = 0; i < ndays; i++) {
        newsnow = 0.0;
        snowmelt = 0.0;
        if (data->s_tmin[i] <= param.MTCLIM_SNOW_TCRIT) {
            newsnow = data->s_prcp[i];
        }
        else {
            snowmelt = param.MTCLIM_SNOW_TRATE *
                       (data->s_tmin[i] - param.MTCLIM_SNOW_TCRIT);
        }
        snowpack += newsnow - snowmelt;
        if (snowpack < 0.0) {
            snowpack = 0.0;
        }
        data->s_swe[i] = snowpack;
    }

    /* use the first pass to set the initial snowpack conditions for the
       first day of data */
    start_yday = data->yday[0];
    if (start_yday == 1) {
        prev_yday = DAYS_PER_YEAR;
    }
    else {
        prev_yday = start_yday - 1;
    }
    count = 0;
    sum = 0.0;
    for (i = 1; i < ndays; i++) {
        if (data->yday[i] == start_yday || data->yday[i] == prev_yday) {
            count++;
            sum += data->s_swe[i];
        }
    }

    /* Proceed with correction if there are valid days to reinitialize
       the snowpack estiamtes. Otherwise use the first-pass estimate. */
    if (count) {
        snowpack = sum / (double)count;
        for (i = 0; i < ndays; i++) {
            newsnow = 0.0;
            snowmelt = 0.0;
            if (data->s_tmin[i] <= param.MTCLIM_SNOW_TCRIT) {
                newsnow = data->s_prcp[i];
            }
            else {
                snowmelt = param.MTCLIM_SNOW_TRATE *
                           (data->s_tmin[i] - param.MTCLIM_SNOW_TCRIT);
            }
            snowpack += newsnow - snowmelt;
            if (snowpack < 0.0) {
                snowpack = 0.0;
            }
            data->s_swe[i] = snowpack;
        }
    }

    return (!ok);
}

/******************************************************************************
 * @brief    iterative estimation of shortwave radiation and humidity
 * @note     Note: too many changes to maintain the start/end vic change
 *           comments
 *****************************************************************************/
int
calc_srad_humidity_iterative(const control_struct   *ctrl,
                             const parameter_struct *p,
                             data_struct            *data,
                             double                **tiny_radfract)
{
    extern option_struct     options;
    extern parameters_struct param;

    int                      ok = 1;
    int                      i, j, ndays;
    int                      start_yday, end_yday, isloop;
    int                      ami, yday;
    double                   ttmax0[DAYS_PER_LYEAR];
    double                   flat_potrad[DAYS_PER_LYEAR];
    double                   slope_potrad[DAYS_PER_LYEAR];
    double                   daylength[DAYS_PER_LYEAR];
    double                  *dtr, *sm_dtr;
    double                  *parray, *window, *t_fmax, *tdew;
    double                  *pet;
    double                   sum_prcp, ann_prcp, effann_prcp;
    double                   sum_pet, ann_pet;
    double                   tmax, tmin;
    double                   t1, t2;
    double                   pratio;
    double                   lat, coslat, sinlat, dt, h, dh;
    double                   cosslp, sinslp, cosasp, sinasp;
    double                   bsg1, bsg2, bsg3;
    double                   decl, cosdecl, sindecl, cosegeom, sinegeom, coshss,
                             hss;
    double                   sc, dir_beam_topa;
    double                   sum_flat_potrad, sum_slope_potrad, sum_trans;
    double                   cosh, sinh;
    double                   cza, cbsa, coszeh, coszwh;
    double                   dir_flat_topa, am;
    double                   b;
    double                   pvs, vpd;
    double                   trans1, trans2;
    double                   pa;
    double                   sky_prop;
    double                   avg_horizon, slope_excess;
    double                   horizon_scalar, slope_scalar;

    /* optical airmass by degrees */
    double                   optam[21] = {
        2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37, 4.72, 5.12, 5.60,
        6.18, 6.88, 7.77, 8.90, 10.39, 12.44, 15.36, 19.79, 26.96, 30.00
    };

    /* start vic_change */
    int                      tinystep;
    int                      tinystepspday;
    /* end vic_change */

    int                      iter;
    int                      max_iter;
    double                   rmse_tdew;
    double                   tol;
    double                  *pva;
    double                  *tdew_save;
    double                  *pva_save;

    /* number of simulation days */
    ndays = ctrl->ndays;

    /* local array memory allocation */
    /* allocate space for DTR and smoothed DTR arrays */
    if (!(dtr = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for DTR array");
        ok = 0;
    }
    if (!(sm_dtr = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for smoothed DTR array");
        ok = 0;
    }
    /* allocate space for effective annual precip array */
    if (!(parray = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for effective annual precip array");
        ok = 0;
    }
    /* allocate space for the prcp totaling array */
    if (!(window = (double*) malloc((ndays + 90) * sizeof(double)))) {
        log_err("Error allocating for prcp totaling array");
        ok = 0;
    }
    /* allocate space for t_fmax */
    if (!(t_fmax = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for p_tt_max array");
        ok = 0;
    }
    /* allocate space for Tdew array */
    if (!(tdew = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for Tdew array");
        ok = 0;
    }
    /* allocate space for pet array */
    if (!(pet = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for pet array");
        ok = 0;
    }
    /* allocate space for pva array */
    if (!(pva = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for pva array");
        ok = 0;
    }
    /* allocate space for tdew_save array */
    if (!(tdew_save = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for tdew_save array");
        ok = 0;
    }
    /* allocate space for pva_save array */
    if (!(pva_save = (double*) malloc(ndays * sizeof(double)))) {
        log_err("Error allocating for pva_save array");
        ok = 0;
    }

    /* calculate diurnal temperature range for transmittance calculations */
    for (i = 0; i < ndays; i++) {
        tmax = data->tmax[i];
        tmin = data->tmin[i];
        if (tmax < tmin) {
            tmax = tmin;
        }
        dtr[i] = tmax - tmin;
    }

    /* smooth dtr array: After Bristow and Campbell, 1984 */
    if (ndays >= 30) { /* use 30-day antecedent smoothing window */
        if (pulled_boxcar(dtr, sm_dtr, ndays, 30, 0)) {
            log_err("Error in boxcar smoothing, calc_srad_humidity()");
            ok = 0;
        }
    }
    else { /* smoothing window width = ndays */
        if (pulled_boxcar(dtr, sm_dtr, ndays, ndays, 0)) {
            log_err("Error in boxcar smoothing, calc_srad_humidity()");
            ok = 0;
        }
    }

    /* calculate the annual total precip */
    sum_prcp = 0.0;
    for (i = 0; i < ndays; i++) {
        sum_prcp += data->s_prcp[i];
    }
    ann_prcp = (sum_prcp / (double)ndays) * CONST_DDAYS_PER_YEAR;
    if (ann_prcp == 0.0) {
        ann_prcp = 1.0;
    }

    /* Generate the effective annual precip, based on a 3-month
       moving-window. Requires some special case handling for the
       beginning of the record and for short records. */

    /* check if there are at least 90 days in this input file, if not,
       use a simple total scaled to effective annual precip */
    if (ndays < 90) {
        sum_prcp = 0.0;
        for (i = 0; i < ndays; i++) {
            sum_prcp += data->s_prcp[i];
        }
        effann_prcp = (sum_prcp / (double)ndays) * CONST_DDAYS_PER_YEAR;

        /* if the effective annual precip for this period
           is less than 8 cm, set the effective annual precip to 8 cm
           to reflect an arid condition, while avoiding possible
           division-by-zero errors and very large ratios (PET/Pann) */
        if (effann_prcp < 8.0) {
            effann_prcp = 8.0;
        }
        for (i = 0; i < ndays; i++) {
            parray[i] = effann_prcp;
        }
    }
    else {
        /* Check if the yeardays at beginning and the end of this input file
           match up. If so, use parts of the three months at the end
           of the input file to generate effective annual precip for
           the first 3-months. Otherwise, duplicate the first 90 days
           of the record. */
        start_yday = data->yday[0];
        end_yday = data->yday[ndays - 1];
        if (start_yday != 1) {
            isloop = (end_yday == start_yday - 1) ? 1 : 0;
        }
        else {
            isloop =
                (end_yday == DAYS_PER_YEAR || end_yday ==
                 DAYS_PER_LYEAR) ? 1 : 0;
        }

        /* fill the first 90 days of window */
        for (i = 0; i < 90; i++) {
            if (isloop) {
                window[i] = data->s_prcp[ndays - 90 + i];
            }
            else {
                window[i] = data->s_prcp[i];
            }
        }
        /* fill the rest of the window array */
        for (i = 0; i < ndays; i++) {
            window[i + 90] = data->s_prcp[i];
        }

        /* for each day, calculate the effective annual precip from
           scaled 90-day total */
        for (i = 0; i < ndays; i++) {
            sum_prcp = 0.0;
            for (j = 0; j < 90; j++) {
                sum_prcp += window[i + j];
            }
            sum_prcp = (sum_prcp / 90.0) * CONST_DDAYS_PER_YEAR;

            /* if the effective annual precip for this 90-day period
               is less than 8 cm, set the effective annual precip to 8 cm
               to reflect an arid condition, while avoiding possible
               division-by-zero errors and very large ratios (PET/Pann) */
            parray[i] = (sum_prcp < 8.0) ? 8.0 : sum_prcp;
        }
    } /* end if ndays >= 90 */

    /*****************************************
    *                                       *
    * start of the main radiation algorithm *
    *                                       *
    *****************************************/

    /* before starting the iterative algorithm between humidity and
       radiation, calculate all the variables that don't depend on
       humidity so they only get done once. */

    /* STEP (1) calculate pressure ratio (site/reference) = f(elevation) */
    t1 = 1.0 - (-1 * param.LAPSE_RATE * p->site_elev) / CONST_TSTD;
    t2 = CONST_G / (-1 * param.LAPSE_RATE * CONST_RDAIR);
    pratio = pow(t1, t2);

    /* STEP (2) correct initial transmittance for elevation */
    trans1 = pow(param.MTCLIM_TBASE, pratio);

    /* STEP (3) build 366-day array of ttmax0, potential rad, and daylength */

    /* precalculate the transcendentals */
    lat = p->site_lat;
    /* check for (+/-) 90 degrees latitude, throws off daylength calc */
    lat *= RAD_PER_DEG;
    if (lat > 1.5707) {
        lat = 1.5707;
    }
    if (lat < -1.5707) {
        lat = -1.5707;
    }
    coslat = cos(lat);
    sinlat = sin(lat);
    cosslp = cos(p->site_slp * RAD_PER_DEG);
    sinslp = sin(p->site_slp * RAD_PER_DEG);
    cosasp = cos(p->site_asp * RAD_PER_DEG);
    sinasp = sin(p->site_asp * RAD_PER_DEG);
    /* cosine of zenith angle for east and west horizons */
    coszeh = cos(1.570796 - (p->site_ehoriz * RAD_PER_DEG));
    coszwh = cos(1.570796 - (p->site_whoriz * RAD_PER_DEG));

    /* sub-daily time and angular increment information */
    dt = param.MTCLIM_SRADDT;              /* set timestep */
    dh = dt / CONST_SECPERRAD;      /* calculate hour-angle step */
    /* start vic_change */
    tinystepspday = CONST_CDAY / param.MTCLIM_SRADDT;
    /* end vic_change */

    /* begin loop through yeardays */
    for (i = 0; i < DAYS_PER_YEAR; i++) {
        /* calculate cos and sin of declination */
        decl = CONST_MINDECL *
               cos(((double)i + CONST_DAYSOFF) * CONST_RADPERDAY);
        cosdecl = cos(decl);
        sindecl = sin(decl);

        /* do some precalculations for beam-slope geometry (bsg) */
        bsg1 = -sinslp * sinasp * cosdecl;
        bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl;
        bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl;

        /* calculate daylength as a function of lat and decl */
        cosegeom = coslat * cosdecl;
        sinegeom = sinlat * sindecl;
        coshss = -(sinegeom) / cosegeom;
        if (coshss < -1.0) {
            coshss = -1.0; /* 24-hr daylight */
        }
        if (coshss > 1.0) {
            coshss = 1.0; /* 0-hr daylight */
        }
        hss = acos(coshss);            /* hour angle at sunset (radians) */
        /* daylength (seconds) */
        daylength[i] = 2.0 * hss * CONST_SECPERRAD;

        /* start vic_change */
        if (daylength[i] > CONST_CDAY) {
            daylength[i] = CONST_CDAY;
        }
        /* end vic_change */

        /* solar constant as a function of yearday (W/m^2) */
        sc = param.MTCLIM_SOLAR_CONSTANT + 45.5 *
             sin((2.0 * CONST_PI * (double)i / CONST_DDAYS_PER_YEAR) + 1.7);

        /* extraterrestrial radiation perpendicular to beam, total over
           the timestep (J) */
        dir_beam_topa = sc * dt;

        sum_trans = 0.0;
        sum_flat_potrad = 0.0;
        sum_slope_potrad = 0.0;

        /* begin sub-daily hour-angle loop, from -hss to hss */
        for (h = -hss; h < hss; h += dh) {
            /* precalculate cos and sin of hour angle */
            cosh = cos(h);
            sinh = sin(h);

            /* calculate cosine of solar zenith angle */
            cza = cosegeom * cosh + sinegeom;

            /* calculate cosine of beam-slope angle */
            cbsa = sinh * bsg1 + cosh * bsg2 + bsg3;

            /* check if sun is above a flat horizon */
            if (cza > 0.0) {
                /* when sun is above the ideal (flat) horizon, do all the
                   flat-surface calculations to determine daily total
                   transmittance, and save flat-surface potential radiation
                   for later calculations of diffuse radiation */

                /* potential radiation for this time period, flat surface,
                   top of atmosphere */
                dir_flat_topa = dir_beam_topa * cza;

                /* determine optical air mass */
                am = 1.0 / (cza + 0.0000001);
                if (am > 2.9) {
                    ami = (int)(acos(cza) / RAD_PER_DEG) - 69;
                    if (ami < 0) {
                        ami = 0;
                    }
                    if (ami > 20) {
                        ami = 20;
                    }
                    am = optam[ami];
                }

                /* correct instantaneous transmittance for this optical
                   air mass */
                trans2 = pow(trans1, am);

                /* instantaneous transmittance is weighted by potential
                   radiation for flat surface at top of atmosphere to get
                   daily total transmittance */
                sum_trans += trans2 * dir_flat_topa;

                /* keep track of total potential radiation on a flat
                   surface for ideal horizons */
                sum_flat_potrad += dir_flat_topa;

                /* keep track of whether this time step contributes to
                   component 1 (direct on slope) */
                if ((h < 0.0 && cza > coszeh && cbsa > 0.0) ||
                    (h >= 0.0 && cza > coszwh && cbsa > 0.0)) {
                    /* sun between east and west horizons, and direct on
                       slope. this period contributes to component 1 */
                    sum_slope_potrad += dir_beam_topa * cbsa;
                }
            } /* end if sun above ideal horizon */
            else {
                dir_flat_topa = -1;
            }

            /* start vic_change */
            tinystep =
                (MONTHS_PER_YEAR * SEC_PER_HOUR + h *
                 CONST_SECPERRAD) / param.MTCLIM_SRADDT;
            if (tinystep < 0) {
                tinystep = 0;
            }
            if (tinystep > tinystepspday - 1) {
                tinystep = tinystepspday - 1;
            }
            if (dir_flat_topa > 0) {
                tiny_radfract[i][tinystep] = dir_flat_topa;
            }
            else {
                tiny_radfract[i][tinystep] = 0;
            }
            /* end vic_change */
        } /* end of sub-daily hour-angle loop */

        /* start vic_change */
        if (daylength[i] && sum_flat_potrad > 0) {
            for (j = 0; j < tinystepspday; j++) {
                tiny_radfract[i][j] /= sum_flat_potrad;
            }
        }
        /* end vic_change */

        /* calculate maximum daily total transmittance and daylight average
           flux density for a flat surface and the slope */
        if (daylength[i]) {
            ttmax0[i] = sum_trans / sum_flat_potrad;
            flat_potrad[i] = sum_flat_potrad / daylength[i];
            slope_potrad[i] = sum_slope_potrad / daylength[i];
        }
        else {
            ttmax0[i] = 0.0;
            flat_potrad[i] = 0.0;
            slope_potrad[i] = 0.0;
        }
    } /* end of i=365 days loop */

    /* force yearday 366 = yearday 365 */
    ttmax0[DAYS_PER_YEAR] = ttmax0[DAYS_PER_YEAR - 1];
    flat_potrad[DAYS_PER_YEAR] = flat_potrad[DAYS_PER_YEAR - 1];
    slope_potrad[DAYS_PER_YEAR] = slope_potrad[DAYS_PER_YEAR - 1];
    daylength[DAYS_PER_YEAR] = daylength[DAYS_PER_YEAR - 1];

    /* start vic_change */
    for (j = 0; j < tinystepspday; j++) {
        tiny_radfract[DAYS_PER_YEAR][j] = tiny_radfract[DAYS_PER_YEAR - 1][j];
    }
    /* end vic_change */

    /* STEP (4)  calculate the sky proportion for diffuse radiation */

    /* uses the product of spherical cap defined by average horizon angle
       and the great-circle truncation of a hemisphere. this factor does not
       vary by yearday. */
    avg_horizon = (p->site_ehoriz + p->site_whoriz) / 2.0;
    horizon_scalar = 1.0 - sin(avg_horizon * RAD_PER_DEG);
    if (p->site_slp > avg_horizon) {
        slope_excess = p->site_slp - avg_horizon;
    }
    else {
        slope_excess = 0.0;
    }
    if (2.0 * avg_horizon > 180.0) {
        slope_scalar = 0.0;
    }
    else {
        slope_scalar = 1.0 - (slope_excess / (180.0 - 2.0 * avg_horizon));
        if (slope_scalar < 0.0) {
            slope_scalar = 0.0;
        }
    }
    sky_prop = horizon_scalar * slope_scalar;

    /* b parameter, and t_fmax not varying with Tdew, so these can be
       calculated once, outside the iteration between radiation and humidity
       estimates. Requires storing t_fmax in an array. */
    for (i = 0; i < ndays; i++) {
        /* b parameter from 30-day average of DTR */
        b = param.MTCLIM_B0 + param.MTCLIM_B1 *
            exp(-param.MTCLIM_B2 * sm_dtr[i]);

        /* proportion of daily maximum transmittance */
        t_fmax[i] = 1.0 - 0.9 * exp(-b * pow(dtr[i], param.MTCLIM_C));

        /* correct for precipitation if this is a rain day */
        if (data->prcp[i] > param.MTCLIM_SW_PREC_THRESH) {
            t_fmax[i] *= param.MTCLIM_RAIN_SCALAR;
        }
        data->s_tfmax[i] = t_fmax[i];
    }

    /* Initial values of vapor pressure, etc */
    if (ctrl->indewpt) {
        /* Observed Tdew supplied */
        for (i = 0; i < ndays; i++) {
            tdew[i] = data->tdew[i];
        }
    }
    else {
        /* Estimate Tdew */
        for (i = 0; i < ndays; i++) {
            tdew[i] = data->s_tmin[i];
        }
    }
    if (ctrl->invp) {
        /* Observed vapor pressure supplied */
        for (i = 0; i < ndays; i++) {
            pva[i] = data->s_hum[i];
        }
    }
    else {
        /* convert dewpoint to vapor pressure */
        for (i = 0; i < ndays; i++) {
            /* start vic_change */
            /* pva[i] = 610.7 * exp(17.38 * tdew[i] / (239.0 + tdew[i])); */
            pva[i] = svp(tdew[i]);
            /* end vic_change */
        }
    }

    /* Other values needed for srad_humidity calculation */
    pa = atm_pres(p->site_elev);
    for (i = 0; i < ndays; i++) {
        yday = data->yday[i] - 1;
        data->s_dayl[i] = daylength[yday];
        tdew_save[i] = tdew[i];
        pva_save[i] = pva[i];
    }

    /* Initial estimates of solar radiation, cloud fraction, etc. */
    compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0,
                                  flat_potrad, slope_potrad, sky_prop,
                                  daylength, pet, parray, pa, dtr);

    /* estimate annual PET */
    sum_pet = 0.0;
    for (i = 0; i < ndays; i++) {
        sum_pet += pet[i];
    }
    ann_pet = (sum_pet / (double)ndays) * CONST_DDAYS_PER_YEAR;

    /* Reset humidity terms if no iteration desired */
    if (ctrl->indewpt || ctrl->invp ||
        (options.VP_ITER == VP_ITER_ANNUAL && ann_pet / ann_prcp >= 2.5)) {
        for (i = 0; i < ndays; i++) {
            tdew[i] = tdew_save[i];
            pva[i] = pva_save[i];
        }
    }

    /* Set up srad-humidity iterations */
    if (options.VP_ITER == VP_ITER_ALWAYS ||
        (options.VP_ITER == VP_ITER_ANNUAL &&
         ann_pet / ann_prcp >= 2.5) || options.VP_ITER == VP_ITER_CONVERGE) {
        if (options.VP_ITER == VP_ITER_CONVERGE) {
            max_iter = 100;
        }
        else {
            max_iter = 2;
        }
    }
    else {
        max_iter = 1;
    }

    /* srad-humidity iterations */
    tol = 0.01;
    iter = 1;
    rmse_tdew = tol + 1;
    while (rmse_tdew > tol && iter < max_iter) {
        for (i = 0; i < ndays; i++) {
            tdew_save[i] = tdew[i];
        }
        compute_srad_humidity_onetime(ndays, ctrl, data, tdew, pva, ttmax0,
                                      flat_potrad, slope_potrad, sky_prop,
                                      daylength, pet, parray, pa, dtr);
        rmse_tdew = 0;
        for (i = 0; i < ndays; i++) {
            rmse_tdew += (tdew[i] - tdew_save[i]) * (tdew[i] - tdew_save[i]);
        }
        rmse_tdew /= ndays;
        rmse_tdew = pow(rmse_tdew, 0.5);
        iter++;
    }

    /* save humidity in output data structure */
    if (ctrl->outhum) {
        if (!ctrl->invp) {
            for (i = 0; i < ndays; i++) {
                data->s_hum[i] = pva[i];
            }
        }
    }
    else {
        /* output humidity as vapor pressure deficit (Pa) */
        for (i = 0; i < ndays; i++) {
            /* calculate saturated VP at tday */
            /* start vic_change */
            pvs = svp(data->s_tday[i]);
            /* end vic_change */
            vpd = pvs - pva[i];
            if (vpd < 0.0) {
                vpd = 0.0;
            }
            data->s_hum[i] = vpd;
        }
    }

    /* free local array memory */
    free(dtr);
    free(sm_dtr);
    free(parray);
    free(window);
    free(t_fmax);
    free(tdew);
    free(pet);
    free(pva);
    free(tdew_save);
    free(pva_save);

    return (!ok);
}

/******************************************************************************
 * @brief    calculate solar radiation and humidity without iterating.
 *****************************************************************************/
void
compute_srad_humidity_onetime(int                   ndays,
                              const control_struct *ctrl,
                              data_struct          *data,
                              double               *tdew,
                              double               *pva,
                              double               *ttmax0,
                              double               *flat_potrad,
                              double               *slope_potrad,
                              double                sky_prop,
                              double               *daylength,
                              double               *pet,
                              double               *parray,
                              double                pa,
                              double               *dtr)
{
    extern option_struct     options;
    extern parameters_struct param;

    int                      i;
    int                      yday;
    double                   t_tmax;
    double                   t_final;
    double                   pdif;
    double                   pdir;
    double                   srad1;
    double                   srad2;
    double                   sc;
    double                   potrad;
    double                   tmink;
    double                   ratio;
    double                   ratio2;
    double                   ratio3;
    double                   tdewk;

    for (i = 0; i < ndays; i++) {
        yday = data->yday[i] - 1;

        /*** Compute SW radiation ***/

        t_tmax = ttmax0[yday] + param.MTCLIM_ABASE * pva[i];
        if (t_tmax < 0.0001) {
            t_tmax = 0.0001;              // this is mainly for the case of observed VP supplied, for which t_tmax sometimes ends up being negative (when potential radiation is low and VP is high)
        }
        data->s_ttmax[i] = t_tmax;

        /* final daily total transmittance */
        t_final = t_tmax * data->s_tfmax[i];

        /* estimate fraction of radiation that is diffuse, on an
           instantaneous basis, from relationship with daily total
           transmittance in Jones (Plants and Microclimate, 1992)
           Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
           Fig 6.14, p. 122. */
        pdif = -1.25 * t_final + 1.25;
        if (pdif > 1.0) {
            pdif = 1.0;
        }
        if (pdif < 0.0) {
            pdif = 0.0;
        }

        /* estimate fraction of radiation that is direct, on an
           instantaneous basis */
        pdir = 1.0 - pdif;

        /* the daily total radiation is estimated as the sum of the
           following two components:
           1. The direct radiation arriving during the part of
           the day when there is direct beam on the slope.
           2. The diffuse radiation arriving over the entire daylength
           (when sun is above ideal horizon).
         */

        /* component 1 */
        srad1 = slope_potrad[yday] * t_final * pdir;

        /* component 2 (diffuse) */

        /* includes the effect of surface albedo in raising the diffuse
           radiation for obstructed horizons */
        srad2 = flat_potrad[yday] * t_final * pdif *
                (sky_prop + param.MTCLIM_DIF_ALB * (1.0 - sky_prop));

        /* snow pack influence on radiation */
        if (options.MTCLIM_SWE_CORR && data->s_swe[i] > 0.0) {
            /* snow correction in J/m2/day */
            sc = (1.32 + 0.096 * data->s_swe[i]) * 1e6;
            /* convert to W/m2 and check for zero daylength */
            if (daylength[yday] > 0.0) {
                sc /= daylength[yday];
            }
            else {
                sc = 0.0;
            }
            /* set a maximum correction of 100 W/m2 */
            if (sc > 100.0) {
                sc = 100.0;
            }
        }
        else {
            sc = 0.0;
        }

        /* save daily radiation */
        /* save cloud transmittance when rad is an input */
        if (ctrl->insw) {
            potrad =
                (srad1 + srad2 + sc) * daylength[yday] / t_final / CONST_CDAY;
            if (potrad > 0 && data->s_srad[i] > 0 && daylength[yday] > 0) {
                data->s_tfmax[i] = (data->s_srad[i]) / (potrad * t_tmax); // both of these are 24hr mean rad. here
                if (data->s_tfmax[i] > 1.0) {
                    data->s_tfmax[i] = 1.0;
                }
            }
            else {
                data->s_tfmax[i] = 1.0;
            }
        }
        else {
            data->s_srad[i] = srad1 + srad2 + sc;
        }

        /* start vic_change */
        if (options.LW_CLOUD == LW_CLOUD_DEARDORFF) {
            data->s_tskc[i] = (1. - data->s_tfmax[i]);
        }
        else {
            data->s_tskc[i] = sqrt((1. - data->s_tfmax[i]) / 0.65);
        }
        data->s_fdir[i] = pdir;
        /* end vic_change */
    }

    /*** Compute PET using SW radiation estimate, and update Tdew, pva ***/
    for (i = 0; i < ndays; i++) {
        tmink = data->s_tmin[i] + CONST_TKFRZ;
        pet[i] =
            calc_pet(data->s_srad[i], data->s_tday[i], pa, data->s_dayl[i]);

        /* calculate ratio (PET/effann_prcp) and correct the dewpoint */
        ratio = pet[i] / parray[i];
        data->s_ppratio[i] = ratio * CONST_DDAYS_PER_YEAR;
        ratio2 = ratio * ratio;
        ratio3 = ratio2 * ratio;
        tdewk = tmink *
                (-0.127 + 1.121 * (1.003 - 1.444 * ratio + 12.312 * ratio2 -
                                   32.766 *
                                   ratio3) + 0.0006 * (dtr[i]));
        tdew[i] = tdewk - CONST_TKFRZ;

        /* start vic_change */
        pva[i] = svp(tdew[i]);
        /* end vic_change */
    }

    return;
}

/******************************************************************************
 * @brief    frees the memory previously allocated by data_alloc()
 *****************************************************************************/
int
data_free(const control_struct *ctrl,
          data_struct          *data)
{
    int ok = 1;
    if (ctrl->inyear) {
        free(data->year);
    }
    free(data->yday);
    free(data->tmax);
    free(data->tmin);
    free(data->prcp);
    if (ctrl->indewpt) {
        free(data->tdew);
    }
    free(data->s_tmax);
    free(data->s_tmin);
    free(data->s_tday);
    free(data->s_prcp);
    free(data->s_hum);
    free(data->s_srad);
    free(data->s_dayl);
    free(data->s_swe);
    /* start vic_change */
    free(data->s_fdir);
    free(data->s_tskc);
    free(data->s_ppratio);
    free(data->s_tfmax);
    free(data->s_ttmax);
    /* end vic_change */
    return (!ok);
}

/******************************************************************************
 * @brief    calculates the potential evapotranspiration for aridity
 *           corrections in calc_vpd(), according to Kimball et al., 1997
 *****************************************************************************/
double
calc_pet(double rad,
         double ta,
         double pa,
         double dayl)
{
    /* input parameters and units :
       double rad      (W/m2)  daylight average incident shortwave radiation
       double ta       (deg C) daylight average air temperature
       double pa       (Pa)    air pressure
       double dayl     (s)     daylength
     */

    double rnet;     /* (W m-2) absorbed shortwave radiation avail. for ET */
    double lhvap;    /* (J kg-1) latent heat of vaporization of water */
    double gamma;    /* (Pa K-1) psychrometer parameter */
    double dt = 0.2; /* offset for saturation vapor pressure calculation */
    double t1, t2;   /* (deg C) air temperatures */
    double pvs1, pvs2; /* (Pa)   saturated vapor pressures */
    double pet;      /* (kg m-2 day-1) potential evapotranspiration */
    double s;        /* (Pa K-1) slope of saturated vapor pressure curve */

    /* calculate absorbed radiation, assuming albedo = 0.2  and ground
       heat flux = 10% of absorbed radiation during daylight */
    rnet = rad * 0.72;

    /* calculate latent heat of vaporization as a function of ta */
    lhvap = 2.5023e6 - 2430.54 * ta;

    /* calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
       where:
       cp       (J/kg K)   specific heat of air
       epsilon  (unitless) ratio of molecular weights of water and air
     */
    gamma = CONST_CPDAIR * pa / (lhvap * CONST_EPS);

    /* estimate the slope of the saturation vapor pressure curve at ta */
    /* temperature offsets for slope estimate */
    t1 = ta + dt;
    t2 = ta - dt;

    /* calculate saturation vapor pressures at t1 and t2, using formula from
       Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
       Meteorol. Mag., 114:49-56.
     */
    /* start vic_change */
    pvs1 = svp(t1);
    pvs2 = svp(t2);
    /* end vic_change */

    /* calculate slope of pvs vs. T curve near ta */
    s = (pvs1 - pvs2) / (t1 - t2);

    /* calculate PET using Priestly-Taylor approximation, with coefficient
       set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day */
    pet = (1.26 * (s / (s + gamma)) * rnet * dayl) / lhvap;

    /* return a value in centimeters/day, because this value is used in a ratio
       to annual total precip, and precip units are centimeters */
    return (pet / 10.0);
}

/******************************************************************************
 * @brief    calculates the atmospheric pressure as a function of elevation
 *****************************************************************************/
double
atm_pres(double elev)
{
    /* daily atmospheric pressure (Pa) as a function of elevation (m) */

    /* From the discussion on atmospheric statics in:
       Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
       Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
       (p. 168)
     */

    double                   t1, t2;
    double                   pa;

    extern parameters_struct param;

    t1 = 1.0 - (-1 * param.LAPSE_RATE * elev) / CONST_TSTD;
    t2 = CONST_G / (-1 * param.LAPSE_RATE * CONST_RDAIR);
    pa = CONST_PSTD * pow(t1, t2);

    return(pa);
}

/******************************************************************************
 * @brief    calculates a moving average of antecedent values in an array,
 *           using either a ramped (w_flag=1) or a flat (w_flag=0) weighting
 *****************************************************************************/
int
pulled_boxcar(double *input,
              double *output,
              int     n,
              int     w,
              int     w_flag)
{
    int     ok = 1;
    int     i, j;
    double *wt;
    double  total, sum_wt;

    if (w > n) {
        log_err("Boxcar longer than array...Resize boxcar and try again");
        ok = 0;
    }

    if (ok && !(wt = (double*) malloc(w * sizeof(double)))) {
        log_err("Allocation error in boxcar()");
        ok = 0;
    }

    if (ok) {
        /* when w_flag != 0, use linear ramp to weight tails,
           otherwise use constant weight */
        sum_wt = 0.0;
        if (w_flag) {
            for (i = 0; i < w; i++) {
                wt[i] = (double)(i + 1);
                sum_wt += wt[i];
            }
        }
        else {
            for (i = 0; i < w; i++) {
                wt[i] = 1.0;
                sum_wt += wt[i];
            }
        }

        /* fill the output array, starting with the point where a full
           boxcar can be calculated */
        for (i = w - 1; i < n; i++) {
            total = 0.0;
            for (j = 0; j < w; j++) {
                total += input[i - w + j + 1] * wt[j];
            }
            output[i] = total / sum_wt;
        }

        /* fill the first w elements of the output array with the value from
           the first full boxcar */
        for (i = 0; i < w - 1; i++) {
            output[i] = output[w - 1];
        }

        free(wt);
    } /* end if ok */

    return (!ok);
}
