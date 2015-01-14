/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize and call mtclim routines to estimate meteorological variables
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

/******************************************************************************/
/*			    PREPROCESSOR DIRECTIVES                           */
/******************************************************************************/

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>
#include <mtclim.h>

/******************************************************************************/
/*			TYPE DEFINITIONS, GLOBALS, ETC.                       */
/******************************************************************************/

/******************************************************************************/
/*			      FUNCTION PROTOTYPES                             */
/******************************************************************************/
void mtclim_init(int have_dewpt, int have_shortwave, double elevation,
                 double slope, double aspect, double ehoriz, double whoriz,
                 double annual_prcp, double lat, size_t Ndays, dmy_struct *dmy,
                 double *prec, double *tmax, double *tmin, double *vp,
                 double *subdailyrad, double **tiny_radfract,
                 control_struct *ctrl, parameter_struct *p,
                 data_struct *mtclim_data);

void mtclim_to_vic(double sec_offset, dmy_struct *dmy, double **tiny_radfract,
                   control_struct *ctrl, data_struct *mtclim_data, double *tskc,
                   double *vp, double *subdailyrad, double *fdir);

/******************************************************************************
 * @brief    interface between VIC and MTCLIM.
 *****************************************************************************/
void
mtclim_wrapper(int         have_dewpt,
               int         have_shortwave,
               double      sec_offset,
               double      elevation,
               double      slope,
               double      aspect,
               double      ehoriz,
               double      whoriz,
               double      annual_prcp,
               double      lat,
               size_t      Ndays,
               dmy_struct *dmy,
               double     *prec,
               double     *tmax,
               double     *tmin,
               double     *tskc,
               double     *vp,
               double     *subdailyrad,
               double     *fdir)
{
    control_struct   ctrl;
    parameter_struct p;
    data_struct      mtclim_data;
    double         **tiny_radfract;
    int              i;

    /* allocate space for the tiny_radfract array */
    tiny_radfract = (double **) calloc(DAYS_PER_LYEAR, sizeof(double*));
    if (tiny_radfract == NULL) {
        log_err("Memory allocation error in mtclim_init() ...");
    }
    for (i = 0; i < DAYS_PER_LYEAR; i++) {
        tiny_radfract[i] = (double *) calloc((CONST_CDAY), sizeof(double));
        if (tiny_radfract[i] == NULL) {
            log_err("Memory allocation error in mtclim_init() ...");
        }
    }

    /* initialize the mtclim data structures */
    mtclim_init(have_dewpt, have_shortwave, elevation, slope, aspect, ehoriz,
                whoriz,
                annual_prcp, lat, Ndays, dmy, prec,
                tmax, tmin, vp, subdailyrad, tiny_radfract, &ctrl, &p,
                &mtclim_data);

    /* calculate daily air temperatures */
    if (calc_tair(&ctrl, &p, &mtclim_data)) {
        log_err("Error in calc_tair()... exiting");
    }

    /* calculate daily precipitation */
    if (calc_prcp(&ctrl, &p, &mtclim_data)) {
        log_err("Error in calc_prcp()... exiting");
    }

    /* calculate daily snowpack using simple model (this is only for radiation
       correction, *not* the same as the VIC snowpack estimate) */
    if (snowpack(&ctrl, &mtclim_data)) {
        log_err("Error in snowpack()... exiting");
    }

    /* calculate srad and humidity with iterative algorithm */
    if (calc_srad_humidity_iterative(&ctrl, &p, &mtclim_data, tiny_radfract)) {
        log_err("Error in calc_srad_humidity_iterative()... exiting");
    }

    /* translate the mtclim structures back to the VIC data structures */
    mtclim_to_vic(sec_offset, dmy, tiny_radfract, &ctrl, &mtclim_data, tskc,
                  vp, subdailyrad, fdir);

    /* clean up */
    if (data_free(&ctrl, &mtclim_data)) {
        log_err("Error in data_free()... exiting");
    }
    for (i = 0; i < DAYS_PER_LYEAR; i++) {
        free(tiny_radfract[i]);
    }
    free(tiny_radfract);
}

/******************************************************************************
 * @brief    Initialize mtclim control structures
 *****************************************************************************/
void
mtclim_init(int               have_dewpt,
            int               have_shortwave,
            double            elevation,
            double            slope,
            double            aspect,
            double            ehoriz,
            double            whoriz,
            double            annual_prcp,
            double            lat,
            size_t            Ndays,
            dmy_struct       *dmy,
            double           *prec,
            double           *tmax,
            double           *tmin,
            double           *vp,
            double           *subdailyrad,
            double          **tiny_radfract,
            control_struct   *ctrl,
            parameter_struct *p,
            data_struct      *mtclim_data)
{
    extern parameters_struct   param;
    extern global_param_struct global_param;

    size_t                     i, j;
    size_t                     tinystepspday;
    size_t                     atmos_steps_per_day;

    atmos_steps_per_day = global_param.model_steps_per_day;
    if (atmos_steps_per_day < (HOURS_PER_DAY)) {
        atmos_steps_per_day = (HOURS_PER_DAY);
    }

    /* initialize the control structure */

    ctrl->ndays = Ndays;

    ctrl->indewpt = 0;
    ctrl->invp = 0;
    if (have_dewpt) {
        if (have_dewpt == 1) {
            log_err("have_dewpt not yet implemented for tdew; however you can "
                    "supply observed vapor pressure and set have_dewpt to 2");
        }
        else if (have_dewpt == 2) {
            ctrl->invp = 1;
        }
    }
    if (have_shortwave) {
        ctrl->insw = 1;
    }
    else {
        ctrl->insw = 0;
    }
    ctrl->outhum = 1;           /* output vapor pressure */

    /* initialize the parameter structure.  Meteorological variables are only
       calculated for the mean grid cell elevation.  The temperatures are lapsed
       outside of the mtclim code.  Therefore p->base_elev and p->site_elev are
       set to the same value.  The same is true for p->base_isoh and
       p->site_isoh. */
    p->base_elev = elevation;
    p->base_isoh = annual_prcp / MM_PER_CM; /* MTCLIM prcp in cm */
    p->site_lat = lat;
    p->site_elev = elevation;
    p->site_slp = slope;
    p->site_asp = aspect;
    p->site_isoh = annual_prcp / MM_PER_CM; /* MTCLIM prcp in cm */
    p->site_ehoriz = ehoriz;
    p->site_whoriz = whoriz;
    p->tmax_lr = -1 * param.LAPSE_RATE * M_PER_KM;  /* not used since site_elev == base_elev */
    p->tmin_lr = -1 * param.LAPSE_RATE * M_PER_KM;  /* not used since site_elev == base_elev */

    /* allocate space in the data arrays for input and output data */
    if (data_alloc(ctrl, mtclim_data)) {
        log_err("Error in data_alloc()... exiting");
    }

    /* initialize the data arrays with the vic input data */
    for (i = 0; i < ctrl->ndays; i++) {
        mtclim_data->yday[i] = dmy[i * atmos_steps_per_day].day_in_year;
        mtclim_data->tmax[i] = tmax[i];
        mtclim_data->tmin[i] = tmin[i];
        if (ctrl->insw) {
            mtclim_data->s_srad[i] = 0;
            for (j = 0; j < atmos_steps_per_day; j++) {
                mtclim_data->s_srad[i] +=
                    subdailyrad[i * atmos_steps_per_day + j];
            }
            mtclim_data->s_srad[i] /= atmos_steps_per_day;
        }
        if (ctrl->invp) {
            mtclim_data->s_hum[i] = vp[i];
        }
        /* MTCLIM prcp in cm */
        mtclim_data->prcp[i] = prec[i] / MM_PER_CM;
        if (have_dewpt == 1) {
            log_err("have_dewpt not yet implemented ...");
        }
    }
    tinystepspday = (size_t)((CONST_CDAY) / param.MTCLIM_SRADDT);
    for (i = 0; i < DAYS_PER_LYEAR; i++) {
        for (j = 0; j < tinystepspday; j++) {
            tiny_radfract[i][j] = 0;
        }
    }
}

/******************************************************************************
 * @brief    Store MTCLIM variables in VIC arrays.
 *****************************************************************************/
void
mtclim_to_vic(double          sec_offset,
              dmy_struct     *dmy,
              double        **tiny_radfract,
              control_struct *ctrl,
              data_struct    *mtclim_data,
              double         *tskc,
              double         *vp,
              double         *subdailyrad,
              double         *fdir)
{
    extern parameters_struct   param;
    extern global_param_struct global_param;

    size_t                     i, j, k;
    size_t                     tinystepspstep;
    int                        tinystep;
    int                        tiny_offset;
    double                     tmp_rad;
    size_t                     atmos_steps_per_day;

    atmos_steps_per_day = global_param.model_steps_per_day;
    if (atmos_steps_per_day < (HOURS_PER_DAY)) {
        atmos_steps_per_day = (HOURS_PER_DAY);
    }

    tinystepspstep =
        (size_t)((CONST_CDAY) / param.MTCLIM_SRADDT) / atmos_steps_per_day;

    tiny_offset = (int)((double)tinystepspstep * sec_offset);
    for (i = 0; i < ctrl->ndays; i++) {
        if (ctrl->insw) {
            tmp_rad = mtclim_data->s_srad[i] * atmos_steps_per_day;
        }
        else {
            tmp_rad = mtclim_data->s_srad[i] * mtclim_data->s_dayl[i] /
                      (double) (SEC_PER_DAY);
        }
        for (j = 0; j < atmos_steps_per_day; j++) {
            subdailyrad[i * atmos_steps_per_day + j] = 0;
            for (k = 0; k < tinystepspstep; k++) {
                tinystep = j * tinystepspstep + k - tiny_offset;
                if (tinystep < 0) {
                    tinystep += atmos_steps_per_day * tinystepspstep;
                }
                if (tinystep >
                    (int)(atmos_steps_per_day * tinystepspstep - 1)) {
                    tinystep -= atmos_steps_per_day * tinystepspstep;
                }
                subdailyrad[i * atmos_steps_per_day + j] +=
                    tiny_radfract[dmy[i * atmos_steps_per_day +
                                      j].day_in_year - 1][tinystep];
            }
            subdailyrad[i * atmos_steps_per_day + j] *= tmp_rad;
        }
    }

    for (i = 0; i < ctrl->ndays; i++) {
        fdir[i] = mtclim_data->s_fdir[i];
        tskc[i] = mtclim_data->s_tskc[i];
        vp[i] = mtclim_data->s_hum[i];
    }
}
