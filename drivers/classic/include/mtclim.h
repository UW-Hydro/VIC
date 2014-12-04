/******************************************************************************
 * @section DESCRIPTION
 *
 * MTCLIM Constants header file.
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

#ifndef MTCLIM_CONSTANTS_VIC_H
#define MTCLIM_CONSTANTS_VIC_H

/******************************************************************************
 * @brief    MTCLIM control structure
 *****************************************************************************/
typedef struct
{
    int ndays;           /**< number of days of data in input file */
    int insw;            /**< input shortwave radiation flag (0=NO, 1=YES) */
    int indewpt;         /**< input dewpoint temperature flag (0=NO, 1=YES) */
    int invp;            /**< input vapor pressure flag (0=NO, 1=YES) */
    int outhum;          /**< output humidity flag            (0=VPD, 1=VP) */
    int inyear;          /**< input year flag                 (0=NO, 1=YES) */
} control_struct;

/******************************************************************************
 * @brief    MTCLIM parameter structure
 *****************************************************************************/
typedef struct
{
    double base_elev;    /**< base elevation, meters */
    double base_isoh;    /**< base annual precip isohyet, cm */
    double site_lat;     /**< site latitude, dec. degrees (- for south) */
    double site_elev;    /**< site elevation, meters */
    double site_slp;     /**< site slope, degrees */
    double site_asp;     /**< site aspect, degrees */
    double site_isoh;    /**< site annual precip isohyet, cm */
    double site_ehoriz;  /**< site east horizon, degrees */
    double site_whoriz;  /**< site west horizon, degrees */
    double tmax_lr;      /**< maximum temperature lapse rate, deg C/1000m */
    double tmin_lr;      /**< minimum temperature lapse rate, deg C/1000m */
} parameter_struct;

/******************************************************************************
 * @brief    MTCLIM data structure
 *****************************************************************************/
typedef struct
{
    int *year;           /**< array of year values */
    int *yday;           /**< array of yearday values */
    double *tmax;        /**< array of base maximum temperature values */
    double *tmin;        /**< array of base minimum temperature values */
    double *prcp;        /**< array of base daily precipitation values */
    double *tdew;        /**< array of base dewpoint temperature values */
    double *s_tmax;      /**< array of site tmax values */
    double *s_tmin;      /**< array of site tmin values */
    double *s_tday;      /**< array of site daylight temperature values */
    double *s_prcp;      /**< array of site prcp values */
    double *s_hum;       /**< array of site humidity values (VPD or VP, Pa) */
    double *s_srad;      /**< array of site shortwave radiation values */
    double *s_dayl;      /**< array of site daylength values */
    double *s_swe;       /**< array of site snowpack values */
    /* start vic_change */
    double *s_fdir;            /**< array of site values of direct fraction of shortwave radiation */
    double *s_tskc;            /**< array of cloudiness values */
    double *s_ppratio;   /**< array of pet/prcp ratio values */
    double *s_ttmax;     /**< array of clear sky transmittance values */
    double *s_tfmax;     /**< array of cloud transmittance factor values */
    /* end vic_change */
} data_struct;

int calc_tair(const control_struct *ctrl, const parameter_struct *p,
              data_struct *data);
int calc_prcp(const control_struct *ctrl, const parameter_struct *p,
              data_struct *data);
/* start vic_change */
int calc_srad_humidity_iterative(const control_struct *ctrl,
                                 const parameter_struct *p, data_struct *data,
                                 double **tiny_radfract);
int snowpack(const control_struct *ctrl, data_struct *data);
void compute_srad_humidity_onetime(int ndays, const control_struct *ctrl,
                                   data_struct *data, double *tdew, double *pva,
                                   double *ttmax0, double *flat_potrad,
                                   double *slope_potrad, double sky_prop,
                                   double *daylength, double *pet,
                                   double *parray, double pa, double *dtr);
/* end vic_change */
int data_alloc(const control_struct *ctrl, data_struct *data);
int data_free(const control_struct *ctrl, data_struct *data);
double calc_pet(double rad, double ta, double pa, double dayl);
double atm_pres(double elev);
int pulled_boxcar(double *input, double *output, int n, int w, int w_flag);

#endif
