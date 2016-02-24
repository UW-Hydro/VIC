/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the atmos data struct
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(atmos_data_struct *atmos)
{
    extern option_struct options;

    atmos->air_temp = calloc(NR + 1, sizeof(*(atmos->air_temp)));
    if (atmos->air_temp == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->density = calloc(NR + 1, sizeof(*(atmos->density)));
    if (atmos->density == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->longwave = calloc(NR + 1, sizeof(*(atmos->longwave)));
    if (atmos->longwave == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->prec = calloc(NR + 1, sizeof(*(atmos->prec)));
    if (atmos->prec == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->pressure = calloc(NR + 1, sizeof(*(atmos->pressure)));
    if (atmos->pressure == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->shortwave = calloc(NR + 1, sizeof(*(atmos->shortwave)));
    if (atmos->shortwave == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->snowflag = calloc(NR + 1, sizeof(*(atmos->snowflag)));
    if (atmos->snowflag == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->vp = calloc(NR + 1, sizeof(*(atmos->vp)));
    if (atmos->vp == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->vpd = calloc(NR + 1, sizeof(*(atmos->vpd)));
    if (atmos->vpd == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->wind = calloc(NR + 1, sizeof(*(atmos->wind)));
    if (atmos->wind == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    if (options.LAKES) {
        atmos->channel_in = calloc(NR + 1, sizeof(*(atmos->channel_in)));
        if (atmos->channel_in == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
    }
    if (options.CARBON) {
        atmos->Catm = calloc(NR + 1, sizeof(*(atmos->Catm)));
        if (atmos->Catm == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        atmos->fdir = calloc(NR + 1, sizeof(*(atmos->fdir)));
        if (atmos->fdir == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        atmos->par = calloc(NR + 1, sizeof(*(atmos->par)));
        if (atmos->par == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
    }
}

/******************************************************************************
 * @brief    Free memory for the atmos data structure.
 *****************************************************************************/
void
free_atmos(atmos_data_struct *atmos)
{
    extern option_struct options;

    if (atmos == NULL) {
        return;
    }

    free(atmos->air_temp);
    free(atmos->prec);
    free(atmos->pressure);
    free(atmos->shortwave);
    free(atmos->snowflag);
    free(atmos->vp);
    free(atmos->vpd);
    free(atmos->wind);
    if (options.LAKES) {
        free(atmos->channel_in);
    }
    if (options.CARBON) {
        free(atmos->Catm);
        free(atmos->fdir);
        free(atmos->par);
    }
}
