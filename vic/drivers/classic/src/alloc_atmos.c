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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(int                 nrecs,
            atmos_data_struct **atmos)
{
    extern option_struct options;
    int                  i;

    *atmos = calloc(nrecs, sizeof(atmos_data_struct));
    if (*atmos == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }

    for (i = 0; i < nrecs; i++) {
        (*atmos)[i].air_temp = calloc(NR + 1, sizeof(*(*atmos)[i].air_temp));
        if ((*atmos)[i].air_temp == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].density = calloc(NR + 1, sizeof(*(*atmos)[i].density));
        if ((*atmos)[i].density == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].longwave = calloc(NR + 1, sizeof(*(*atmos)[i].longwave));
        if ((*atmos)[i].longwave == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].prec = calloc(NR + 1, sizeof(*(*atmos)[i].prec));
        if ((*atmos)[i].prec == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].pressure = calloc(NR + 1, sizeof(*(*atmos)[i].pressure));
        if ((*atmos)[i].pressure == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].shortwave = calloc(NR + 1, sizeof(*(*atmos)[i].shortwave));
        if ((*atmos)[i].shortwave == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].snowflag = calloc(NR + 1, sizeof(*(*atmos)[i].snowflag));
        if ((*atmos)[i].snowflag == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].vp = calloc(NR + 1, sizeof(*(*atmos)[i].vp));
        if ((*atmos)[i].vp == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].vpd = calloc(NR + 1, sizeof(*(*atmos)[i].vpd));
        if ((*atmos)[i].vpd == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].wind = calloc(NR + 1, sizeof(*(*atmos)[i].wind));
        if ((*atmos)[i].wind == NULL) {
            log_err("Memory allocation error in alloc_atmos().");
        }
        if (options.LAKES) {
            (*atmos)[i].channel_in =
                calloc(NR + 1, sizeof(*(*atmos)[i].channel_in));
            if ((*atmos)[i].channel_in == NULL) {
                log_err("Memory allocation error in alloc_atmos().");
            }
        }
        if (options.CARBON) {
            (*atmos)[i].Catm = calloc(NR + 1, sizeof(*(*atmos)[i].Catm));
            if ((*atmos)[i].Catm == NULL) {
                log_err("Memory allocation error in alloc_atmos().");
            }
            (*atmos)[i].fdir = calloc(NR + 1, sizeof(*(*atmos)[i].fdir));
            if ((*atmos)[i].fdir == NULL) {
                log_err("Memory allocation error in alloc_atmos().");
            }
            (*atmos)[i].par = calloc(NR + 1, sizeof(*(*atmos)[i].par));
            if ((*atmos)[i].par == NULL) {
                log_err("Memory allocation error in alloc_atmos().");
            }
        }
    }
}

/******************************************************************************
 * @brief    Free memory for the atmos data structure.
 *****************************************************************************/
void
free_atmos(int                 nrecs,
           atmos_data_struct **atmos)
{
    extern option_struct options;
    int                  i;

    if (*atmos == NULL) {
        return;
    }

    for (i = 0; i < nrecs; i++) {
        free((*atmos)[i].air_temp);
        free((*atmos)[i].density);
        free((*atmos)[i].longwave);
        free((*atmos)[i].prec);
        free((*atmos)[i].pressure);
        free((*atmos)[i].shortwave);
        free((*atmos)[i].snowflag);
        free((*atmos)[i].vp);
        free((*atmos)[i].vpd);
        free((*atmos)[i].wind);
        if (options.LAKES) {
            free((*atmos)[i].channel_in);
        }
        if (options.CARBON) {
            free((*atmos)[i].Catm);
            free((*atmos)[i].fdir);
            free((*atmos)[i].par);
        }
    }

    free(*atmos);
}
