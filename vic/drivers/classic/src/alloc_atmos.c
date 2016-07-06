/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the atmos data struct
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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
    check_alloc_status(*atmos, "Memory allocation error.");

    for (i = 0; i < nrecs; i++) {
        (*atmos)[i].air_temp = calloc(NR + 1, sizeof(*(*atmos)[i].air_temp));
        check_alloc_status((*atmos)[i].air_temp, "Memory allocation error.");
        (*atmos)[i].density = calloc(NR + 1, sizeof(*(*atmos)[i].density));
        check_alloc_status((*atmos)[i].density, "Memory allocation error.");
        (*atmos)[i].longwave = calloc(NR + 1, sizeof(*(*atmos)[i].longwave));
        check_alloc_status((*atmos)[i].longwave, "Memory allocation error.");
        (*atmos)[i].prec = calloc(NR + 1, sizeof(*(*atmos)[i].prec));
        check_alloc_status((*atmos)[i].prec, "Memory allocation error.");
        (*atmos)[i].pressure = calloc(NR + 1, sizeof(*(*atmos)[i].pressure));
        check_alloc_status((*atmos)[i].pressure, "Memory allocation error.");
        (*atmos)[i].shortwave = calloc(NR + 1, sizeof(*(*atmos)[i].shortwave));
        check_alloc_status((*atmos)[i].shortwave, "Memory allocation error.");
        (*atmos)[i].snowflag = calloc(NR + 1, sizeof(*(*atmos)[i].snowflag));
        check_alloc_status((*atmos)[i].snowflag, "Memory allocation error.");
        (*atmos)[i].vp = calloc(NR + 1, sizeof(*(*atmos)[i].vp));
        check_alloc_status((*atmos)[i].vp, "Memory allocation error.");
        (*atmos)[i].vpd = calloc(NR + 1, sizeof(*(*atmos)[i].vpd));
        check_alloc_status((*atmos)[i].vpd, "Memory allocation error.");
        (*atmos)[i].wind = calloc(NR + 1, sizeof(*(*atmos)[i].wind));
        check_alloc_status((*atmos)[i].wind, "Memory allocation error.");
        if (options.LAKES) {
            (*atmos)[i].channel_in =
                calloc(NR + 1, sizeof(*(*atmos)[i].channel_in));
            check_alloc_status((*atmos)[i].channel_in,
                               "Memory allocation error.");
        }
        if (options.CARBON) {
            (*atmos)[i].Catm = calloc(NR + 1, sizeof(*(*atmos)[i].Catm));
            check_alloc_status((*atmos)[i].Catm, "Memory allocation error.");
            (*atmos)[i].coszen = calloc(NR + 1, sizeof(*(*atmos)[i].coszen));
            check_alloc_status((*atmos)[i].coszen, "Memory allocation error.");
            (*atmos)[i].fdir = calloc(NR + 1, sizeof(*(*atmos)[i].fdir));
            check_alloc_status((*atmos)[i].fdir, "Memory allocation error.");
            (*atmos)[i].par = calloc(NR + 1, sizeof(*(*atmos)[i].par));
            check_alloc_status((*atmos)[i].par, "Memory allocation error.");
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
            free((*atmos)[i].coszen);
            free((*atmos)[i].fdir);
            free((*atmos)[i].par);
        }
    }

    free(*atmos);
}
