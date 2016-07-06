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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(atmos_data_struct *atmos)
{
    extern option_struct options;

    atmos->air_temp = calloc(NR + 1, sizeof(*(atmos->air_temp)));
    check_alloc_status(atmos->air_temp, "Memory allocation error.");

    atmos->density = calloc(NR + 1, sizeof(*(atmos->density)));
    check_alloc_status(atmos->density, "Memory allocation error.");

    atmos->longwave = calloc(NR + 1, sizeof(*(atmos->longwave)));
    check_alloc_status(atmos->longwave, "Memory allocation error.");

    atmos->prec = calloc(NR + 1, sizeof(*(atmos->prec)));
    check_alloc_status(atmos->prec, "Memory allocation error.");

    atmos->pressure = calloc(NR + 1, sizeof(*(atmos->pressure)));
    check_alloc_status(atmos->pressure, "Memory allocation error.");

    atmos->shortwave = calloc(NR + 1, sizeof(*(atmos->shortwave)));
    check_alloc_status(atmos->shortwave, "Memory allocation error.");

    atmos->snowflag = calloc(NR + 1, sizeof(*(atmos->snowflag)));
    check_alloc_status(atmos->snowflag, "Memory allocation error.");

    atmos->vp = calloc(NR + 1, sizeof(*(atmos->vp)));
    check_alloc_status(atmos->vp, "Memory allocation error.");

    atmos->vpd = calloc(NR + 1, sizeof(*(atmos->vpd)));
    check_alloc_status(atmos->vpd, "Memory allocation error.");

    atmos->wind = calloc(NR + 1, sizeof(*(atmos->wind)));
    check_alloc_status(atmos->wind, "Memory allocation error.");

    if (options.LAKES) {
        atmos->channel_in = calloc(NR + 1, sizeof(*(atmos->channel_in)));
        check_alloc_status(atmos->channel_in, "Memory allocation error.");
    }
    if (options.CARBON) {
        atmos->Catm = calloc(NR + 1, sizeof(*(atmos->Catm)));
        check_alloc_status(atmos->Catm, "Memory allocation error.");

        atmos->coszen = calloc(NR + 1, sizeof(*(atmos->coszen)));
        check_alloc_status(atmos->coszen, "Memory allocation error.");

        atmos->fdir = calloc(NR + 1, sizeof(*(atmos->fdir)));
        check_alloc_status(atmos->fdir, "Memory allocation error.");

        atmos->par = calloc(NR + 1, sizeof(*(atmos->par)));
        check_alloc_status(atmos->par, "Memory allocation error.");
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
    free(atmos->density);
    free(atmos->longwave);
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
        free(atmos->coszen);
        free(atmos->fdir);
        free(atmos->par);
    }
}
