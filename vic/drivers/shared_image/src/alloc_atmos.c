/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the force data struct
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
 * @brief    Allocate memory for the force data structure.
 *****************************************************************************/
void
alloc_force(force_data_struct *force)
{
    extern option_struct options;

    force->air_temp = calloc(NR + 1, sizeof(*(force->air_temp)));
    check_alloc_status(force->air_temp, "Memory allocation error.");

    force->density = calloc(NR + 1, sizeof(*(force->density)));
    check_alloc_status(force->density, "Memory allocation error.");

    force->longwave = calloc(NR + 1, sizeof(*(force->longwave)));
    check_alloc_status(force->longwave, "Memory allocation error.");

    force->prec = calloc(NR + 1, sizeof(*(force->prec)));
    check_alloc_status(force->prec, "Memory allocation error.");

    force->pressure = calloc(NR + 1, sizeof(*(force->pressure)));
    check_alloc_status(force->pressure, "Memory allocation error.");

    force->shortwave = calloc(NR + 1, sizeof(*(force->shortwave)));
    check_alloc_status(force->shortwave, "Memory allocation error.");

    force->snowflag = calloc(NR + 1, sizeof(*(force->snowflag)));
    check_alloc_status(force->snowflag, "Memory allocation error.");

    force->vp = calloc(NR + 1, sizeof(*(force->vp)));
    check_alloc_status(force->vp, "Memory allocation error.");

    force->vpd = calloc(NR + 1, sizeof(*(force->vpd)));
    check_alloc_status(force->vpd, "Memory allocation error.");

    force->wind = calloc(NR + 1, sizeof(*(force->wind)));
    check_alloc_status(force->wind, "Memory allocation error.");

    if (options.LAKES) {
        force->channel_in = calloc(NR + 1, sizeof(*(force->channel_in)));
        check_alloc_status(force->channel_in, "Memory allocation error.");
    }
    if (options.CARBON) {
        force->Catm = calloc(NR + 1, sizeof(*(force->Catm)));
        check_alloc_status(force->Catm, "Memory allocation error.");

        force->coszen = calloc(NR + 1, sizeof(*(force->coszen)));
        check_alloc_status(force->coszen, "Memory allocation error.");

        force->fdir = calloc(NR + 1, sizeof(*(force->fdir)));
        check_alloc_status(force->fdir, "Memory allocation error.");

        force->par = calloc(NR + 1, sizeof(*(force->par)));
        check_alloc_status(force->par, "Memory allocation error.");
    }
}

/******************************************************************************
 * @brief    Free memory for the force data structure.
 *****************************************************************************/
void
free_force(force_data_struct *force)
{
    extern option_struct options;

    if (force == NULL) {
        return;
    }

    free(force->air_temp);
    free(force->density);
    free(force->longwave);
    free(force->prec);
    free(force->pressure);
    free(force->shortwave);
    free(force->snowflag);
    free(force->vp);
    free(force->vpd);
    free(force->wind);
    if (options.LAKES) {
        free(force->channel_in);
    }
    if (options.CARBON) {
        free(force->Catm);
        free(force->coszen);
        free(force->fdir);
        free(force->par);
    }
}
