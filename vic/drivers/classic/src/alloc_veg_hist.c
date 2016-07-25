/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the veg_hist data struct
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
 * @brief    Allocate memory for veg his structure.
 *****************************************************************************/
void
alloc_veg_hist(int                nrecs,
               int                nveg,
               veg_hist_struct ***veg_hist)
{
    int i, j;

    *veg_hist = calloc(nrecs, sizeof(*(*veg_hist)));
    check_alloc_status((*veg_hist), "Memory allocation error.");

    for (i = 0; i < nrecs; i++) {
        (*veg_hist)[i] = calloc(nveg + 1, sizeof(*((*veg_hist)[i])));
        check_alloc_status((*veg_hist)[i], "Memory allocation error.");

        for (j = 0; j < nveg + 1; j++) {
            (*veg_hist)[i][j].albedo =
                calloc(NR + 1, sizeof(*((*veg_hist)[i][j].albedo)));
            check_alloc_status((*veg_hist)[i][j].albedo,
                               "Memory allocation error.");
            (*veg_hist)[i][j].displacement = calloc(NR + 1,
                                                    sizeof(*((*veg_hist)[i][j].
                                                             displacement)));
            check_alloc_status((*veg_hist)[i][j].displacement,
                               "Memory allocation error.");
            (*veg_hist)[i][j].fcanopy = calloc(NR + 1,
                                               sizeof(*((*veg_hist)[i][j].
                                                        fcanopy)));
            check_alloc_status((*veg_hist)[i][j].fcanopy,
                               "Memory allocation error.");
            (*veg_hist)[i][j].LAI =
                calloc(NR + 1, sizeof(*((*veg_hist)[i][j].LAI)));
            check_alloc_status((*veg_hist)[i][j].LAI,
                               "Memory allocation error.");
            (*veg_hist)[i][j].roughness = calloc(NR + 1,
                                                 sizeof(*((*veg_hist)[i][j].
                                                          roughness)));
            check_alloc_status((*veg_hist)[i][j].roughness,
                               "Memory allocation error.");
        }
    }
}

/******************************************************************************
 * @brief    Free veg hist structure.
 *****************************************************************************/
void
free_veg_hist(int                nrecs,
              int                nveg,
              veg_hist_struct ***veg_hist)
{
    int i, j;

    if (*veg_hist == NULL) {
        return;
    }

    for (i = 0; i < nrecs; i++) {
        for (j = 0; j < nveg + 1; j++) {
            free((*veg_hist)[i][j].albedo);
            free((*veg_hist)[i][j].displacement);
            free((*veg_hist)[i][j].fcanopy);
            free((*veg_hist)[i][j].LAI);
            free((*veg_hist)[i][j].roughness);
        }
        free((*veg_hist)[i]);
    }

    free(*veg_hist);
}
