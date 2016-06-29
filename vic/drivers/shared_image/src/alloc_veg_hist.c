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

 #include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief
 *****************************************************************************/
void
alloc_veg_hist(veg_hist_struct *veg_hist)
{
    veg_hist->albedo = calloc(NR + 1, sizeof(*(veg_hist->albedo)));
    check_alloc_status(veg_hist->albedo, "Memory allocation error.");

    veg_hist->displacement = calloc(NR + 1, sizeof(*(veg_hist->displacement)));
    check_alloc_status(veg_hist->displacement, "Memory allocation error.");

    veg_hist->fcanopy = calloc(NR + 1, sizeof(*(veg_hist->fcanopy)));
    check_alloc_status(veg_hist->fcanopy, "Memory allocation error.");

    veg_hist->LAI = calloc(NR + 1, sizeof(*(veg_hist->LAI)));
    check_alloc_status(veg_hist->LAI, "Memory allocation error.");

    veg_hist->roughness = calloc(NR + 1, sizeof(*(veg_hist->roughness)));
    check_alloc_status(veg_hist->roughness, "Memory allocation error.");
}

/******************************************************************************
 * @brief    Free veg hist structure.
 *****************************************************************************/
void
free_veg_hist(veg_hist_struct *veg_hist)
{
    if (veg_hist == NULL) {
        return;
    }

    free(veg_hist->albedo);
    free(veg_hist->displacement);
    free(veg_hist->fcanopy);
    free(veg_hist->LAI);
    free(veg_hist->roughness);
}
