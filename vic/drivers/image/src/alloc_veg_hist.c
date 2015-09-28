/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the veg_hist data struct
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

/******************************************************************************
 * @brief
 *****************************************************************************/
void
alloc_veg_hist(veg_hist_struct *veg_hist)
{
    veg_hist->albedo = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->albedo == NULL) {
        log_err("Memory allocation error in alloc_veg_hist().");
    }
    veg_hist->LAI = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->LAI == NULL) {
        log_err("Memory allocation error in alloc_veg_hist().");
    }
    veg_hist->vegcover = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->vegcover == NULL) {
        log_err("Memory allocation error in alloc_veg_hist().");
    }
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
    free(veg_hist->LAI);
    free(veg_hist->vegcover);
}
