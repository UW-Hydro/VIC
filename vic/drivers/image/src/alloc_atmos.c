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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

/******************************************************************************
 * @brief    Allocate memory for the atmos data structure.
 *****************************************************************************/
void
alloc_atmos(atmos_data_struct *atmos)
{
    atmos->air_temp = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->air_temp == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->Catm = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->Catm == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->channel_in = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->channel_in == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->coszen = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->coszen == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->density = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->density == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->fdir = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->fdir == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->longwave = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->longwave == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->par = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->par == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->prec = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->prec == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->pressure = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->pressure == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->shortwave = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->shortwave == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->snowflag = (bool *) calloc(NR + 1, sizeof(bool));
    if (atmos->snowflag == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->tskc = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->tskc == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->vp = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->vp == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->vpd = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->vpd == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
    atmos->wind = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->wind == NULL) {
        log_err("Memory allocation error in alloc_atmos().");
    }
}

/******************************************************************************
 * @brief    Free memory for the atmos data structure.
 *****************************************************************************/
void
free_atmos(atmos_data_struct *atmos)
{
    if (atmos == NULL) {
        return;
    }

    free(atmos->air_temp);
    free(atmos->Catm);
    free(atmos->channel_in);
    free(atmos->coszen);
    free(atmos->density);
    free(atmos->fdir);
    free(atmos->longwave);
    free(atmos->par);
    free(atmos->prec);
    free(atmos->pressure);
    free(atmos->shortwave);
    free(atmos->snowflag);
    free(atmos->tskc);
    free(atmos->vp);
    free(atmos->vpd);
    free(atmos->wind);
}
