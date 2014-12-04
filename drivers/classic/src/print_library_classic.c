/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
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
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Print atmos data structure.
 *****************************************************************************/
void
print_atmos_data(atmos_data_struct *atmos,
                 size_t             nr)
{
    size_t i;

    printf("atmos_data  :\n");
    printf("\tair_temp  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->air_temp[i]);
    }
    printf("\n");
    printf("\tCatm      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->Catm[i]);
    }
    printf("\n");
    printf("\tchannel_in:");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->channel_in[i]);
    }
    printf("\n");
    printf("\tcoszen    :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->coszen[i]);
    }
    printf("\n");
    printf("\tdensity   :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->density[i]);
    }
    printf("\n");
    printf("\tfdir      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->fdir[i]);
    }
    printf("\n");
    printf("\tlongwave  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->longwave[i]);
    }
    printf("\n");
    printf("\tout_prec  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_prec);
    }
    printf("\n");
    printf("\tout_rain  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_rain);
    }
    printf("\n");
    printf("\tout_snow  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_snow);
    }
    printf("\n");
    printf("\tpar       :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->par[i]);
    }
    printf("\n");
    printf("\tprec      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->prec[i]);
    }
    printf("\n");
    printf("\tpressure  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->pressure[i]);
    }
    printf("\n");
    printf("\tshortwave :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->shortwave[i]);
    }
    printf("\n");
    printf("\tsnowflag  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%d\n", atmos->snowflag[i]);
    }
    printf("\n");
    printf("\ttskc      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->tskc[i]);
    }
    printf("\n");
    printf("\tvp        :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->vp[i]);
    }
    printf("\n");
    printf("\tvpd       :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->vpd[i]);
    }
    printf("\n");
    printf("\twind      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->wind[i]);
    }
    printf("\n");
}
