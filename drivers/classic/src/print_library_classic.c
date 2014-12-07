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

    fprintf(LOG_DEST, "atmos_data  :\n");
    fprintf(LOG_DEST, "\tair_temp  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->air_temp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCatm      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->Catm[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tchannel_in:");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->channel_in[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tcoszen    :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->coszen[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdensity   :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfdir      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->fdir[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tlongwave  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->longwave[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_prec  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->out_prec);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_rain  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->out_rain);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_snow  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->out_snow);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tpar       :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->par[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tprec      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->prec[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tpressure  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->pressure[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tshortwave :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->shortwave[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsnowflag  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%d\n", atmos->snowflag[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttskc      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->tskc[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tvp        :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->vp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tvpd       :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->vpd[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\twind      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4lf", atmos->wind[i]);
    }
    fprintf(LOG_DEST, "\n");
}
