/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
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
 * @brief    Print atmos data structure.
 *****************************************************************************/
void
print_atmos_data(force_data_struct *force,
                 size_t             nr)
{
    extern option_struct options;
    size_t               i;

    fprintf(LOG_DEST, "atmos_data  :\n");
    fprintf(LOG_DEST, "\tair_temp  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->air_temp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdensity   :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tlongwave  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->longwave[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_prec  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->out_prec);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_rain  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->out_rain);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tout_snow  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->out_snow);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tprec      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->prec[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tpressure  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->pressure[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tshortwave :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->shortwave[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsnowflag  :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%d\n", force->snowflag[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tvp        :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->vp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tvpd       :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->vpd[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\twind      :");
    for (i = 0; i <= nr; i++) {
        fprintf(LOG_DEST, "\t%.4f", force->wind[i]);
    }
    fprintf(LOG_DEST, "\n");
    if (options.LAKES) {
        fprintf(LOG_DEST, "\tchannel_in:");
        for (i = 0; i <= nr; i++) {
            fprintf(LOG_DEST, "\t%.4f", force->channel_in[i]);
        }
        fprintf(LOG_DEST, "\n");
    }
    if (options.CARBON) {
        fprintf(LOG_DEST, "\tCatm      :");
        for (i = 0; i <= nr; i++) {
            fprintf(LOG_DEST, "\t%.4f", force->Catm[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\tfdir      :");
        for (i = 0; i <= nr; i++) {
            fprintf(LOG_DEST, "\t%.4f", force->fdir[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\tpar       :");
        for (i = 0; i <= nr; i++) {
            fprintf(LOG_DEST, "\t%.4f", force->par[i]);
        }
        fprintf(LOG_DEST, "\n");
    }
}

/******************************************************************************
 * @brief    Print filenames structure.
 *****************************************************************************/
void
print_filenames(filenames_struct *fnames)
{
    fprintf(LOG_DEST, "filenames:\n");
    fprintf(LOG_DEST, "\tforcing[0]   : %s\n", fnames->forcing[0]);
    fprintf(LOG_DEST, "\tforcing[1]   : %s\n", fnames->forcing[1]);
    fprintf(LOG_DEST, "\tf_path_pfx[0]: %s\n", fnames->f_path_pfx[0]);
    fprintf(LOG_DEST, "\tf_path_pfx[1]: %s\n", fnames->f_path_pfx[1]);
    fprintf(LOG_DEST, "\tglobal       : %s\n", fnames->global);
    fprintf(LOG_DEST, "\tconstants    : %s\n", fnames->constants);
    fprintf(LOG_DEST, "\tinit_state   : %s\n", fnames->init_state);
    fprintf(LOG_DEST, "\tlakeparam    : %s\n", fnames->lakeparam);
    fprintf(LOG_DEST, "\tresult_dir   : %s\n", fnames->result_dir);
    fprintf(LOG_DEST, "\tsnowband     : %s\n", fnames->snowband);
    fprintf(LOG_DEST, "\tsoil         : %s\n", fnames->soil);
    fprintf(LOG_DEST, "\tstatefile    : %s\n", fnames->statefile);
    fprintf(LOG_DEST, "\tveg          : %s\n", fnames->veg);
    fprintf(LOG_DEST, "\tveglib       : %s\n", fnames->veglib);
    fprintf(LOG_DEST, "\tlog_path     : %s\n", fnames->log_path);
}

/******************************************************************************
 * @brief    Print file path structure.
 *****************************************************************************/
void
print_filep(filep_struct *fp)
{
    fprintf(LOG_DEST, "filep:\n");
    fprintf(LOG_DEST, "\tforcing[0] : %p\n", fp->forcing[0]);
    fprintf(LOG_DEST, "\tforcing[1] : %p\n", fp->forcing[1]);
    fprintf(LOG_DEST, "\tglobalparam: %p\n", fp->globalparam);
    fprintf(LOG_DEST, "\tconstants  : %p\n", fp->constants);
    fprintf(LOG_DEST, "\tinit_state : %p\n", fp->init_state);
    fprintf(LOG_DEST, "\tlakeparam  : %p\n", fp->lakeparam);
    fprintf(LOG_DEST, "\tsnowband   : %p\n", fp->snowband);
    fprintf(LOG_DEST, "\tsoilparam  : %p\n", fp->soilparam);
    fprintf(LOG_DEST, "\tstatefile  : %p\n", fp->statefile);
    fprintf(LOG_DEST, "\tveglib     : %p\n", fp->veglib);
    fprintf(LOG_DEST, "\tvegparam   : %p\n", fp->vegparam);
    fprintf(LOG_DEST, "\tlogfile    : %p\n", fp->logfile);
}
