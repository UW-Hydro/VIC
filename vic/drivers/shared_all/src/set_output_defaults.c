/******************************************************************************
 * @section DESCRIPTION
 *
 * Set the output_stream and out_data structures to default values. These can
 * be overridden by the user in the global control file.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Set the output_stream and out_data structures to default values.
             These can be overridden by the user in the global control file.
 *****************************************************************************/
void
get_default_nstreams_nvars(size_t *nstreams,
                           size_t  nvars[])
{
    extern option_struct options;

    size_t               streamnum;

    // Output files
    (*nstreams) = 2;
    if (options.FROZEN_SOIL) {
        (*nstreams)++;
    }
    if (options.SNOW_BAND) {
        (*nstreams)++;
    }
    if (options.LAKES) {
        (*nstreams)++;
    }

    streamnum = 0;
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        nvars[streamnum] = 26;
    }
    else {
        nvars[streamnum] = 20;
    }

    streamnum++;
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        nvars[streamnum] = 14;
    }
    else {
        nvars[streamnum] = 4;
    }
    if (options.BLOWING) {
        nvars += 3;
    }

    if (options.FROZEN_SOIL) {
        streamnum++;
        nvars[streamnum] = 4;
    }
    if (options.SNOW_BAND) {
        streamnum++;
        if (options.FULL_ENERGY) {
            nvars[streamnum] = 13;
        }
        else {
            nvars[streamnum] = 9;
        }
    }
    if (options.LAKES) {
        streamnum++;
        nvars[streamnum] = 8;
    }
}

/******************************************************************************
 * @brief    Set the output_stream and out_data structures to default values.
             These can be overridden by the user in the global control file.
 *****************************************************************************/
void
set_output_defaults(stream_struct **streams,
                    dmy_struct     *dmy_current,
                    unsigned short  default_file_format)
{
    extern option_struct options;

    size_t               streamnum;
    size_t               varnum;
    alarm_struct         default_alarm;
    int                  default_freq_n = 1;


    set_alarm(dmy_current, FREQ_NDAYS, &default_freq_n, &default_alarm);

    for (streamnum = 0; streamnum < options.Noutstreams; streamnum++) {
        (*streams)[streamnum].agg_alarm = default_alarm;
        (*streams)[streamnum].file_format = default_file_format;
    }

    // Variables in first file
    streamnum = 0;
    varnum = 0;
    strcpy((*streams)[streamnum].prefix, "fluxes");
    set_output_var(&((*streams)[streamnum]), "OUT_PREC", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_RUNOFF", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_BASEFLOW", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_WDEW", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SOIL_LIQ", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_RAD_TEMP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_SWNET", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_R_NET", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_LATENT", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP_CANOP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_TRANSP_VEG", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP_BARE", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SUB_CANOP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SUB_SNOW", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_SENSIBLE", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_GRND_FLUX", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_DELTAH", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_FUSION", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_AERO_RESIST", varnum++,
                   "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SURF_TEMP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_ALBEDO", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_REL_HUMID", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_IN_LONG", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_AIR_TEMP", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_WIND", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);

    // Variables in second file
    streamnum++;
    varnum = 0;
    strcpy((*streams)[streamnum].prefix, "snow");
    set_output_var(&((*streams)[streamnum]), "OUT_SWE", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_DEPTH", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_CANOPY", varnum++,
                   "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_COVER", varnum++, "%.4f",
                   OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_ADVECTION", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_DELTACC", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_FLUX", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_RFRZ_ENERGY", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_MELT_ENERGY", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_ADV_SENS", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LATENT_SUB", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_SURF_TEMP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_PACK_TEMP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_MELT", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    if (options.BLOWING) {
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_BLOWING", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_SURFACE", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_SNOW", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }

    // Variables in other files
    if (options.FROZEN_SOIL) {
        streamnum++;
        varnum = 0;
        strcpy((*streams)[streamnum].prefix, "fdepth");
        set_output_var(&((*streams)[streamnum]), "OUT_FDEPTH", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_TDEPTH", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SOIL_MOIST", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SURF_FROST_FRAC",
                       varnum++, "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    if (options.SNOW_BAND) {
        streamnum++;
        varnum = 0;
        strcpy((*streams)[streamnum].prefix, "snowband");
        set_output_var(&((*streams)[streamnum]), "OUT_SWE_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_DEPTH_BAND",
                       varnum++, "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_CANOPY_BAND",
                       varnum++, "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        if (options.FULL_ENERGY) {
            set_output_var(&((*streams)[streamnum]), "OUT_ADVECTION_BAND",
                           varnum++, "%.4f", OUT_TYPE_FLOAT, 1,
                           AGG_TYPE_DEFAULT);
            set_output_var(&((*streams)[streamnum]), "OUT_DELTACC_BAND",
                           varnum++, "%.4f", OUT_TYPE_FLOAT, 1,
                           AGG_TYPE_DEFAULT);
            set_output_var(&((*streams)[streamnum]), "OUT_SNOW_FLUX_BAND",
                           varnum++, "%.4f", OUT_TYPE_FLOAT, 1,
                           AGG_TYPE_DEFAULT);
            set_output_var(&((*streams)[streamnum]), "OUT_RFRZ_ENERGY_BAND",
                           varnum++, "%.4f", OUT_TYPE_FLOAT, 1,
                           AGG_TYPE_DEFAULT);
        }
        set_output_var(&((*streams)[streamnum]), "OUT_SWNET_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LWNET_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_ALBEDO_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LATENT_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SENSIBLE_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_GRND_FLUX_BAND", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
    if (options.LAKES) {
        streamnum++;
        varnum = 0;
        strcpy((*streams)[streamnum].prefix, "lake");
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_TEMP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_HEIGHT",
                       varnum++, "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_FRACT", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_DEPTH", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_SURF_AREA", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_VOLUME", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_SURF_TEMP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_EVAP", varnum++,
                       "%.4f", OUT_TYPE_FLOAT, 1, AGG_TYPE_DEFAULT);
    }
}
