/******************************************************************************
 * @section DESCRIPTION
 *
 * Set the output_stream and out_data structures to default values. These can
 * be overridden by the user in the global control file.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Set the output_stream and out_data structures to default values.
             These can be overridden by the user in the global control file.
 *****************************************************************************/
void
set_output_defaults(stream_struct **streams)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    size_t                     streamnum;
    size_t                     varnum;
    size_t                     nvars;
    unsigned int               nextagg;

    nextagg = global_param.model_steps_per_day;

    // Output files
    options.Noutstreams = 2;
    if (options.FROZEN_SOIL) {
        options.Noutstreams++;
    }
    if (options.LAKES) {
        options.Noutstreams++;
    }

    *streams = calloc(options.Noutstreams, sizeof(*(*streams)));
    if (*streams == NULL) {
        log_err("Memory allocation error in set_output_defaults().");
    }

    streamnum = 0;
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        nvars = 26;
    }
    else {
        nvars = 20;
    }
    setup_stream(&((*streams)[streamnum]), nvars);
    strcpy((*streams)[streamnum].prefix, "fluxes");
    (*streams)[streamnum].file_format = ASCII;

    streamnum++;
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        nvars = 14;
    }
    else {
        nvars = 4;
    }
    if (options.BLOWING) {
        nvars += 3;
    }
    setup_stream(&((*streams)[streamnum]), nvars);
    strcpy((*streams)[streamnum].prefix, "snow");
    (*streams)[streamnum].file_format = ASCII;

    if (options.FROZEN_SOIL) {
        streamnum++;
        nvars = 4;

        setup_stream(&((*streams)[streamnum]), nvars);
        strcpy((*streams)[streamnum].prefix, "fdepth");
        (*streams)[streamnum].file_format = ASCII;
    }
    if (options.LAKES) {
        streamnum++;
        nvars = 8;
        setup_stream(&((*streams)[streamnum]), nvars);
        strcpy((*streams)[streamnum].prefix, "lake");
        (*streams)[streamnum].file_format = ASCII;
    }

    // Variables in first file
    streamnum = 0;
    varnum = 0;
    set_output_var(&((*streams)[streamnum]), "OUT_PREC", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_RUNOFF", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_BASEFLOW", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_WDEW", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SOIL_LIQ", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_RAD_TEMP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_SWNET", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_R_NET", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_LATENT", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP_CANOP", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_TRANSP_VEG", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_EVAP_BARE", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SUB_CANOP", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SUB_SNOW", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_SENSIBLE", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_GRND_FLUX", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_DELTAH", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_FUSION", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }
    set_output_var(&((*streams)[streamnum]), "OUT_AERO_RESIST", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SURF_TEMP", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_ALBEDO", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_REL_HUMID", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_IN_LONG", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_AIR_TEMP", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_WIND", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);

    // Variables in second file
    streamnum++;
    varnum = 0;
    set_output_var(&((*streams)[streamnum]), "OUT_SWE", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_DEPTH", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_CANOPY", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    set_output_var(&((*streams)[streamnum]), "OUT_SNOW_COVER", varnum++,
                   OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE, OUT_MULT_DEFAULT,
                   AGG_TYPE_DEFAULT);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(&((*streams)[streamnum]), "OUT_ADVECTION", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_DELTACC", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_FLUX", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_RFRZ_ENERGY", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_MELT_ENERGY", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_ADV_SENS", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LATENT_SUB", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_SURF_TEMP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_PACK_TEMP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SNOW_MELT", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }
    if (options.BLOWING) {
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_BLOWING", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_SURFACE", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SUB_SNOW", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }

    // Variables in other files
    if (options.FROZEN_SOIL) {
        streamnum++;
        varnum = 0;
        set_output_var(&((*streams)[streamnum]), "OUT_FDEPTH", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_TDEPTH", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SOIL_MOIST", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_SURF_FROST_FRAC",
                       varnum++, OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }

    if (options.LAKES) {
        streamnum++;
        varnum = 0;
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_TEMP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_HEIGHT",
                       varnum++, OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_ICE_FRACT", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_DEPTH", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_SURF_AREA", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_VOLUME", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_SURF_TEMP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
        set_output_var(&((*streams)[streamnum]), "OUT_LAKE_EVAP", varnum++,
                       OUT_ASCII_FORMAT_DEFAULT, OUT_TYPE_DOUBLE,
                       OUT_MULT_DEFAULT, AGG_TYPE_DEFAULT);
    }
}
