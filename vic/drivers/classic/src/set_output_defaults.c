/******************************************************************************
 * @section DESCRIPTION
 *
 * Set the out_data_files and out_data structures to default values. These can
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Set the out_data_files and out_data structures to default values.
             These can be overridden by the user in the global control file.
 *****************************************************************************/
out_data_file_struct *
set_output_defaults(out_data_struct *out_data)
{
    extern option_struct  options;
    out_data_file_struct *out_data_files;
    unsigned int          filenum;
    unsigned int          varnum;

    // Output files
    options.Noutfiles = 2;
    if (options.FROZEN_SOIL) {
        options.Noutfiles++;
    }
    if (options.PRT_SNOW_BAND) {
        options.Noutfiles++;
    }
    if (options.LAKES) {
        options.Noutfiles++;
    }
    out_data_files = calloc(options.Noutfiles, sizeof(*out_data_files));
    filenum = 0;
    strcpy(out_data_files[filenum].prefix, "fluxes");
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        out_data_files[filenum].nvars = 26;
    }
    else {
        out_data_files[filenum].nvars = 20;
    }
    filenum++;
    strcpy(out_data_files[filenum].prefix, "snow");
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        out_data_files[filenum].nvars = 14;
    }
    else {
        out_data_files[filenum].nvars = 4;
    }
    if (options.BLOWING) {
        out_data_files[filenum].nvars += 3;
    }
    if (options.FROZEN_SOIL) {
        filenum++;
        strcpy(out_data_files[filenum].prefix, "fdepth");
        out_data_files[filenum].nvars = 4;
    }
    if (options.PRT_SNOW_BAND) {
        filenum++;
        strcpy(out_data_files[filenum].prefix, "snowband");
        if (options.FULL_ENERGY) {
            out_data_files[filenum].nvars = 13;
        }
        else {
            out_data_files[filenum].nvars = 9;
        }
    }
    if (options.LAKES) {
        filenum++;
        strcpy(out_data_files[filenum].prefix, "lake");
        out_data_files[filenum].nvars = 8;
    }
    for (filenum = 0; filenum < options.Noutfiles; filenum++) {
        out_data_files[filenum].varid = calloc(
            out_data_files[filenum].nvars,
            sizeof(*(out_data_files[filenum].varid)));
    }

    // Variables in first file
    filenum = 0;
    varnum = 0;
    set_output_var(out_data_files, true, filenum, out_data, "OUT_PREC",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_EVAP",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_RUNOFF",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_BASEFLOW",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_WDEW",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SOIL_LIQ",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_RAD_TEMP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SWNET",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_R_NET",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LATENT", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_EVAP_CANOP", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_TRANSP_VEG", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_EVAP_BARE",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SUB_CANOP",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SUB_SNOW",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SENSIBLE", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_GRND_FLUX", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_DELTAH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_FUSION", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_AERO_RESIST", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SURF_TEMP",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_ALBEDO",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_REL_HUMID",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_IN_LONG",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_AIR_TEMP",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data, "OUT_WIND",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);

    // Variables in second file
    filenum++;
    varnum = 0;
    set_output_var(out_data_files, true, filenum, out_data, "OUT_SWE",
                   varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_SNOW_DEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_SNOW_CANOPY", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    set_output_var(out_data_files, true, filenum, out_data,
                   "OUT_SNOW_COVER", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    if (options.FULL_ENERGY || options.FROZEN_SOIL) {
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_ADVECTION", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_DELTACC", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_FLUX", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_RFRZ_ENERGY", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_MELT_ENERGY", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_ADV_SENS", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LATENT_SUB", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_SURF_TEMP", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_PACK_TEMP", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_MELT", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
    }
    if (options.BLOWING) {
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SUB_BLOWING", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SUB_SURFACE", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SUB_SNOW", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
    }

    // Variables in other files
    if (options.FROZEN_SOIL) {
        filenum++;
        varnum = 0;
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_FDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_TDEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SOIL_MOIST", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SURF_FROST_FRAC", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
    }
    if (options.PRT_SNOW_BAND) {
        filenum++;
        varnum = 0;
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SWE_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_DEPTH_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SNOW_CANOPY_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        if (options.FULL_ENERGY) {
            set_output_var(out_data_files, true, filenum, out_data,
                           "OUT_ADVECTION_BAND", varnum++, "%.4f",
                           OUT_TYPE_FLOAT, 1);
            set_output_var(out_data_files, true, filenum, out_data,
                           "OUT_DELTACC_BAND", varnum++, "%.4f",
                           OUT_TYPE_FLOAT, 1);
            set_output_var(out_data_files, true, filenum, out_data,
                           "OUT_SNOW_FLUX_BAND", varnum++, "%.4f",
                           OUT_TYPE_FLOAT, 1);
            set_output_var(out_data_files, true, filenum, out_data,
                           "OUT_RFRZ_ENERGY_BAND", varnum++, "%.4f",
                           OUT_TYPE_FLOAT, 1);
        }
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SWNET_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LWNET_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_ALBEDO_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LATENT_BAND", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_SENSIBLE_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_GRND_FLUX_BAND", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
    }
    if (options.LAKES) {
        filenum++;
        varnum = 0;
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_ICE_TEMP", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_ICE_HEIGHT", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_ICE_FRACT", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_DEPTH", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_SURF_AREA", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_VOLUME", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_SURF_TEMP", varnum++, "%.4f",
                       OUT_TYPE_FLOAT, 1);
        set_output_var(out_data_files, true, filenum, out_data,
                       "OUT_LAKE_EVAP", varnum++, "%.4f", OUT_TYPE_FLOAT,
                       1);
    }

    return out_data_files;
}
