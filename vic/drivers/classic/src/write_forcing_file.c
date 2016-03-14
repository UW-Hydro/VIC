/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine writes the complete forcing data files for use in future
 * simulations.
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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Write forcing data files for use in future simulations.
 *****************************************************************************/
void
write_forcing_file(atmos_data_struct    *atmos,
                   int                   nrecs,
                   out_data_file_struct *out_data_files,
                   out_data_struct      *out_data,
                   dmy_struct           *dmy)
{
    extern global_param_struct global_param;
    extern option_struct       options;
    extern parameters_struct   param;

    int                        rec, v;
    unsigned                   i;
    size_t                     j;
    double                     dt_sec;

    dt_sec = global_param.dt;

    for (rec = 0; rec < nrecs; rec++) {
        for (j = 0; j < NF; j++) {
            out_data[OUT_AIR_TEMP].data[0] = atmos[rec].air_temp[j];
            out_data[OUT_DENSITY].data[0] = atmos[rec].density[j];
            out_data[OUT_LWDOWN].data[0] = atmos[rec].longwave[j];
            out_data[OUT_PREC].data[0] = atmos[rec].prec[j];
            out_data[OUT_PRESSURE].data[0] = atmos[rec].pressure[j] /
                                             PA_PER_KPA;
            out_data[OUT_QAIR].data[0] = CONST_EPS * atmos[rec].vp[j] /
                                         atmos[rec].pressure[j];
            out_data[OUT_REL_HUMID].data[0] = FRACT_TO_PERCENT *
                                              atmos[rec].vp[j] /
                                              (atmos[rec].vp[j] +
                                               atmos[rec].vpd[j]);
            out_data[OUT_SWDOWN].data[0] = atmos[rec].shortwave[j];
            out_data[OUT_VP].data[0] = atmos[rec].vp[j] / PA_PER_KPA;
            out_data[OUT_VPD].data[0] = atmos[rec].vpd[j] / PA_PER_KPA;
            out_data[OUT_WIND].data[0] = atmos[rec].wind[j];
            if (out_data[OUT_AIR_TEMP].data[0] >= param.SNOW_MAX_SNOW_TEMP) {
                out_data[OUT_RAINF].data[0] = out_data[OUT_PREC].data[0];
                out_data[OUT_SNOWF].data[0] = 0;
            }
            else if (out_data[OUT_AIR_TEMP].data[0] <=
                     param.SNOW_MIN_RAIN_TEMP) {
                out_data[OUT_RAINF].data[0] = 0;
                out_data[OUT_SNOWF].data[0] = out_data[OUT_PREC].data[0];
            }
            else {
                out_data[OUT_RAINF].data[0] =
                    ((out_data[OUT_AIR_TEMP].data[0] -
                      param.SNOW_MIN_RAIN_TEMP) /
                     (param.SNOW_MAX_SNOW_TEMP -
                      param.SNOW_MIN_RAIN_TEMP)) * out_data[OUT_PREC].data[0];
                out_data[OUT_SNOWF].data[0] = out_data[OUT_PREC].data[0] -
                                              out_data[OUT_RAINF].data[0];
            }
            if (options.CARBON) {
                out_data[OUT_CATM].data[0] = atmos[rec].Catm[j] /
                                             PPM_to_MIXRATIO;
                out_data[OUT_FDIR].data[0] = atmos[rec].fdir[j];
                out_data[OUT_PAR].data[0] = atmos[rec].par[j];
            }
            else {
                out_data[OUT_CATM].data[0] = MISSING;
                out_data[OUT_FDIR].data[0] = MISSING;
                out_data[OUT_PAR].data[0] = MISSING;
            }

            for (v = 0; v < N_OUTVAR_TYPES; v++) {
                for (i = 0; i < out_data[v].nelem; i++) {
                    out_data[v].aggdata[i] = out_data[v].data[i];
                }
            }

            if (options.ALMA_OUTPUT) {
                out_data[OUT_PREC].aggdata[0] /= dt_sec;
                out_data[OUT_RAINF].aggdata[0] /= dt_sec;
                out_data[OUT_SNOWF].aggdata[0] /= dt_sec;
                out_data[OUT_AIR_TEMP].aggdata[0] += CONST_TKFRZ;
                out_data[OUT_PRESSURE].aggdata[0] *= PA_PER_KPA;
                out_data[OUT_VP].aggdata[0] *= PA_PER_KPA;
                out_data[OUT_VPD].aggdata[0] *= PA_PER_KPA;
            }

            if (options.BINARY_OUTPUT) {
                for (v = 0; v < N_OUTVAR_TYPES; v++) {
                    for (i = 0; i < out_data[v].nelem; i++) {
                        out_data[v].aggdata[i] *= out_data[v].mult;
                    }
                }
            }
            write_data(out_data_files, out_data, &dmy[rec], global_param.dt);
        }
    }
}
