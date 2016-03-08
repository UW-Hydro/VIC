/******************************************************************************
 * @section DESCRIPTION
 *
 * Setup netCDF output variables.
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
 * @brief    Setup netCDF output variables.
 *****************************************************************************/
void
vic_nc_info(nc_file_struct   *nc_hist_file,
            out_data_struct **out_data,
            nc_var_struct    *nc_vars)
{
    extern option_struct options;

    size_t               i;
    size_t               j;

    // default is a variable of type NC_DOUBLE with only a single field per
    // timestep
    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        strncpy(nc_vars[i].nc_var_name,
                out_data[0][i].varname, MAXSTRING);
        nc_vars[i].nc_write = out_data[0][i].write;
        nc_vars[i].nc_aggtype = out_data[0][i].aggtype;
        nc_vars[i].nc_type = NC_DOUBLE;
        for (j = 0; j < MAXDIMS; j++) {
            nc_vars[i].nc_dimids[j] = -1;
            nc_vars[i].nc_counts[j] = -1;
        }
        nc_vars[i].nc_dims = 3;
        nc_vars[i].nc_dimids[0] = nc_hist_file->time_dimid;
        nc_vars[i].nc_dimids[1] = nc_hist_file->nj_dimid;
        nc_vars[i].nc_counts[1] = nc_hist_file->nj_size;
        nc_vars[i].nc_dimids[2] = nc_hist_file->ni_dimid;
        nc_vars[i].nc_counts[2] = nc_hist_file->ni_size;
    }

    // Set the non-default ones

    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        switch (i) {
        case OUT_FDEPTH:
        case OUT_TDEPTH:
            nc_vars[i].nc_dims = 4;
            nc_vars[i].nc_dimids[0] = nc_hist_file->time_dimid;
            nc_vars[i].nc_dimids[1] = nc_hist_file->front_dimid;
            nc_vars[i].nc_counts[1] = nc_hist_file->front_size;
            nc_vars[i].nc_dimids[2] = nc_hist_file->nj_dimid;
            nc_vars[i].nc_counts[2] = nc_hist_file->nj_size;
            nc_vars[i].nc_dimids[3] = nc_hist_file->ni_dimid;
            nc_vars[i].nc_counts[3] = nc_hist_file->ni_size;
            break;
        case OUT_SMLIQFRAC:
        case OUT_SMFROZFRAC:
        case OUT_SOIL_ICE:
        case OUT_SOIL_LIQ:
        case OUT_SOIL_MOIST:
        case OUT_SOIL_TEMP:
            nc_vars[i].nc_dims = 4;
            nc_vars[i].nc_dimids[0] = nc_hist_file->time_dimid;
            nc_vars[i].nc_dimids[1] = nc_hist_file->layer_dimid;
            nc_vars[i].nc_counts[1] = nc_hist_file->layer_size;
            nc_vars[i].nc_dimids[2] = nc_hist_file->nj_dimid;
            nc_vars[i].nc_counts[2] = nc_hist_file->nj_size;
            nc_vars[i].nc_dimids[3] = nc_hist_file->ni_dimid;
            nc_vars[i].nc_counts[3] = nc_hist_file->ni_size;
            break;
        case OUT_SOIL_TNODE:
        case OUT_SOIL_TNODE_WL:
        case OUT_SOILT_FBFLAG:
            nc_vars[i].nc_dims = 4;
            nc_vars[i].nc_dimids[0] = nc_hist_file->time_dimid;
            nc_vars[i].nc_dimids[1] = nc_hist_file->node_dimid;
            nc_vars[i].nc_counts[1] = nc_hist_file->node_size;
            nc_vars[i].nc_dimids[2] = nc_hist_file->nj_dimid;
            nc_vars[i].nc_counts[2] = nc_hist_file->nj_size;
            nc_vars[i].nc_dimids[3] = nc_hist_file->ni_dimid;
            nc_vars[i].nc_counts[3] = nc_hist_file->ni_size;
            break;
        case OUT_ADV_SENS_BAND:
        case OUT_ADVECTION_BAND:
        case OUT_ALBEDO_BAND:
        case OUT_DELTACC_BAND:
        case OUT_GRND_FLUX_BAND:
        case OUT_IN_LONG_BAND:
        case OUT_LATENT_BAND:
        case OUT_LATENT_SUB_BAND:
        case OUT_MELT_ENERGY_BAND:
        case OUT_LWNET_BAND:
        case OUT_RFRZ_ENERGY_BAND:
        case OUT_SENSIBLE_BAND:
        case OUT_SNOW_CANOPY_BAND:
        case OUT_SNOW_COVER_BAND:
        case OUT_SNOW_DEPTH_BAND:
        case OUT_SNOW_FLUX_BAND:
        case OUT_SNOW_MELT_BAND:
        case OUT_SNOW_PACKT_BAND:
        case OUT_SNOW_SURFT_BAND:
        case OUT_SWE_BAND:
            nc_vars[i].nc_dims = 4;
            nc_vars[i].nc_dimids[0] = nc_hist_file->time_dimid;
            nc_vars[i].nc_dimids[1] = nc_hist_file->band_dimid;
            nc_vars[i].nc_counts[1] = nc_hist_file->band_size;
            nc_vars[i].nc_dimids[2] = nc_hist_file->nj_dimid;
            nc_vars[i].nc_counts[2] = nc_hist_file->nj_size;
            nc_vars[i].nc_dimids[3] = nc_hist_file->ni_dimid;
            nc_vars[i].nc_counts[3] = nc_hist_file->ni_size;
            break;
        }
    }


    // set the units
    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        switch (i) {
        // fraction or unitless
        case OUT_ASAT:
        case OUT_LAKE_AREA_FRAC:
        case OUT_LAKE_ICE_FRACT:
        case OUT_SMFROZFRAC:
        case OUT_SMLIQFRAC:
        case OUT_SNOW_COVER:
        case OUT_SOIL_WET:
        case OUT_SURF_FROST_FRAC:
        case OUT_ALBEDO:
        case OUT_SALBEDO:
        case OUT_SNOWT_FBFLAG:
        case OUT_SOILT_FBFLAG:
        case OUT_SURFT_FBFLAG:
        case OUT_TCAN_FBFLAG:
        case OUT_TFOL_FBFLAG:
        case OUT_FDIR:
        case OUT_REL_HUMID:
        case OUT_ALBEDO_BAND:
        case OUT_SNOW_COVER_BAND:
            strncpy(nc_vars[i].nc_units, "-", MAXSTRING);
            break;
        // depth in meters
        case OUT_LAKE_DEPTH:
        case OUT_LAKE_SWE:
            strncpy(nc_vars[i].nc_units, "m", MAXSTRING);
            break;
        // depth in cm
        case OUT_LAKE_ICE_HEIGHT:
        case OUT_SNOW_DEPTH:
        case OUT_ZWT:
        case OUT_ZWT_LUMPED:
        case OUT_SNOW_DEPTH_BAND:
            strncpy(nc_vars[i].nc_units, "cm", MAXSTRING);
            break;
        // depth in cm, changed to m for ALMA_OUTPUT
        case OUT_FDEPTH:
        case OUT_TDEPTH:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "m", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "cm", MAXSTRING);
            }
            break;
        // depth in mm
        case OUT_LAKE_ICE:
        case OUT_LAKE_MOIST:
        case OUT_ROOTMOIST:
        case OUT_SNOW_CANOPY:
        case OUT_SOIL_ICE:
        case OUT_SOIL_LIQ:
        case OUT_SOIL_MOIST:
        case OUT_SURFSTOR:
        case OUT_SWE:
        case OUT_WDEW:
        case OUT_DELINTERCEPT:
        case OUT_DELSOILMOIST:
        case OUT_DELSURFSTOR:
        case OUT_DELSWE:
        case OUT_WATER_ERROR:
        case OUT_SNOW_CANOPY_BAND:
        case OUT_SWE_BAND:
            strncpy(nc_vars[i].nc_units, "mm", MAXSTRING);
            break;
        // rate in mm per timestep, changed to mm/s for ALMA_OUTPUT
        case OUT_BASEFLOW:
        case OUT_EVAP:
        case OUT_EVAP_BARE:
        case OUT_EVAP_CANOP:
        case OUT_INFLOW:
        case OUT_LAKE_BF_IN:
        case OUT_LAKE_BF_OUT:
        case OUT_LAKE_CHAN_IN:
        case OUT_LAKE_CHAN_OUT:
        case OUT_LAKE_DSTOR:
        case OUT_LAKE_DSWE:
        case OUT_LAKE_EVAP:
        case OUT_LAKE_RCHRG:
        case OUT_LAKE_RO_IN:
        case OUT_LAKE_VAPFLX:
        case OUT_PET:
        case OUT_PREC:
        case OUT_RAINF:
        case OUT_REFREEZE:
        case OUT_RUNOFF:
        case OUT_SNOW_MELT:
        case OUT_SNOWF:
        case OUT_SUB_BLOWING:
        case OUT_SUB_CANOP:
        case OUT_SUB_SNOW:
        case OUT_SUB_SURFACE:
        case OUT_TRANSP_VEG:
        case OUT_SNOW_MELT_BAND:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "mm s-1", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "mm", MAXSTRING);
            }
            break;
        // area in m2
        case OUT_LAKE_SURF_AREA:
            strncpy(nc_vars[i].nc_units, "m2", MAXSTRING);
            break;
        // volume in m3
        case OUT_LAKE_SWE_V:
        case OUT_LAKE_VOLUME:
            strncpy(nc_vars[i].nc_units, "m3", MAXSTRING);
            break;
        // rate in m3 per timestep, changed to m3/s for ALMA_OUTPUT
        case OUT_LAKE_BF_IN_V:
        case OUT_LAKE_BF_OUT_V:
        case OUT_LAKE_CHAN_IN_V:
        case OUT_LAKE_CHAN_OUT_V:
        case OUT_LAKE_DSTOR_V:
        case OUT_LAKE_DSWE_V:
        case OUT_LAKE_EVAP_V:
        case OUT_LAKE_PREC_V:
        case OUT_LAKE_RCHRG_V:
        case OUT_LAKE_RO_IN_V:
        case OUT_LAKE_VAPFLX_V:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "m3 s-1", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "m3", MAXSTRING);
            }
            break;
        // rate in m/s
        case OUT_AERO_COND:
        case OUT_AERO_COND1:
        case OUT_AERO_COND2:
        case OUT_SURF_COND:
        case OUT_WIND:
            strncpy(nc_vars[i].nc_units, "m s-1", MAXSTRING);
            break;
        // resistance in s/m
        case OUT_AERO_RESIST:
        case OUT_AERO_RESIST1:
        case OUT_AERO_RESIST2:
            strncpy(nc_vars[i].nc_units, "s m-1", MAXSTRING);
            break;
        // temperature in K
        case OUT_RAD_TEMP:
            strncpy(nc_vars[i].nc_units, "K", MAXSTRING);
            break;
        // temperature in C, changed to K for ALMA_OUTPUT
        case OUT_BARESOILT:
        case OUT_LAKE_ICE_TEMP:
        case OUT_LAKE_SURF_TEMP:
        case OUT_SNOW_PACK_TEMP:
        case OUT_SNOW_SURF_TEMP:
        case OUT_SOIL_TEMP:
        case OUT_SOIL_TNODE:
        case OUT_SOIL_TNODE_WL:
        case OUT_SURF_TEMP:
        case OUT_VEGT:
        case OUT_AIR_TEMP:
        case OUT_SNOW_PACKT_BAND:
        case OUT_SNOW_SURFT_BAND:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "K", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "C", MAXSTRING);
            }
            break;
        // energy flux in W/m2
        case OUT_ADV_SENS:
        case OUT_ADVECTION:
        case OUT_ENERGY_ERROR:
        case OUT_FUSION:
        case OUT_GRND_FLUX:
        case OUT_IN_LONG:
        case OUT_LATENT:
        case OUT_LATENT_SUB:
        case OUT_MELT_ENERGY:
        case OUT_LWNET:
        case OUT_SWNET:
        case OUT_R_NET:
        case OUT_RFRZ_ENERGY:
        case OUT_SENSIBLE:
        case OUT_SNOW_FLUX:
        case OUT_LWDOWN:
        case OUT_PAR:
        case OUT_SWDOWN:
        case OUT_ADV_SENS_BAND:
        case OUT_ADVECTION_BAND:
        case OUT_DELTACC_BAND:
        case OUT_GRND_FLUX_BAND:
        case OUT_IN_LONG_BAND:
        case OUT_LATENT_BAND:
        case OUT_LATENT_SUB_BAND:
        case OUT_MELT_ENERGY_BAND:
        case OUT_SWNET_BAND:
        case OUT_RFRZ_ENERGY_BAND:
        case OUT_SENSIBLE_BAND:
        case OUT_SNOW_FLUX_BAND:
        case OUT_APAR:
            strncpy(nc_vars[i].nc_units, "W m-2", MAXSTRING);
            break;
        // energy flux in W/m2, changed to J/m2 for ALMA_OUTPUT
        case OUT_DELTACC:
        case OUT_DELTAH:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "J m-2", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "W m-2", MAXSTRING);
            }
            break;
        // density in kg m-3
        case OUT_DENSITY:
            strncpy(nc_vars[i].nc_units, "kg m-3", MAXSTRING);
            break;
        // concentration in ppm
        case OUT_CATM:
            strncpy(nc_vars[i].nc_units, "ppm", MAXSTRING);
            break;
        // pressure in kPa, changed to Pa for ALMA_OUTPUT
        case OUT_PRESSURE:
        case OUT_VP:
        case OUT_VPD:
            if (options.ALMA_OUTPUT) {
                strncpy(nc_vars[i].nc_units, "Pa", MAXSTRING);
            }
            else {
                strncpy(nc_vars[i].nc_units, "kPa", MAXSTRING);
            }
            break;
        // mixing ratio
        case OUT_QAIR:
            strncpy(nc_vars[i].nc_units, "kg kg-1", MAXSTRING);
            break;
        // carbon fluxes
        case OUT_GPP:
        case OUT_RAUT:
        case OUT_NPP:
        case OUT_LITTERFALL:
        case OUT_RHET:
        case OUT_NEE:
            strncpy(nc_vars[i].nc_units, "g C m-2 d-1", MAXSTRING);
            break;
        // carbon pools
        case OUT_CLITTER:
        case OUT_CINTER:
        case OUT_CSLOW:
            strncpy(nc_vars[i].nc_units, "g C m-2", MAXSTRING);
            break;
        default:
            log_warn("%zd - not defined\n", i);
            break;
        }
    }
}
