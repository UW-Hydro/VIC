/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine sets the metadata structure for VIC output variables
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
 * @brief    Set output met data information
 *****************************************************************************/
void
set_output_met_data_info()
{
    size_t                 v;

    extern option_struct   options;
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    // Build the list of supported output variables

    // Set missing and/or default values
    for (v = 0; v < N_OUTVAR_TYPES; v++) {
        // Set default string values
        strcpy(out_metadata[v].varname, MISSING_S);
        strcpy(out_metadata[v].long_name, MISSING_S);
        strcpy(out_metadata[v].standard_name, MISSING_S);
        strcpy(out_metadata[v].units, MISSING_S);
        strcpy(out_metadata[v].description, MISSING_S);
        // Set default number of elements
        out_metadata[v].nelem = 1;
    }

    // Water Balance Terms - state variables
    /* saturated area fraction */
    strcpy(out_metadata[OUT_ASAT].varname, "OUT_ASAT");
    strcpy(out_metadata[OUT_ASAT].long_name, "asat");
    strcpy(out_metadata[OUT_ASAT].standard_name, "saturated_area_fraction");
    strcpy(out_metadata[OUT_ASAT].units, "1");
    strcpy(out_metadata[OUT_ASAT].description, "saturated area fraction");

    /* lake surface area as fraction of grid cell area [fraction] */
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].varname, "OUT_LAKE_AREA_FRAC");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].long_name, "lake_area_frac");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].standard_name,
           "lake_area_fraction");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].units, "1");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].description,
           "lake surface area as fraction of grid cell area");

    /* lake depth [m] */
    strcpy(out_metadata[OUT_LAKE_DEPTH].varname, "OUT_LAKE_DEPTH");
    strcpy(out_metadata[OUT_LAKE_DEPTH].long_name, "lake_depth");
    strcpy(out_metadata[OUT_LAKE_DEPTH].standard_name, "lake_depth");
    strcpy(out_metadata[OUT_LAKE_DEPTH].units, "m");
    strcpy(out_metadata[OUT_LAKE_DEPTH].description, "lake depth");

    /* moisture stored as lake ice [mm] */
    strcpy(out_metadata[OUT_LAKE_ICE].varname, "OUT_LAKE_ICE");
    strcpy(out_metadata[OUT_LAKE_ICE].long_name, "lake_ice");
    strcpy(out_metadata[OUT_LAKE_ICE].standard_name, "lake_ice");
    strcpy(out_metadata[OUT_LAKE_ICE].units, "mm");
    strcpy(out_metadata[OUT_LAKE_ICE].description,
           "moisture stored as lake ice");

    /* fractional coverage of lake ice [fraction] */
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].varname, "OUT_LAKE_ICE_FRACT");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].long_name, "lake_ice_fract");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].standard_name,
           "lake_ice_area_fraction");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].units, "1");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].description,
           "fractional coverage of lake ice");

    /* thickness of lake ice [cm] */
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].varname, "OUT_LAKE_ICE_HEIGHT");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].long_name, "lake_ice_height");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].standard_name, "lake_ice_height");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].units, "cm");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].description,
           "thickness of lake ice");

    /* liquid water stored in lake [mm over lake area?] */
    strcpy(out_metadata[OUT_LAKE_MOIST].varname, "OUT_LAKE_MOIST");
    strcpy(out_metadata[OUT_LAKE_MOIST].long_name, "lake_moist");
    strcpy(out_metadata[OUT_LAKE_MOIST].standard_name, "lake_moisture");
    strcpy(out_metadata[OUT_LAKE_MOIST].units, "mm");
    strcpy(out_metadata[OUT_LAKE_MOIST].description,
           "liquid water stored in lake");

    /* lake surface area [m2] */
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].varname, "OUT_LAKE_SURF_AREA");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].long_name, "lake_surf_area");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].standard_name, "lake_area");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].units, "m2");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].description, "lake surface area");

    /* liquid water equivalent of snow on top of lake ice [m over lake ice] */
    strcpy(out_metadata[OUT_LAKE_SWE].varname, "OUT_LAKE_SWE");
    strcpy(out_metadata[OUT_LAKE_SWE].long_name, "lake_swe");
    strcpy(out_metadata[OUT_LAKE_SWE].standard_name,
           "lwe_thickness_of_snow_on_lake");
    strcpy(out_metadata[OUT_LAKE_SWE].units, "m");
    strcpy(out_metadata[OUT_LAKE_SWE].description,
           "liquid water equivalent of snow on top of lake ice");

    /* river discharge [m3 s-1] */
    strcpy(out_metadata[OUT_DISCHARGE].varname, "OUT_DISCHARGE");
    strcpy(out_metadata[OUT_DISCHARGE].long_name,
           "water_volume_transport_in_river_channel");
    strcpy(out_metadata[OUT_DISCHARGE].standard_name, "river_discharge");
    strcpy(out_metadata[OUT_DISCHARGE].units, "m3 s-1");
    strcpy(
        out_metadata[OUT_DISCHARGE].description,
        "The water flux or volume transport in rivers is the amount of water flowing in the river channel and flood plain. 'Water' means water in all phases");

    /* volumetric liquid water equivalent of snow on top of lake ice [m3] */
    strcpy(out_metadata[OUT_LAKE_SWE_V].varname, "OUT_LAKE_SWE_V");
    strcpy(out_metadata[OUT_LAKE_SWE_V].long_name, "lake_swe_v");
    strcpy(out_metadata[OUT_LAKE_SWE_V].standard_name,
           "lwe_volume_of_snow_on_lake");
    strcpy(out_metadata[OUT_LAKE_SWE_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_SWE_V].description,
           "liquid water equivalent of snow on top of lake ice");

    /* lake volume [m3] */
    strcpy(out_metadata[OUT_LAKE_VOLUME].varname, "OUT_LAKE_VOLUME");
    strcpy(out_metadata[OUT_LAKE_VOLUME].long_name, "lake_volume");
    strcpy(out_metadata[OUT_LAKE_VOLUME].standard_name, "lake_volume");
    strcpy(out_metadata[OUT_LAKE_VOLUME].units, "m3");
    strcpy(out_metadata[OUT_LAKE_VOLUME].description, "lake volume");

    /* root zone soil moisture [mm] */
    strcpy(out_metadata[OUT_ROOTMOIST].varname, "OUT_ROOTMOIST");
    strcpy(out_metadata[OUT_ROOTMOIST].long_name, "rootmoist");
    strcpy(out_metadata[OUT_ROOTMOIST].standard_name,
           "soil_moisture_in_root_zone");
    strcpy(out_metadata[OUT_ROOTMOIST].units, "mm");
    strcpy(out_metadata[OUT_ROOTMOIST].description, "root zone soil moisture");

    /* fraction of soil moisture (by mass) that is ice, for each soil layer */
    strcpy(out_metadata[OUT_SMFROZFRAC].varname, "OUT_SMFROZFRAC");
    strcpy(out_metadata[OUT_SMFROZFRAC].long_name, "smfrozfrac");
    strcpy(out_metadata[OUT_SMFROZFRAC].standard_name,
           "soil_moisture_ice_fraction");
    strcpy(out_metadata[OUT_SMFROZFRAC].units, "1");
    strcpy(out_metadata[OUT_SMFROZFRAC].description,
           "fraction of soil moisture (by mass) that is ice, for each soil layer");

    /* fraction of soil moisture (by mass) that is liquid, for each soil layer */
    strcpy(out_metadata[OUT_SMLIQFRAC].varname, "OUT_SMLIQFRAC");
    strcpy(out_metadata[OUT_SMLIQFRAC].long_name, "smliqfrac");
    strcpy(out_metadata[OUT_SMLIQFRAC].standard_name,
           "soil_moisture_liquid_fraction");
    strcpy(out_metadata[OUT_SMLIQFRAC].units, "1");
    strcpy(out_metadata[OUT_SMLIQFRAC].description,
           "fraction of soil moisture (by mass) that is liquid, for each soil layer");

    /* snow interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_SNOW_CANOPY].varname, "OUT_SNOW_CANOPY");
    strcpy(out_metadata[OUT_SNOW_CANOPY].long_name, "snow_canopy");
    strcpy(out_metadata[OUT_SNOW_CANOPY].standard_name,
           "lwe_thickness_of_canopy_intercepted_snow");
    strcpy(out_metadata[OUT_SNOW_CANOPY].units, "mm");
    strcpy(out_metadata[OUT_SNOW_CANOPY].description,
           "snow interception storage in canopy");

    /* fractional area of snow cover [fraction] */
    strcpy(out_metadata[OUT_SNOW_COVER].varname, "OUT_SNOW_COVER");
    strcpy(out_metadata[OUT_SNOW_COVER].long_name, "snow_cover");
    strcpy(out_metadata[OUT_SNOW_COVER].standard_name,
           "snow_cover_area_fraction");
    strcpy(out_metadata[OUT_SNOW_COVER].units, "1");
    strcpy(out_metadata[OUT_SNOW_COVER].description,
           "fractional area of snow cover");

    /* depth of snow pack [cm] */
    strcpy(out_metadata[OUT_SNOW_DEPTH].varname, "OUT_SNOW_DEPTH");
    strcpy(out_metadata[OUT_SNOW_DEPTH].long_name, "snow_depth");
    strcpy(out_metadata[OUT_SNOW_DEPTH].standard_name, "thickness_of_snow");
    strcpy(out_metadata[OUT_SNOW_DEPTH].units, "cm");
    strcpy(out_metadata[OUT_SNOW_DEPTH].description, "depth of snow pack");

    /* soil ice content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_ICE].varname, "OUT_SOIL_ICE");
    strcpy(out_metadata[OUT_SOIL_ICE].long_name, "soil_ice");
    strcpy(out_metadata[OUT_SOIL_ICE].standard_name, "soil_moisture_ice_depth");
    strcpy(out_metadata[OUT_SOIL_ICE].units, "mm");
    strcpy(out_metadata[OUT_SOIL_ICE].description,
           "soil ice content for each soil layer");

    /* soil liquid moisture content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_LIQ].varname, "OUT_SOIL_LIQ");
    strcpy(out_metadata[OUT_SOIL_LIQ].long_name, "soil_liq");
    strcpy(out_metadata[OUT_SOIL_LIQ].standard_name,
           "soil_moisture_liquid_depth");
    strcpy(out_metadata[OUT_SOIL_LIQ].units, "mm");
    strcpy(out_metadata[OUT_SOIL_LIQ].description,
           "soil liquid moisture content for each soil layer");

    /* soil ice content [1] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].varname, "OUT_SOIL_ICE_FRAC");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].long_name, "soil_ice_frac");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].standard_name,
           "soil_moisture_ice_depth_fraction");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].units, "1");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].description,
           "soil ice content fraction for each soil layer");

    /* soil liquid moisture content [1] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].varname, "OUT_SOIL_LIQ_FRAC");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].long_name, "soil_liq_frac");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].standard_name,
           "soil_moisture_liquid_depth_fraction");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].units, "1");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].description,
           "soil liquid moisture content fraction for each soil layer");

    /* soil total moisture content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_MOIST].varname, "OUT_SOIL_MOIST");
    strcpy(out_metadata[OUT_SOIL_MOIST].long_name, "soil_moist");
    strcpy(out_metadata[OUT_SOIL_MOIST].standard_name, "soil_moisture");
    strcpy(out_metadata[OUT_SOIL_MOIST].units, "mm");
    strcpy(out_metadata[OUT_SOIL_MOIST].description,
           "soil total moisture content");

    /* vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
    strcpy(out_metadata[OUT_SOIL_WET].varname, "OUT_SOIL_WET");
    strcpy(out_metadata[OUT_SOIL_WET].long_name, "soil_wet");
    strcpy(out_metadata[OUT_SOIL_WET].standard_name,
           "soil_moisture_wetness_fraction");
    strcpy(out_metadata[OUT_SOIL_WET].units, "1");
    strcpy(out_metadata[OUT_SOIL_WET].description,
           "vertical average of (soil moisture - wilting point)/(maximum "
           "soil moisture - wilting point)");

    /* storage of liquid water on surface (ponding) [mm] */
    strcpy(out_metadata[OUT_SURFSTOR].varname, "OUT_SURFSTOR");
    strcpy(out_metadata[OUT_SURFSTOR].long_name, "surfstor");
    strcpy(out_metadata[OUT_SURFSTOR].standard_name,
           "surface_liquid_water_storage");
    strcpy(out_metadata[OUT_SURFSTOR].units, "mm");
    strcpy(out_metadata[OUT_SURFSTOR].description,
           "storage of liquid water on surface (ponding)");

    /* fraction of soil surface that is frozen [fraction] */
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].varname, "OUT_SURF_FROST_FRAC");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].long_name, "surf_frost_frac");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].standard_name,
           "frozen_soil_surface_fraction");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].units, "1");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].description,
           "fraction of soil surface that is frozen");

    /* snow water equivalent in snow pack [mm] */
    strcpy(out_metadata[OUT_SWE].varname, "OUT_SWE");
    strcpy(out_metadata[OUT_SWE].long_name, "swe");
    strcpy(out_metadata[OUT_SWE].standard_name, "lwe_thickness_of_snow");
    strcpy(out_metadata[OUT_SWE].units, "mm");
    strcpy(out_metadata[OUT_SWE].description,
           "snow water equivalent in snow pack");

    /* total moisture interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_WDEW].varname, "OUT_WDEW");
    strcpy(out_metadata[OUT_WDEW].long_name, "wdew");
    strcpy(out_metadata[OUT_WDEW].standard_name,
           "soil_moisture_storage_depth_in_canopy");
    strcpy(out_metadata[OUT_WDEW].units, "mm");
    strcpy(out_metadata[OUT_WDEW].description,
           "total moisture interception storage in canopy");

    /* water table position [cm] (zwt within lowest unsaturated layer) */
    strcpy(out_metadata[OUT_ZWT].varname, "OUT_ZWT");
    strcpy(out_metadata[OUT_ZWT].long_name, "zwt");
    strcpy(out_metadata[OUT_ZWT].standard_name,
           "water_table_position_lowest_layer");
    strcpy(out_metadata[OUT_ZWT].units, "cm");
    strcpy(out_metadata[OUT_ZWT].description,
           "water table position (zwt within lowest unsaturated layer)");

    /* lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
    strcpy(out_metadata[OUT_ZWT_LUMPED].varname, "OUT_ZWT_LUMPED");
    strcpy(out_metadata[OUT_ZWT_LUMPED].long_name, "zwt_lumped");
    strcpy(out_metadata[OUT_ZWT_LUMPED].standard_name,
           "lumped_water_table_position");
    strcpy(out_metadata[OUT_ZWT_LUMPED].units, "cm");
    strcpy(out_metadata[OUT_ZWT_LUMPED].description,
           "lumped water table position (zwt of total moisture across all layers, lumped together)");

    // Water Balance Terms - fluxes
    /* baseflow out of the bottom layer [mm] */
    strcpy(out_metadata[OUT_BASEFLOW].varname, "OUT_BASEFLOW");
    strcpy(out_metadata[OUT_BASEFLOW].long_name, "baseflow");
    strcpy(out_metadata[OUT_BASEFLOW].standard_name, "baseflow_amount");
    strcpy(out_metadata[OUT_BASEFLOW].units, "mm");
    strcpy(out_metadata[OUT_BASEFLOW].description,
           "baseflow out of the bottom layer");

    /* change in canopy interception storage [mm] */
    strcpy(out_metadata[OUT_DELINTERCEPT].varname, "OUT_DELINTERCEPT");
    strcpy(out_metadata[OUT_DELINTERCEPT].long_name, "delintercept");
    strcpy(out_metadata[OUT_DELINTERCEPT].standard_name,
           "change_in_canopy_interception_amount");
    strcpy(out_metadata[OUT_DELINTERCEPT].units, "mm");
    strcpy(out_metadata[OUT_DELINTERCEPT].description,
           "change in canopy interception storage");

    /* change in soil water content [mm] */
    strcpy(out_metadata[OUT_DELSOILMOIST].varname, "OUT_DELSOILMOIST");
    strcpy(out_metadata[OUT_DELSOILMOIST].long_name, "delsoilmoist");
    strcpy(out_metadata[OUT_DELSOILMOIST].standard_name,
           "change_in_soil_water_amount");
    strcpy(out_metadata[OUT_DELSOILMOIST].units, "mm");
    strcpy(out_metadata[OUT_DELSOILMOIST].description,
           "change in soil water content");

    /* change in snow water equivalent [mm] */
    strcpy(out_metadata[OUT_DELSWE].varname, "OUT_DELSWE");
    strcpy(out_metadata[OUT_DELSWE].long_name, "delswe");
    strcpy(out_metadata[OUT_DELSWE].standard_name,
           "change_in_snow_lwe_thickness");
    strcpy(out_metadata[OUT_DELSWE].units, "mm");
    strcpy(out_metadata[OUT_DELSWE].description,
           "change in snow water equivalent");

    /* change in surface liquid water storage  [mm] */
    strcpy(out_metadata[OUT_DELSURFSTOR].varname, "OUT_DELSURFSTOR");
    strcpy(out_metadata[OUT_DELSURFSTOR].long_name, "delsurfstor");
    strcpy(out_metadata[OUT_DELSURFSTOR].standard_name,
           "change_in_surface_liquid_water_storage");
    strcpy(out_metadata[OUT_DELSURFSTOR].units, "mm");
    strcpy(out_metadata[OUT_DELSURFSTOR].description,
           "change in surface liquid water storage");

    /* total net evaporation [mm] */
    strcpy(out_metadata[OUT_EVAP].varname, "OUT_EVAP");
    strcpy(out_metadata[OUT_EVAP].long_name, "evap");
    strcpy(out_metadata[OUT_EVAP].standard_name, "water_evaporation_flux_net");
    strcpy(out_metadata[OUT_EVAP].units, "mm");
    strcpy(out_metadata[OUT_EVAP].description, "total net evaporation");

    /* net evaporation from bare soil [mm] */
    strcpy(out_metadata[OUT_EVAP_BARE].varname, "OUT_EVAP_BARE");
    strcpy(out_metadata[OUT_EVAP_BARE].long_name, "evap_bare");
    strcpy(out_metadata[OUT_EVAP_BARE].standard_name,
           "water_evaporation_from_soil");
    strcpy(out_metadata[OUT_EVAP_BARE].units, "mm");
    strcpy(out_metadata[OUT_EVAP_BARE].description,
           "net evaporation from bare soil");

    /* net evaporation from canopy interception [mm] */
    strcpy(out_metadata[OUT_EVAP_CANOP].varname, "OUT_EVAP_CANOP");
    strcpy(out_metadata[OUT_EVAP_CANOP].long_name, "evap_canop");
    strcpy(out_metadata[OUT_EVAP_CANOP].standard_name,
           "water_evaporation_from_canopy");
    strcpy(out_metadata[OUT_EVAP_CANOP].units, "mm");
    strcpy(out_metadata[OUT_EVAP_CANOP].description,
           "net evaporation from canopy interception");

    /* moisture that reaches top of soil column [mm] */
    strcpy(out_metadata[OUT_INFLOW].varname, "OUT_INFLOW");
    strcpy(out_metadata[OUT_INFLOW].long_name, "inflow");
    strcpy(out_metadata[OUT_INFLOW].standard_name, "soil_column_inflow");
    strcpy(out_metadata[OUT_INFLOW].units, "mm");
    strcpy(out_metadata[OUT_INFLOW].description,
           "moisture that reaches top of soil column");

    /* incoming baseflow from lake catchment [mm] */
    strcpy(out_metadata[OUT_LAKE_BF_IN].varname, "OUT_LAKE_BF_IN");
    strcpy(out_metadata[OUT_LAKE_BF_IN].long_name, "lake_bf_in");
    strcpy(out_metadata[OUT_LAKE_BF_IN].standard_name, "infiltration_amount");
    strcpy(out_metadata[OUT_LAKE_BF_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_BF_IN].description,
           "incoming baseflow from lake catchment");

    /* incoming volumetric baseflow from lake catchment [m3] */
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].varname, "OUT_LAKE_BF_IN_V");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].long_name, "lake_bf_in_v");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].standard_name,
           "baseflow_volume_from_lake_catchment");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].description,
           "incoming volumetric baseflow from lake catchment");

    /* outgoing baseflow lake [mm] */
    strcpy(out_metadata[OUT_LAKE_BF_OUT].varname, "OUT_LAKE_BF_OUT");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].long_name, "lake_bf_out");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].standard_name,
           "baseflow_outgoing_lake");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].units, "mm");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].description,
           "outgoing baseflow lake");

    /* outgoing volumetric baseflow from lake [m3] */
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].varname, "OUT_LAKE_BF_OUT_V");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].long_name, "lake_bf_out_v");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].standard_name,
           "baseflow_outgoing_volume_lake");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].description,
           "outgoing volumetric baseflow from lake");

    /* channel inflow into lake [mm] */
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].varname, "OUT_LAKE_CHAN_IN");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].long_name, "lake_chan_in");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].standard_name,
           "channel_inflow_into_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].description,
           "channel inflow into lake");

    /* volumetric channel inflow into lake [m3] */
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].varname, "OUT_LAKE_CHAN_IN_V");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].long_name, "lake_chan_in_v");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].standard_name,
           "channel_inflow_volume_into_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].description,
           "volumetric channel inflow into lake");

    /* channel outflow from lake [mm] */
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].varname, "OUT_LAKE_CHAN_OUT");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].long_name, "lake_chan_out");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].standard_name,
           "channel_outflow_from_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].units, "mm");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].description,
           "channel outflow from lake");

    /* volumetric channel outflow from lake [m3] */
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].varname, "OUT_LAKE_CHAN_OUT_V");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].long_name, "lake_chan_out_v");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].standard_name,
           "channel_outflow_volume_from_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].description,
           "volumetric channel outflow from lake");

    /* change in lake moisture storage (liquid plus ice cover) [mm] */
    strcpy(out_metadata[OUT_LAKE_DSTOR].varname, "OUT_LAKE_DSTOR");
    strcpy(out_metadata[OUT_LAKE_DSTOR].long_name, "lake_dstor");
    strcpy(out_metadata[OUT_LAKE_DSTOR].standard_name,
           "change_in_lake_moisture_amount");
    strcpy(out_metadata[OUT_LAKE_DSTOR].units, "mm");
    strcpy(out_metadata[OUT_LAKE_DSTOR].description,
           "change in lake moisture storage (liquid plus ice cover)");

    /* volumetric change in lake moisture storage (liquid plus ice cover) [m3] */
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].varname, "OUT_LAKE_DSTOR_V");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].long_name, "lake_dstor_v");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].standard_name,
           "change_in_lake_moisture_volume");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].description,
           "volumetric change in lake moisture storage (liquid plus ice cover)");

    /* change in snowpack on top of lake ice [mm] */
    strcpy(out_metadata[OUT_LAKE_DSWE].varname, "OUT_LAKE_DSWE");
    strcpy(out_metadata[OUT_LAKE_DSWE].long_name, "lake_dswe");
    strcpy(out_metadata[OUT_LAKE_DSWE].standard_name,
           "change_in_snow_lwe_thickness_on_lake_ice");
    strcpy(out_metadata[OUT_LAKE_DSWE].units, "mm");
    strcpy(out_metadata[OUT_LAKE_DSWE].description,
           "change in snowpack on top of lake ice");

    /* volumetric change in snowpack on top of lake ice [m3] */
    strcpy(out_metadata[OUT_LAKE_DSWE_V].varname, "OUT_LAKE_DSWE_V");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].long_name, "lake_dswe_v");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].standard_name,
           "change_in_snow_lwe_volume_on_lake_ice");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].description,
           "volumetric change in snowpack on top of lake ice");

    /* net evaporation from lake surface [mm] */
    strcpy(out_metadata[OUT_LAKE_EVAP].varname, "OUT_LAKE_EVAP");
    strcpy(out_metadata[OUT_LAKE_EVAP].long_name, "lake_evap");
    strcpy(out_metadata[OUT_LAKE_EVAP].standard_name,
           "water_evaporation_from_lake");
    strcpy(out_metadata[OUT_LAKE_EVAP].units, "mm");
    strcpy(out_metadata[OUT_LAKE_EVAP].description,
           "net evaporation from lake surface");

    /* net volumetric evaporation from lake surface [m3] */
    strcpy(out_metadata[OUT_LAKE_EVAP_V].varname, "OUT_LAKE_EVAP_V");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].long_name, "lake_evap_v");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].standard_name,
           "water_evaporation_volume_from_lake");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].description,
           "net volumetric evaporation from lake surface");

    /* volumetric precipitation over lake surface [m3] */
    strcpy(out_metadata[OUT_LAKE_PREC_V].varname, "OUT_LAKE_PREC_V");
    strcpy(out_metadata[OUT_LAKE_PREC_V].long_name, "lake_prec_v");
    strcpy(out_metadata[OUT_LAKE_PREC_V].standard_name,
           "precipitation_over_lake_volume");
    strcpy(out_metadata[OUT_LAKE_PREC_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_PREC_V].description,
           "volumetric precipitation over lake surface");

    /* recharge from lake to surrounding wetland [mm] */
    strcpy(out_metadata[OUT_LAKE_RCHRG].varname, "OUT_LAKE_RCHRG");
    strcpy(out_metadata[OUT_LAKE_RCHRG].long_name, "lake_rchrg");
    strcpy(out_metadata[OUT_LAKE_RCHRG].standard_name,
           "recharge_from_lake_to_wetland");
    strcpy(out_metadata[OUT_LAKE_RCHRG].units, "mm");
    strcpy(out_metadata[OUT_LAKE_RCHRG].description,
           "recharge from lake to surrounding wetland");

    /* volumetric recharge from lake to surrounding wetland [m3] */
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].varname, "OUT_LAKE_RCHRG_V");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].long_name, "lake_rchrg_v");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].standard_name, "");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].description,
           "volumetric recharge from lake to surrounding wetland");

    /* incoming runoff from lake catchment [mm] */
    strcpy(out_metadata[OUT_LAKE_RO_IN].varname, "OUT_LAKE_RO_IN");
    strcpy(out_metadata[OUT_LAKE_RO_IN].long_name, "lake_ro_in");
    strcpy(out_metadata[OUT_LAKE_RO_IN].standard_name,
           "recharge_volume_from_lake_to_wetland");
    strcpy(out_metadata[OUT_LAKE_RO_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_RO_IN].description,
           "incoming runoff from lake catchment");

    /* incoming volumetric runoff from lake catchment [m3] */
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].varname, "OUT_LAKE_RO_IN_V");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].long_name, "lake_ro_in_v");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].standard_name,
           "runoff_volume_from_lake_catchment");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].description,
           "incoming volumetric runoff from lake catchment");

    /* sublimation from lake snow pack [mm] */
    strcpy(out_metadata[OUT_LAKE_VAPFLX].varname, "OUT_LAKE_VAPFLX");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].long_name, "lake_vapflx");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].standard_name,
           "water_sublimation_flux_from_lake_snow");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].units, "mm");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].description,
           "sublimation from lake snow pack");

    /* volumetric sublimation from lake snow pack [m3] */
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].varname, "OUT_LAKE_VAPFLX_V");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].long_name, "lake_vapflx_v");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].standard_name,
           "water_sublimation_flux_volume_from_lake_snow");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].description,
           "volumetric sublimation from lake snow pack");

    /* Potential evapotranspiration (= area-weighted sum of potential
       transpiration and potential soil evaporation). [mm] */
    strcpy(out_metadata[OUT_PET].varname, "OUT_PET");
    strcpy(out_metadata[OUT_PET].long_name, "pet");
    strcpy(out_metadata[OUT_PET].standard_name,
           "water_potential_evaporation_amount");
    strcpy(out_metadata[OUT_PET].units, "mm");
    strcpy(out_metadata[OUT_PET].description,
           "Potential evapotranspiration (= area-weighted sum of potential "
           "transpiration and potential soil evaporation)");

    /* incoming precipitation [mm] */
    strcpy(out_metadata[OUT_PREC].varname, "OUT_PREC");
    strcpy(out_metadata[OUT_PREC].long_name, "prec");
    strcpy(out_metadata[OUT_PREC].standard_name, "precipitation_amount");
    strcpy(out_metadata[OUT_PREC].units, "mm");
    strcpy(out_metadata[OUT_PREC].description, "incoming precipitation");

    /* rainfall [mm] */
    strcpy(out_metadata[OUT_RAINF].varname, "OUT_RAINF");
    strcpy(out_metadata[OUT_RAINF].long_name, "rainf");
    strcpy(out_metadata[OUT_RAINF].standard_name, "rainfall_amount");
    strcpy(out_metadata[OUT_RAINF].units, "mm");
    strcpy(out_metadata[OUT_RAINF].description, "liquid rainfall amount");

    /* refreezing of water in the snow [mm] */
    strcpy(out_metadata[OUT_REFREEZE].varname, "OUT_REFREEZE");
    strcpy(out_metadata[OUT_REFREEZE].long_name, "refreeze");
    strcpy(out_metadata[OUT_REFREEZE].standard_name,
           "lwe_thickness_of_refreezing_water_in_snow");
    strcpy(out_metadata[OUT_REFREEZE].units, "mm");
    strcpy(out_metadata[OUT_REFREEZE].description,
           "refreezing of water in the snow");

    /* surface runoff [mm] */
    strcpy(out_metadata[OUT_RUNOFF].varname, "OUT_RUNOFF");
    strcpy(out_metadata[OUT_RUNOFF].long_name, "runoff");
    strcpy(out_metadata[OUT_RUNOFF].standard_name, "runoff_amount");
    strcpy(out_metadata[OUT_RUNOFF].units, "mm");
    strcpy(out_metadata[OUT_RUNOFF].description, "surface runoff");

    /* snow melt [mm] */
    strcpy(out_metadata[OUT_SNOW_MELT].varname, "OUT_SNOW_MELT");
    strcpy(out_metadata[OUT_SNOW_MELT].long_name, "snow_melt");
    strcpy(out_metadata[OUT_SNOW_MELT].standard_name, "snow_melt_amount");
    strcpy(out_metadata[OUT_SNOW_MELT].units, "mm");
    strcpy(out_metadata[OUT_SNOW_MELT].description, "snow melt");

    /* snowfall [mm] */
    strcpy(out_metadata[OUT_SNOWF].varname, "OUT_SNOWF");
    strcpy(out_metadata[OUT_SNOWF].long_name, "snowf");
    strcpy(out_metadata[OUT_SNOWF].standard_name, "snowfall_lwe_amount");
    strcpy(out_metadata[OUT_SNOWF].units, "mm");
    strcpy(out_metadata[OUT_SNOWF].description, "snowfall");

    /* net sublimation of blowing snow [mm] */
    strcpy(out_metadata[OUT_SUB_BLOWING].varname, "OUT_SUB_BLOWING");
    strcpy(out_metadata[OUT_SUB_BLOWING].long_name, "sub_blowing");
    strcpy(out_metadata[OUT_SUB_BLOWING].standard_name,
           "submlimation_amount_from_blowing_snow");
    strcpy(out_metadata[OUT_SUB_BLOWING].units, "mm");
    strcpy(out_metadata[OUT_SUB_BLOWING].description,
           "net sublimation of blowing snow");

    /* net sublimation from snow stored in canopy [mm] */
    strcpy(out_metadata[OUT_SUB_CANOP].varname, "OUT_SUB_CANOP");
    strcpy(out_metadata[OUT_SUB_CANOP].long_name, "sub_canop");
    strcpy(out_metadata[OUT_SUB_CANOP].standard_name,
           "sublimation_amount_from_canopy_snow");
    strcpy(out_metadata[OUT_SUB_CANOP].units, "mm");
    strcpy(out_metadata[OUT_SUB_CANOP].description,
           "net sublimation from snow stored in canopy");

    /* net sublimation from snow pack (surface and blowing) [mm] */
    strcpy(out_metadata[OUT_SUB_SNOW].varname, "OUT_SUB_SNOW");
    strcpy(out_metadata[OUT_SUB_SNOW].long_name, "sub_snow");
    strcpy(out_metadata[OUT_SUB_SNOW].standard_name,
           "sublimation_amount_from_snow_pack");
    strcpy(out_metadata[OUT_SUB_SNOW].units, "mm");
    strcpy(out_metadata[OUT_SUB_SNOW].description,
           "net sublimation from snow pack (surface and blowing)");

    /* net sublimation from snow pack surface [mm] */
    strcpy(out_metadata[OUT_SUB_SURFACE].varname, "OUT_SUB_SURFACE");
    strcpy(out_metadata[OUT_SUB_SURFACE].long_name, "sub_surface");
    strcpy(out_metadata[OUT_SUB_SURFACE].standard_name,
           "sublimation_amount_from_snow_pack_surface");
    strcpy(out_metadata[OUT_SUB_SURFACE].units, "mm");
    strcpy(out_metadata[OUT_SUB_SURFACE].description,
           "net sublimation from snow pack surface");

    /* net transpiration from vegetation [mm] */
    strcpy(out_metadata[OUT_TRANSP_VEG].varname, "OUT_TRANSP_VEG");
    strcpy(out_metadata[OUT_TRANSP_VEG].long_name, "transp_veg");
    strcpy(out_metadata[OUT_TRANSP_VEG].standard_name, "transpiration_amount");
    strcpy(out_metadata[OUT_TRANSP_VEG].units, "mm");
    strcpy(out_metadata[OUT_TRANSP_VEG].description,
           "net transpiration from vegetation");

    // Energy Balance Terms - state variables
    /* albedo [fraction] */
    strcpy(out_metadata[OUT_ALBEDO].varname, "OUT_ALBEDO");
    strcpy(out_metadata[OUT_ALBEDO].long_name, "albedo");
    strcpy(out_metadata[OUT_ALBEDO].standard_name, "surface_albedo");
    strcpy(out_metadata[OUT_ALBEDO].units, "1");
    strcpy(out_metadata[OUT_ALBEDO].description, "albedo");

    /* bare soil surface temperature [C] */
    strcpy(out_metadata[OUT_BARESOILT].varname, "OUT_BARESOILT");
    strcpy(out_metadata[OUT_BARESOILT].long_name, "baresoilt");
    strcpy(out_metadata[OUT_BARESOILT].standard_name,
           "surface_temperature_bare_soil");
    strcpy(out_metadata[OUT_BARESOILT].units, "C");
    strcpy(out_metadata[OUT_BARESOILT].description,
           "bare soil surface temperature");

    /* depth of freezing fronts [cm] for each freezing front */
    strcpy(out_metadata[OUT_FDEPTH].varname, "OUT_FDEPTH");
    strcpy(out_metadata[OUT_FDEPTH].long_name, "fdepth");
    strcpy(out_metadata[OUT_FDEPTH].standard_name, "freezing_fronts_depth");
    strcpy(out_metadata[OUT_FDEPTH].units, "cm");
    strcpy(out_metadata[OUT_FDEPTH].description,
           "depth of freezing fronts for each freezing front");

    /* lake ice temperature [K] */
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].varname, "OUT_LAKE_ICE_TEMP");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].long_name, "lake_ice_temp");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].standard_name,
           "lake_ice_temperature");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].units, "K");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].description, "lake ice temperature");

    /* lake surface temperature [K] */
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].varname, "OUT_LAKE_SURF_TEMP");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].long_name, "lake_surf_temp");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].standard_name,
           "lake_surface_temperature");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].units, "K");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].description,
           "lake surface temperature");

    /* average radiative surface temperature [K] */
    strcpy(out_metadata[OUT_RAD_TEMP].varname, "OUT_RAD_TEMP");
    strcpy(out_metadata[OUT_RAD_TEMP].long_name, "rad_temp");
    strcpy(out_metadata[OUT_RAD_TEMP].standard_name,
           "surface_radiative_temperature");
    strcpy(out_metadata[OUT_RAD_TEMP].units, "K");
    strcpy(out_metadata[OUT_RAD_TEMP].description,
           "average radiative surface temperature");

    /* snow albedo [fraction] */
    strcpy(out_metadata[OUT_SALBEDO].varname, "OUT_SALBEDO");
    strcpy(out_metadata[OUT_SALBEDO].long_name, "salbedo");
    strcpy(out_metadata[OUT_SALBEDO].standard_name, "snow_albedo");
    strcpy(out_metadata[OUT_SALBEDO].units, "1");
    strcpy(out_metadata[OUT_SALBEDO].description, "snow albedo");

    /* snow pack temperature [C] */
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].varname, "OUT_SNOW_PACK_TEMP");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].long_name, "snow_pack_temp");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].standard_name,
           "snow_pack_temperature");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].units, "C");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].description,
           "snow pack temperature");

    /* snow surface temperature [C] */
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].varname, "OUT_SNOW_SURF_TEMP");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].long_name, "snow_surf_temp");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].standard_name,
           "snow_surface_temperature");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].units, "C");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].description,
           "snow surface temperature");

    /* snow surface temperature flag */
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].varname, "OUT_SNOWT_FBFLAG");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].long_name, "snowt_fbflag");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].standard_name,
           "snow_surface_temperature_flag");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].description,
           "snow surface temperature flag");

    /* soil temperature [C] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_TEMP].varname, "OUT_SOIL_TEMP");
    strcpy(out_metadata[OUT_SOIL_TEMP].long_name, "soil_temp");
    strcpy(out_metadata[OUT_SOIL_TEMP].standard_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TEMP].units, "C");
    strcpy(out_metadata[OUT_SOIL_TEMP].description,
           "soil temperature for each soil layer");

    /* soil temperature [C] for each soil thermal node */
    strcpy(out_metadata[OUT_SOIL_TNODE].varname, "OUT_SOIL_TNODE");
    strcpy(out_metadata[OUT_SOIL_TNODE].long_name, "soil_tnode");
    strcpy(out_metadata[OUT_SOIL_TNODE].standard_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TNODE].units, "C");
    strcpy(out_metadata[OUT_SOIL_TNODE].description,
           "soil temperature for each soil thermal node");

    /* soil temperature [C] for each soil thermal node in the wetland */
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].varname, "OUT_SOIL_TNODE_WL");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].long_name, "soil_tnode_wl");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].standard_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].units, "C");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].description,
           "soil temperature for each soil thermal node in the wetland");

    /* soil temperature flag for each soil thermal node */
    strcpy(out_metadata[OUT_SOILT_FBFLAG].varname, "OUT_SOILT_FBFLAG");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].long_name, "soilt_fbflag");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].standard_name,
           "soil_temperature_flag");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].description,
           "soil temperature flag for each soil thermal node");

    /* average surface temperature [C] */
    strcpy(out_metadata[OUT_SURF_TEMP].varname, "OUT_SURF_TEMP");
    strcpy(out_metadata[OUT_SURF_TEMP].long_name, "surf_temp");
    strcpy(out_metadata[OUT_SURF_TEMP].standard_name, "surface_temperature");
    strcpy(out_metadata[OUT_SURF_TEMP].units, "C");
    strcpy(out_metadata[OUT_SURF_TEMP].description,
           "average surface temperature");

    /* surface temperature flag */
    strcpy(out_metadata[OUT_SURFT_FBFLAG].varname, "OUT_SURFT_FBFLAG");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].long_name, "surft_fbflag");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].standard_name,
           "surface_temperature_flag");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].description,
           "surface temperature flag");

    /* Tcanopy flag */
    strcpy(out_metadata[OUT_TCAN_FBFLAG].varname, "OUT_TCAN_FBFLAG");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].long_name, "tcan_fbflag");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].standard_name,
           "canopy_temperature_flag");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].description,
           "Canopy temperature fallback flag");

    /* depth of thawing fronts [cm] for each thawing front */
    strcpy(out_metadata[OUT_TDEPTH].varname, "OUT_TDEPTH");
    strcpy(out_metadata[OUT_TDEPTH].long_name, "tdepth");
    strcpy(out_metadata[OUT_TDEPTH].standard_name, "depth_of_thawing_fronts");
    strcpy(out_metadata[OUT_TDEPTH].units, "cm");
    strcpy(out_metadata[OUT_TDEPTH].description,
           "depth of thawing fronts for each thawing front");

    /* Tfoliage flag */
    strcpy(out_metadata[OUT_TFOL_FBFLAG].varname, "OUT_TFOL_FBFLAG");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].long_name, "tfol_fbflag");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].standard_name,
           "foliage_temperature_flag");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].description,
           "foilage temperature fallback flag");

    /* average vegetation canopy temperature [C] */
    strcpy(out_metadata[OUT_VEGT].varname, "OUT_VEGT");
    strcpy(out_metadata[OUT_VEGT].long_name, "vegt");
    strcpy(out_metadata[OUT_VEGT].standard_name, "canopy_temperature");
    strcpy(out_metadata[OUT_VEGT].units, "C");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].description,
           "average vegetation canopy temperature");

    // Energy Balance Terms - fluxes
    /* net sensible heat advected to snow pack [W m-2] */
    strcpy(out_metadata[OUT_ADV_SENS].varname, "OUT_ADV_SENS");
    strcpy(out_metadata[OUT_ADV_SENS].long_name, "adv_sens");
    strcpy(out_metadata[OUT_ADV_SENS].standard_name,
           "net_sensible_heat_flux_to_snow_pack");
    strcpy(out_metadata[OUT_ADV_SENS].units, "W m-2");
    strcpy(out_metadata[OUT_ADV_SENS].description,
           "net sensible heat advected to snow pack");

    /* advected energy [W m-2] */
    strcpy(out_metadata[OUT_ADVECTION].varname, "OUT_ADVECTION");
    strcpy(out_metadata[OUT_ADVECTION].long_name, "advection");
    strcpy(out_metadata[OUT_ADVECTION].standard_name, "advected_energy");
    strcpy(out_metadata[OUT_ADVECTION].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTION].description, "advected energy ");

    /* rate of change in cold content in snow pack [W m-2] */
    strcpy(out_metadata[OUT_DELTACC].varname, "OUT_DELTACC");
    strcpy(out_metadata[OUT_DELTACC].long_name, "deltacc");
    strcpy(out_metadata[OUT_DELTACC].standard_name,
           "rate_change_in_snow_pack_cold_content");
    strcpy(out_metadata[OUT_DELTACC].units, "W m-2");
    strcpy(out_metadata[OUT_DELTACC].description,
           "rate of change in cold content in snow pack");

    /* rate of change in heat storage [W m-2] */
    strcpy(out_metadata[OUT_DELTAH].varname, "OUT_DELTAH");
    strcpy(out_metadata[OUT_DELTAH].long_name, "deltah");
    strcpy(out_metadata[OUT_DELTAH].standard_name,
           "rate_change_in_heat_storage");
    strcpy(out_metadata[OUT_DELTAH].units, "W m-2");
    strcpy(out_metadata[OUT_DELTAH].description,
           "rate of change in heat storage");

    /* energy budget error [W m-2] */
    strcpy(out_metadata[OUT_ENERGY_ERROR].varname, "OUT_ENERGY_ERROR");
    strcpy(out_metadata[OUT_ENERGY_ERROR].long_name, "energy_error");
    strcpy(out_metadata[OUT_ENERGY_ERROR].standard_name, "energy_budget_error");
    strcpy(out_metadata[OUT_ENERGY_ERROR].units, "W m-2");
    strcpy(out_metadata[OUT_ENERGY_ERROR].description, "energy budget error");

    /* water budget error [mm] */
    strcpy(out_metadata[OUT_WATER_ERROR].varname, "OUT_WATER_ERROR");
    strcpy(out_metadata[OUT_WATER_ERROR].long_name, "water_error");
    strcpy(out_metadata[OUT_WATER_ERROR].standard_name, "water_budget_error");
    strcpy(out_metadata[OUT_WATER_ERROR].units, "mm");
    strcpy(out_metadata[OUT_WATER_ERROR].description, "water budget error");

    /* net energy used to melt/freeze soil moisture [W m-2] */
    strcpy(out_metadata[OUT_FUSION].varname, "OUT_FUSION");
    strcpy(out_metadata[OUT_FUSION].long_name, "fusion");
    strcpy(out_metadata[OUT_FUSION].standard_name,
           "energy_of_fusion_in_soil_moisture");
    strcpy(out_metadata[OUT_FUSION].units, "W m-2");
    strcpy(out_metadata[OUT_FUSION].description,
           "net energy used to melt/freeze soil moisture");

    /* net heat flux into ground [W m-2] */
    strcpy(out_metadata[OUT_GRND_FLUX].varname, "OUT_GRND_FLUX");
    strcpy(out_metadata[OUT_GRND_FLUX].long_name, "grnd_flux");
    strcpy(out_metadata[OUT_GRND_FLUX].standard_name,
           "downward_heat_flux_in_soil");
    strcpy(out_metadata[OUT_GRND_FLUX].units, "W m-2");
    strcpy(out_metadata[OUT_GRND_FLUX].description,
           "net heat flux into ground");

    /* incoming longwave flux at surface (under veg) [W m-2] */
    strcpy(out_metadata[OUT_IN_LONG].varname, "OUT_IN_LONG");
    strcpy(out_metadata[OUT_IN_LONG].long_name, "in_long");
    strcpy(out_metadata[OUT_IN_LONG].standard_name,
           "downwelling_longwave_flux_at_ground_surface");
    strcpy(out_metadata[OUT_IN_LONG].units, "W m-2");
    strcpy(out_metadata[OUT_IN_LONG].description,
           "incoming longwave flux at surface (under veg)");

    /* net upward latent heat flux [W m-2] */
    strcpy(out_metadata[OUT_LATENT].varname, "OUT_LATENT");
    strcpy(out_metadata[OUT_LATENT].long_name, "latent");
    strcpy(out_metadata[OUT_LATENT].standard_name,
           "surface_upward_latent_heat_flux");
    strcpy(out_metadata[OUT_LATENT].units, "W m-2");
    strcpy(out_metadata[OUT_LATENT].description, "net upward latent heat flux");

    /* net upward latent heat flux from sublimation [W m-2] */
    strcpy(out_metadata[OUT_LATENT_SUB].varname, "OUT_LATENT_SUB");
    strcpy(out_metadata[OUT_LATENT_SUB].long_name, "latent_sub");
    strcpy(out_metadata[OUT_LATENT_SUB].standard_name,
           "surface_net_latent_heat_flux_from_sublimation");
    strcpy(out_metadata[OUT_LATENT_SUB].units, "W m-2");
    strcpy(out_metadata[OUT_LATENT_SUB].description,
           "net upward latent heat flux from sublimation");

    /* energy of fusion (melting) [W m-2] */
    strcpy(out_metadata[OUT_MELT_ENERGY].varname, "OUT_MELT_ENERGY");
    strcpy(out_metadata[OUT_MELT_ENERGY].long_name, "melt_energy");
    strcpy(out_metadata[OUT_MELT_ENERGY].standard_name, "energy_of_fusion");
    strcpy(out_metadata[OUT_MELT_ENERGY].units, "W m-2");
    strcpy(out_metadata[OUT_MELT_ENERGY].description,
           "energy of fusion (melting)");

    /* net downward longwave flux [W m-2] */
    strcpy(out_metadata[OUT_LWNET].varname, "OUT_LWNET");
    strcpy(out_metadata[OUT_LWNET].long_name, "lwnet");
    strcpy(out_metadata[OUT_LWNET].standard_name,
           "net_downward_longwave_flux_at_surface");
    strcpy(out_metadata[OUT_LWNET].units, "W m-2");
    strcpy(out_metadata[OUT_LWNET].description, "net downward longwave flux");

    /* net downward shortwave flux [W m-2] */
    strcpy(out_metadata[OUT_SWNET].varname, "OUT_SWNET");
    strcpy(out_metadata[OUT_SWNET].long_name, "swnet");
    strcpy(out_metadata[OUT_SWNET].standard_name,
           "net_downward_shortwave_flux_at_surface");
    strcpy(out_metadata[OUT_SWNET].units, "W m-2");
    strcpy(out_metadata[OUT_SWNET].description, "net downward shortwave flux");

    /* net downward radiation flux [W m-2] */
    strcpy(out_metadata[OUT_R_NET].varname, "OUT_R_NET");
    strcpy(out_metadata[OUT_R_NET].long_name, "r_net");
    strcpy(out_metadata[OUT_R_NET].standard_name,
           "net_downward_radiation_flux_at_surface");
    strcpy(out_metadata[OUT_R_NET].units, "W m-2");
    strcpy(out_metadata[OUT_R_NET].description, "net downward radiation flux");

    /* net energy used to refreeze liquid water in snowpack [W m-2] */
    strcpy(out_metadata[OUT_RFRZ_ENERGY].varname, "OUT_RFRZ_ENERGY");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].long_name, "rfrz_energy");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].standard_name,
           "net_energy_used_to_refreeze_water_in_snowpack");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].units, "W m-2");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].description,
           "net energy used to refreeze liquid water in snowpack");

    /* net upward sensible heat flux [W m-2] */
    strcpy(out_metadata[OUT_SENSIBLE].varname, "OUT_SENSIBLE");
    strcpy(out_metadata[OUT_SENSIBLE].long_name, "sensible");
    strcpy(out_metadata[OUT_SENSIBLE].standard_name,
           "surface_upward_net_sensible_heat_flux");
    strcpy(out_metadata[OUT_SENSIBLE].units, "W m-2");
    strcpy(out_metadata[OUT_SENSIBLE].description,
           "net upward sensible heat flux");

    /* energy flux through snow pack [W m-2] */
    strcpy(out_metadata[OUT_SNOW_FLUX].varname, "OUT_SNOW_FLUX");
    strcpy(out_metadata[OUT_SNOW_FLUX].long_name, "snow_flux");
    strcpy(out_metadata[OUT_SNOW_FLUX].standard_name,
           "energy_flux_through_snow_pack");
    strcpy(out_metadata[OUT_SNOW_FLUX].units, "W m-2");
    strcpy(out_metadata[OUT_SNOW_FLUX].description,
           "energy flux through snow pack");

    // Miscellaneous Terms

    /* "scene" aerodynamic conductance [m/s] (tiles with overstory contribute
        overstory conductance; others contribue surface conductance) */
    strcpy(out_metadata[OUT_AERO_COND].varname, "OUT_AERO_COND");
    strcpy(out_metadata[OUT_AERO_COND].long_name, "aero_cond");
    strcpy(out_metadata[OUT_AERO_COND].standard_name,
           "aerodynamic_conductance");
    strcpy(out_metadata[OUT_AERO_COND].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND].description,
           "scene aerodynamic conductance (tiles with overstory contribute "
           "overstory conductance; others contribue surface conductance)");

    /* surface aerodynamic conductance [m/s] */
    strcpy(out_metadata[OUT_AERO_COND1].varname, "OUT_AERO_COND1");
    strcpy(out_metadata[OUT_AERO_COND1].long_name, "aero_cond1");
    strcpy(out_metadata[OUT_AERO_COND1].standard_name,
           "aerodynamic_conductance_surface");
    strcpy(out_metadata[OUT_AERO_COND1].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND1].description,
           "surface aerodynamic conductance");

    /* overstory aerodynamic conductance [m/s] */
    strcpy(out_metadata[OUT_AERO_COND2].varname, "OUT_AERO_COND2");
    strcpy(out_metadata[OUT_AERO_COND2].long_name, "aero_cond2");
    strcpy(out_metadata[OUT_AERO_COND2].standard_name,
           "aerodynamic_conductance_overstory");
    strcpy(out_metadata[OUT_AERO_COND2].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND2].description,
           "overstory aerodynamic conductance");

    /* "scene" aerodynamic resistance [s m-1] (tiles with overstory contribute overstory resistance; others contribue surface resistance)*/
    strcpy(out_metadata[OUT_AERO_RESIST].varname, "OUT_AERO_RESIST");
    strcpy(out_metadata[OUT_AERO_RESIST].long_name, "aero_resist");
    strcpy(out_metadata[OUT_AERO_RESIST].standard_name,
           "aerodynamic_resistance");
    strcpy(out_metadata[OUT_AERO_RESIST].units, "s m-1");
    strcpy(out_metadata[OUT_AERO_RESIST].description,
           "scene aerodynamic resistance (tiles with overstory contribute overstory resistance; others contribue surface resistance)");

    /* surface aerodynamic resistance [m/s] */
    strcpy(out_metadata[OUT_AERO_RESIST1].varname, "OUT_AERO_RESIST1");
    strcpy(out_metadata[OUT_AERO_RESIST1].long_name, "aero_resist1");
    strcpy(out_metadata[OUT_AERO_RESIST1].standard_name,
           "aerodynamic_resistance_surface");
    strcpy(out_metadata[OUT_AERO_RESIST1].units, "s m-1");
    strcpy(out_metadata[OUT_AERO_RESIST1].description,
           "surface aerodynamic resistance");

    /* overstory aerodynamic resistance [m/s] */
    strcpy(out_metadata[OUT_AERO_RESIST2].varname, "OUT_AERO_RESIST2");
    strcpy(out_metadata[OUT_AERO_RESIST2].long_name, "aero_resist2");
    strcpy(out_metadata[OUT_AERO_RESIST2].standard_name,
           "aerodynamic_resistance_overstory");
    strcpy(out_metadata[OUT_AERO_RESIST2].units, "s m-1");
    strcpy(out_metadata[OUT_AERO_RESIST2].description,
           "overstory aerodynamic resistance");

    /* air temperature [C] */
    strcpy(out_metadata[OUT_AIR_TEMP].varname, "OUT_AIR_TEMP");
    strcpy(out_metadata[OUT_AIR_TEMP].long_name, "air_temp");
    strcpy(out_metadata[OUT_AIR_TEMP].standard_name, "air_temperature");
    strcpy(out_metadata[OUT_AIR_TEMP].units, "C");
    strcpy(out_metadata[OUT_AIR_TEMP].description, "air temperature");

    /* atmospheric CO2 concentration [ppm] */
    strcpy(out_metadata[OUT_CATM].varname, "OUT_CATM");
    strcpy(out_metadata[OUT_CATM].long_name, "catm");
    strcpy(out_metadata[OUT_CATM].standard_name,
           "concentration_of_carbon_dioxide_in_air");
    strcpy(out_metadata[OUT_CATM].units, "ppm");
    strcpy(out_metadata[OUT_CATM].description, "atmospheric CO2 concentration");

    /* near-surface atmospheric density [kg m-3] */
    strcpy(out_metadata[OUT_DENSITY].varname, "OUT_DENSITY");
    strcpy(out_metadata[OUT_DENSITY].long_name, "density");
    strcpy(out_metadata[OUT_DENSITY].standard_name, "air_density");
    strcpy(out_metadata[OUT_DENSITY].units, "kg m-3");
    strcpy(out_metadata[OUT_DENSITY].description,
           "near-surface atmospheric density");

    /* fractional area covered by plant canopy [fraction] */
    strcpy(out_metadata[OUT_FCANOPY].varname, "OUT_FCANOPY");
    strcpy(out_metadata[OUT_FCANOPY].long_name, "fcanopy");
    strcpy(out_metadata[OUT_FCANOPY].standard_name,
           "canopy_cover_area_fraction");
    strcpy(out_metadata[OUT_FCANOPY].units, "1");
    strcpy(out_metadata[OUT_FCANOPY].description,
           "fractional area covered by plant canopy");

    /* fraction of incoming shortwave that is direct [fraction] */
    strcpy(out_metadata[OUT_FDIR].varname, "OUT_FDIR");
    strcpy(out_metadata[OUT_FDIR].long_name, "fdir");
    strcpy(out_metadata[OUT_FDIR].standard_name,
           "fraction_of_incoming_shorwave_radiation_that_is_direct");
    strcpy(out_metadata[OUT_FDIR].units, "1");
    strcpy(out_metadata[OUT_FDIR].description,
           "fraction of incoming shortwave that is direct");

    /* leaf area index [1] */
    strcpy(out_metadata[OUT_LAI].varname, "OUT_LAI");
    strcpy(out_metadata[OUT_LAI].long_name, "lai");
    strcpy(out_metadata[OUT_LAI].standard_name, "leaf_area_index");
    strcpy(out_metadata[OUT_LAI].units, "1");
    strcpy(out_metadata[OUT_LAI].description, "leaf area index");

    /* incoming longwave [W m-2] */
    strcpy(out_metadata[OUT_LWDOWN].varname, "OUT_LWDOWN");
    strcpy(out_metadata[OUT_LWDOWN].long_name, "lwdown");
    strcpy(out_metadata[OUT_LWDOWN].standard_name,
           "downwelling_longwave_flux_in_air");
    strcpy(out_metadata[OUT_LWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_LWDOWN].description, "incoming longwave");

    /* incoming photosynthetically active radiation [W m-2] */
    strcpy(out_metadata[OUT_PAR].varname, "OUT_PAR");
    strcpy(out_metadata[OUT_PAR].long_name, "par");
    strcpy(out_metadata[OUT_PAR].standard_name,
           "surface_downwelling_photosynthetic_radiative_flux_in_air");
    strcpy(out_metadata[OUT_PAR].units, "W m-2");
    strcpy(out_metadata[OUT_PAR].description,
           "incoming photosynthetically active radiation");

    /* near surface atmospheric pressure [kPa] */
    strcpy(out_metadata[OUT_PRESSURE].varname, "OUT_PRESSURE");
    strcpy(out_metadata[OUT_PRESSURE].long_name, "pressure");
    strcpy(out_metadata[OUT_PRESSURE].standard_name, "surface_air_pressure");
    strcpy(out_metadata[OUT_PRESSURE].units, "kPa");
    strcpy(out_metadata[OUT_PRESSURE].description,
           "near surface atmospheric pressure");

    /* specific humidity [kg/kg] */
    strcpy(out_metadata[OUT_QAIR].varname, "OUT_QAIR");
    strcpy(out_metadata[OUT_QAIR].long_name, "qair");
    strcpy(out_metadata[OUT_QAIR].standard_name, "specific_humidity");
    strcpy(out_metadata[OUT_QAIR].units, "1");
    strcpy(out_metadata[OUT_QAIR].description, "specific humidity");

    /* relative humidity [fraction]*/
    strcpy(out_metadata[OUT_REL_HUMID].varname, "OUT_REL_HUMID");
    strcpy(out_metadata[OUT_REL_HUMID].long_name, "rel_humid");
    strcpy(out_metadata[OUT_REL_HUMID].standard_name, "relative_humidity");
    strcpy(out_metadata[OUT_REL_HUMID].units, "1");
    strcpy(out_metadata[OUT_REL_HUMID].description, "relative humidity");

    /* incoming shortwave [W m-2] */
    strcpy(out_metadata[OUT_SWDOWN].varname, "OUT_SWDOWN");
    strcpy(out_metadata[OUT_SWDOWN].long_name, "swdown");
    strcpy(out_metadata[OUT_SWDOWN].standard_name, "incoming shortwave");
    strcpy(out_metadata[OUT_SWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_SWDOWN].description, "incoming shortwave");

    /* surface conductance [m/s] */
    strcpy(out_metadata[OUT_SURF_COND].varname, "OUT_SURF_COND");
    strcpy(out_metadata[OUT_SURF_COND].long_name, "surf_cond");
    strcpy(out_metadata[OUT_SURF_COND].standard_name,
           "surface_conductance");
    strcpy(out_metadata[OUT_SURF_COND].units, "m s-1");
    strcpy(out_metadata[OUT_SURF_COND].description, "surface conductance");

    /* near surface vapor pressure [kPa] */
    strcpy(out_metadata[OUT_VP].varname, "OUT_VP");
    strcpy(out_metadata[OUT_VP].long_name, "vp");
    strcpy(out_metadata[OUT_VP].standard_name, "water_vapor_pressure");
    strcpy(out_metadata[OUT_VP].units, "kPa");
    strcpy(out_metadata[OUT_VP].description, "near surface vapor pressure");

    /* near surface vapor pressure deficit [kPa] */
    strcpy(out_metadata[OUT_VPD].varname, "OUT_VPD");
    strcpy(out_metadata[OUT_VPD].long_name, "vpd");
    strcpy(out_metadata[OUT_VPD].standard_name,
           "water_vapor_saturation_deficit");
    strcpy(out_metadata[OUT_VPD].units, "kPa");
    strcpy(out_metadata[OUT_VPD].description,
           "near surface vapor pressure deficit");

    /* near surface wind speed [m/s] */
    strcpy(out_metadata[OUT_WIND].varname, "OUT_WIND");
    strcpy(out_metadata[OUT_WIND].long_name, "wind");
    strcpy(out_metadata[OUT_WIND].standard_name, "wind_speed");
    strcpy(out_metadata[OUT_WIND].units, "m s-1");
    strcpy(out_metadata[OUT_WIND].description, "near surface wind speed");

    // Carbon-cycling Terms
    /* absorbed PAR [W m-2] */
    strcpy(out_metadata[OUT_APAR].varname, "OUT_APAR");
    strcpy(out_metadata[OUT_APAR].long_name, "apar");
    strcpy(out_metadata[OUT_APAR].standard_name,
           "absorbed_surface_diffuse_downwelling_photosynthetic_radiative_flux");
    strcpy(out_metadata[OUT_APAR].units, "W m-2");
    strcpy(out_metadata[OUT_APAR].description, "absorbed PAR");

    /* gross primary productivity [g C/m2s] */
    strcpy(out_metadata[OUT_GPP].varname, "OUT_GPP");
    strcpy(out_metadata[OUT_GPP].long_name, "gpp");
    strcpy(out_metadata[OUT_GPP].standard_name,
           "gross_primary_productivity_of_biomass_expressed_as_carbon");
    strcpy(out_metadata[OUT_GPP].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_GPP].description, "gross primary productivity");

    /* autotrophic respiration [g C/m2s] */
    strcpy(out_metadata[OUT_RAUT].varname, "OUT_RAUT");
    strcpy(out_metadata[OUT_RAUT].long_name, "raut");
    strcpy(out_metadata[OUT_RAUT].standard_name,
           "autotrophic_respiration_carbon_flux");
    strcpy(out_metadata[OUT_RAUT].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_RAUT].description, "autotrophic respiration");

    /* net primary productivity [g C/m2s] */
    strcpy(out_metadata[OUT_NPP].varname, "OUT_NPP");
    strcpy(out_metadata[OUT_NPP].long_name, "npp");
    strcpy(out_metadata[OUT_NPP].standard_name,
           "net_primary_productivity_of_biomass_expressed_as_carbon");
    strcpy(out_metadata[OUT_NPP].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_NPP].description, "et primary productivity");

    /* flux of carbon from living biomass into soil [g C/m2d] */
    strcpy(out_metadata[OUT_LITTERFALL].varname, "OUT_LITTERFALL");
    strcpy(out_metadata[OUT_LITTERFALL].long_name, "litterfall");
    strcpy(out_metadata[OUT_LITTERFALL].standard_name,
           "carbon_mass_flux_into_soil_from_litter");
    strcpy(out_metadata[OUT_LITTERFALL].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_LITTERFALL].description,
           "flux of carbon from living biomass into soil");

    /* heterotrophic respiration [g C/m2d] */
    strcpy(out_metadata[OUT_RHET].varname, "OUT_RHET");
    strcpy(out_metadata[OUT_RHET].long_name, "rhet");
    strcpy(out_metadata[OUT_RHET].standard_name, "heterotrophic_respiration");
    strcpy(out_metadata[OUT_RHET].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_RHET].description, "heterotrophic respiration");

    /* net ecosystem exchange [g C/m2d] */
    strcpy(out_metadata[OUT_NEE].varname, "OUT_NEE");
    strcpy(out_metadata[OUT_NEE].long_name, "nee");
    strcpy(out_metadata[OUT_NEE].standard_name,
           "net_ecosystem_exhanged_expressed_as_carbon");
    strcpy(out_metadata[OUT_NEE].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_NEE].description, "net ecosystem exchange");

    /* litter pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CLITTER].varname, "OUT_CLITTER");
    strcpy(out_metadata[OUT_CLITTER].long_name, "clitter");
    strcpy(out_metadata[OUT_CLITTER].standard_name, "litter_carbon_content");
    strcpy(out_metadata[OUT_CLITTER].units, "g m-2");
    strcpy(out_metadata[OUT_CLITTER].description, "litter pool carbon density");

    /* intermediate pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CINTER].varname, "OUT_CINTER");
    strcpy(out_metadata[OUT_CINTER].long_name, "cinter");
    strcpy(out_metadata[OUT_CINTER].standard_name,
           "intermediate_pool_carbon_content");
    strcpy(out_metadata[OUT_CINTER].units, "g m-2");
    strcpy(out_metadata[OUT_CINTER].description,
           "intermediate pool carbon density");

    /* slow pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CSLOW].varname, "OUT_CSLOW");
    strcpy(out_metadata[OUT_CSLOW].long_name, "cslow");
    strcpy(out_metadata[OUT_CSLOW].standard_name, "slow_pool_carbon_content");
    strcpy(out_metadata[OUT_CSLOW].units, "g m-2");
    strcpy(out_metadata[OUT_CSLOW].description, "slow pool carbon density");

    // Band-specific quantities
    /* net sensible heat flux advected to snow pack [W m-2] */
    strcpy(out_metadata[OUT_ADV_SENS_BAND].varname, "OUT_ADV_SENS_BAND");
    strcpy(out_metadata[OUT_ADV_SENS_BAND].long_name, "adv_sens_band");
    strcpy(out_metadata[OUT_ADV_SENS_BAND].standard_name,
           out_metadata[OUT_ADV_SENS].standard_name);
    strcpy(out_metadata[OUT_ADV_SENS_BAND].units,
           out_metadata[OUT_ADV_SENS].units);
    strcpy(out_metadata[OUT_ADV_SENS_BAND].description,
           out_metadata[OUT_ADV_SENS].description);

    strcpy(out_metadata[OUT_ADV_SENS].varname, "OUT_ADV_SENS");
    strcpy(out_metadata[OUT_ADV_SENS].long_name, "adv_sens");
    strcpy(out_metadata[OUT_ADV_SENS].standard_name,
           "net_sensible_heat_flux_to_snow_pack");
    strcpy(out_metadata[OUT_ADV_SENS].units, "W m-2");
    strcpy(out_metadata[OUT_ADV_SENS].description,
           "net sensible heat advected to snow pack");

    /* advected energy [W m-2] */
    strcpy(out_metadata[OUT_ADVECTION_BAND].varname, "OUT_ADVECTION_BAND");
    strcpy(out_metadata[OUT_ADVECTION_BAND].long_name, "advection_band");
    strcpy(out_metadata[OUT_ADVECTION_BAND].standard_name,
           out_metadata[OUT_ADVECTION].standard_name);
    strcpy(out_metadata[OUT_ADVECTION_BAND].units,
           out_metadata[OUT_ADVECTION].units);
    strcpy(out_metadata[OUT_ADVECTION_BAND].description,
           out_metadata[OUT_ADVECTION].description);

    /* albedo [fraction] */
    strcpy(out_metadata[OUT_ALBEDO_BAND].varname, "OUT_ALBEDO_BAND");
    strcpy(out_metadata[OUT_ALBEDO_BAND].long_name, "albedo_band");
    strcpy(out_metadata[OUT_ALBEDO_BAND].standard_name,
           out_metadata[OUT_ALBEDO].standard_name);
    strcpy(out_metadata[OUT_ALBEDO_BAND].units, out_metadata[OUT_ALBEDO].units);
    strcpy(out_metadata[OUT_ALBEDO_BAND].description,
           out_metadata[OUT_ALBEDO].description);

    /* change in cold content in snow pack [W m-2] */
    strcpy(out_metadata[OUT_DELTACC_BAND].varname, "OUT_DELTACC_BAND");
    strcpy(out_metadata[OUT_DELTACC_BAND].long_name, "deltacc_band");
    strcpy(out_metadata[OUT_DELTACC_BAND].standard_name,
           out_metadata[OUT_DELTACC].standard_name);
    strcpy(out_metadata[OUT_DELTACC_BAND].units,
           out_metadata[OUT_DELTACC].units);
    strcpy(out_metadata[OUT_DELTACC_BAND].description,
           out_metadata[OUT_DELTACC].description);

    /* net heat flux into ground [W m-2] */
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].varname, "OUT_GRND_FLUX_BAND");
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].long_name, "grnd_flux_band");
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].standard_name,
           out_metadata[OUT_GRND_FLUX].standard_name);
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].units,
           out_metadata[OUT_GRND_FLUX].units);
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].description,
           out_metadata[OUT_GRND_FLUX].description);

    /* incoming longwave flux at surface (under veg) [W m-2] */
    strcpy(out_metadata[OUT_IN_LONG_BAND].varname, "OUT_IN_LONG_BAND");
    strcpy(out_metadata[OUT_IN_LONG_BAND].long_name, "in_long_band");
    strcpy(out_metadata[OUT_IN_LONG_BAND].standard_name,
           out_metadata[OUT_IN_LONG].standard_name);
    strcpy(out_metadata[OUT_IN_LONG_BAND].units,
           out_metadata[OUT_IN_LONG].units);
    strcpy(out_metadata[OUT_IN_LONG_BAND].description,
           out_metadata[OUT_IN_LONG].description);

    /* net upward latent heat flux [W m-2] */
    strcpy(out_metadata[OUT_LATENT_BAND].varname, "OUT_LATENT_BAND");
    strcpy(out_metadata[OUT_LATENT_BAND].long_name, "latent_band");
    strcpy(out_metadata[OUT_LATENT_BAND].standard_name,
           out_metadata[OUT_LATENT].standard_name);
    strcpy(out_metadata[OUT_LATENT_BAND].units, out_metadata[OUT_LATENT].units);
    strcpy(out_metadata[OUT_LATENT_BAND].description,
           out_metadata[OUT_LATENT].description);

    /* net upward latent heat flux from sublimation [W m-2] */
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].varname, "OUT_LATENT_SUB_BAND");
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].long_name, "latent_sub_band");
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].standard_name,
           out_metadata[OUT_LATENT_SUB].standard_name);
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].units,
           out_metadata[OUT_LATENT_SUB].units);
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].description,
           out_metadata[OUT_LATENT_SUB].description);

    /* energy of fusion (melting) [W m-2] */
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].varname, "OUT_MELT_ENERGY_BAND");
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].long_name, "melt_energy_band");
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].standard_name,
           out_metadata[OUT_MELT_ENERGY].standard_name);
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].units,
           out_metadata[OUT_MELT_ENERGY].units);
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].description,
           out_metadata[OUT_MELT_ENERGY].description);

    /* net downward longwave flux [W m-2] */
    strcpy(out_metadata[OUT_LWNET_BAND].varname, "OUT_LWNET_BAND");
    strcpy(out_metadata[OUT_LWNET_BAND].long_name, "lwnet_band");
    strcpy(out_metadata[OUT_LWNET_BAND].standard_name,
           out_metadata[OUT_LWNET].standard_name);
    strcpy(out_metadata[OUT_LWNET_BAND].units, out_metadata[OUT_LWNET].units);
    strcpy(out_metadata[OUT_LWNET_BAND].description,
           out_metadata[OUT_LWNET].description);

    /* net downward shortwave flux [W m-2] */
    strcpy(out_metadata[OUT_SWNET_BAND].varname, "OUT_SWNET_BAND");
    strcpy(out_metadata[OUT_SWNET_BAND].long_name, "swnet_band");
    strcpy(out_metadata[OUT_SWNET_BAND].standard_name,
           out_metadata[OUT_SWNET].standard_name);
    strcpy(out_metadata[OUT_SWNET_BAND].units, out_metadata[OUT_SWNET].units);
    strcpy(out_metadata[OUT_SWNET_BAND].description,
           out_metadata[OUT_SWNET].description);

    /* net energy used to refreeze liquid water in snowpack [W m-2] */
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].varname, "OUT_RFRZ_ENERGY_BAND");
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].long_name, "rfrz_energy_band");
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].standard_name,
           out_metadata[OUT_RFRZ_ENERGY].standard_name);
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].units,
           out_metadata[OUT_RFRZ_ENERGY].units);
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].description,
           out_metadata[OUT_RFRZ_ENERGY].description);

    /* net upward sensible heat flux [W m-2] */
    strcpy(out_metadata[OUT_SENSIBLE_BAND].varname, "OUT_SENSIBLE_BAND");
    strcpy(out_metadata[OUT_SENSIBLE_BAND].long_name, "sensible_band");
    strcpy(out_metadata[OUT_SENSIBLE_BAND].standard_name,
           out_metadata[OUT_SENSIBLE].standard_name);
    strcpy(out_metadata[OUT_SENSIBLE_BAND].units,
           out_metadata[OUT_SENSIBLE].units);
    strcpy(out_metadata[OUT_SENSIBLE_BAND].description,
           out_metadata[OUT_SENSIBLE].description);

    /* snow interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].varname, "OUT_SNOW_CANOPY_BAND");
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].long_name, "snow_canopy_band");
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].standard_name,
           out_metadata[OUT_SNOW_CANOPY].standard_name);
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].units,
           out_metadata[OUT_SNOW_CANOPY].units);
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].description,
           out_metadata[OUT_SNOW_CANOPY].description);

    /* fractional area of snow cover [fraction] */
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].varname, "OUT_SNOW_COVER_BAND");
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].long_name, "snow_cover_band");
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].standard_name,
           out_metadata[OUT_SNOW_COVER].standard_name);
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].units,
           out_metadata[OUT_SNOW_COVER].units);
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].description,
           out_metadata[OUT_SNOW_COVER].description);

    /* depth of snow pack [cm] */
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].varname, "OUT_SNOW_DEPTH_BAND");
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].long_name, "snow_depth_band");
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].standard_name,
           out_metadata[OUT_SNOW_DEPTH].standard_name);
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].units,
           out_metadata[OUT_SNOW_DEPTH].units);
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].description,
           out_metadata[OUT_SNOW_DEPTH].description);

    /* energy flux through snow pack [W m-2] */
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].varname, "OUT_SNOW_FLUX_BAND");
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].long_name, "snow_flux_band");
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].standard_name,
           out_metadata[OUT_SNOW_FLUX].standard_name);
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].units,
           out_metadata[OUT_SNOW_FLUX].units);
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].description,
           out_metadata[OUT_SNOW_FLUX].description);

    /* snow melt [mm] */
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].varname, "OUT_SNOW_MELT_BAND");
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].long_name, "snow_melt_band");
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].standard_name,
           out_metadata[OUT_SNOW_MELT].standard_name);
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].units,
           out_metadata[OUT_SNOW_MELT].units);
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].description,
           out_metadata[OUT_SNOW_MELT].description);

    /* snow pack temperature [C] */
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].varname, "OUT_SNOW_PACKT_BAND");
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].long_name, "snow_packt_band");
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].standard_name,
           out_metadata[OUT_SNOW_PACK_TEMP].standard_name);
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].units,
           out_metadata[OUT_SNOW_PACK_TEMP].units);
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].description,
           out_metadata[OUT_SNOW_PACK_TEMP].description);

    /* snow surface temperature [C] */
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].varname, "OUT_SNOW_SURFT_BAND");
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].long_name, "snow_surft_band");
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].standard_name,
           out_metadata[OUT_SNOW_SURF_TEMP].standard_name);
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].units,
           out_metadata[OUT_SNOW_SURF_TEMP].units);
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].description,
           out_metadata[OUT_SNOW_SURF_TEMP].description);

    /* snow water equivalent in snow pack [mm] */
    strcpy(out_metadata[OUT_SWE_BAND].varname, "OUT_SWE_BAND");
    strcpy(out_metadata[OUT_SWE_BAND].long_name, "swe_band");
    strcpy(out_metadata[OUT_SWE_BAND].standard_name,
           out_metadata[OUT_SWE].standard_name);
    strcpy(out_metadata[OUT_SWE_BAND].units, out_metadata[OUT_SWE].units);
    strcpy(out_metadata[OUT_SWE_BAND].description,
           out_metadata[OUT_SWE].description);

    /* Wall time spent inside vic_run [seconds] */
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].varname, "OUT_TIME_VICRUN_WALL");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].long_name, "time_vicrun_wall");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].standard_name,
           "vic_run_wall_time");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].units, "seconds");
    strcpy(out_metadata[OUT_TIME_VICRUN_WALL].description,
           "Wall time spent inside vic_run");

    /* CPU time spent inside vic_run [seconds] */
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].varname, "OUT_TIME_VICRUN_CPU");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].long_name, "time_vicrun_cpu");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].standard_name, "vic_run_cpu_time");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].units, "seconds");
    strcpy(out_metadata[OUT_TIME_VICRUN_CPU].description,
           "CPU time spent inside vic_run");

    if (options.FROZEN_SOIL) {
        out_metadata[OUT_FDEPTH].nelem = MAX_FRONTS;
        out_metadata[OUT_TDEPTH].nelem = MAX_FRONTS;
    }

    out_metadata[OUT_SMLIQFRAC].nelem = options.Nlayer;
    out_metadata[OUT_SMFROZFRAC].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_ICE].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_LIQ].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_ICE_FRAC].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_LIQ_FRAC].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_MOIST].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_TEMP].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_TNODE].nelem = options.Nnode;
    out_metadata[OUT_SOIL_TNODE_WL].nelem = options.Nnode;
    out_metadata[OUT_SOILT_FBFLAG].nelem = options.Nnode;
    out_metadata[OUT_ADV_SENS_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_ADVECTION_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_ALBEDO_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_DELTACC_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_GRND_FLUX_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_IN_LONG_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_LATENT_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_LATENT_SUB_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_MELT_ENERGY_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_LWNET_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SWNET_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_RFRZ_ENERGY_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SENSIBLE_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_CANOPY_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_COVER_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_DEPTH_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_FLUX_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_MELT_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_PACKT_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SNOW_SURFT_BAND].nelem = options.SNOW_BAND;
    out_metadata[OUT_SWE_BAND].nelem = options.SNOW_BAND;
}
