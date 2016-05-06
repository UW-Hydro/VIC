/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine creates the list of output variables.
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

#include <vic_driver_shared_all.h>

void
set_output_met_data_info()
{
    size_t v;

    extern option_struct        options;
    extern out_metadata_struct  out_metadata[N_OUTVAR_TYPES];

    // Build the list of supported output variables

    // Set missing and/or default values
    for (v = 0; v < N_OUTVAR_TYPES; v++) {
        // Set default string values
       strcpy(out_metadata[v].varname, MISSING_S);
       strcpy(out_metadata[v].long_name, MISSING_S);
       strcpy(out_metadata[v].units, MISSING_S);
       strcpy(out_metadata[v].description, MISSING_S);
        // Set default number of elements
        out_metadata[v].nelem = 1;
    }

    // Water Balance Terms - state variables
    /* saturated area fraction */
    strcpy(out_metadata[OUT_ASAT].varname, "OUT_ASAT");
    strcpy(out_metadata[OUT_ASAT].long_name, "saturated_area_fraction");
    strcpy(out_metadata[OUT_ASAT].units, "1");
    strcpy(out_metadata[OUT_ASAT].description, "saturated area fraction");

    /* lake surface area as fraction of grid cell area [fraction] */
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].varname, "OUT_LAKE_AREA_FRAC");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].long_name, "lake_area_fraction");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].units, "1");
    strcpy(out_metadata[OUT_LAKE_AREA_FRAC].description, "lake surface area as fraction of grid cell area");

    /* lake depth [m] */
    strcpy(out_metadata[OUT_LAKE_DEPTH].varname, "OUT_LAKE_DEPTH");
    strcpy(out_metadata[OUT_LAKE_DEPTH].long_name, "lake_depth");
    strcpy(out_metadata[OUT_LAKE_DEPTH].units, "m");
    strcpy(out_metadata[OUT_LAKE_DEPTH].description, "lake depth");

    /* moisture stored as lake ice [mm] */
    strcpy(out_metadata[OUT_LAKE_ICE].varname, "OUT_LAKE_ICE");
    strcpy(out_metadata[OUT_LAKE_ICE].long_name, "lake_ice");
    strcpy(out_metadata[OUT_LAKE_ICE].units, "mm");
    strcpy(out_metadata[OUT_LAKE_ICE].description, "moisture stored as lake ice");

    /* fractional coverage of lake ice [fraction] */
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].varname, "OUT_LAKE_ICE_FRACT");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].long_name, "lake_ice_area_fraction");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].units, "1");
    strcpy(out_metadata[OUT_LAKE_ICE_FRACT].description, "fractional coverage of lake ice");

    /* thickness of lake ice [cm] */
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].varname, "OUT_LAKE_ICE_HEIGHT");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].long_name, "lake_ice_height");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].units, "cm");
    strcpy(out_metadata[OUT_LAKE_ICE_HEIGHT].description, "something");

    /* liquid water stored in lake [mm over lake area?] */
    strcpy(out_metadata[OUT_LAKE_MOIST].varname, "OUT_LAKE_MOIST");
    strcpy(out_metadata[OUT_LAKE_MOIST].long_name, "lake_moisture");
    strcpy(out_metadata[OUT_LAKE_MOIST].units, "mm");
    strcpy(out_metadata[OUT_LAKE_MOIST].description, "thickness of lake ice");

    /* lake surface area [m2] */
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].varname, "OUT_LAKE_SURF_AREA");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].long_name, "lake_area");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].units, "m2");
    strcpy(out_metadata[OUT_LAKE_SURF_AREA].description, "lake surface area");

    /* liquid water equivalent of snow on top of lake ice [m over lake ice] */
    strcpy(out_metadata[OUT_LAKE_SWE].varname, "OUT_LAKE_SWE");
    strcpy(out_metadata[OUT_LAKE_SWE].long_name, "lwe_thickness_of_snow_on_lake");
    strcpy(out_metadata[OUT_LAKE_SWE].units, "m");
    strcpy(out_metadata[OUT_LAKE_SWE].description, "something");

    /* volumetric liquid water equivalent of snow on top of lake ice [m3] */
    strcpy(out_metadata[OUT_LAKE_SWE_V].varname, "OUT_LAKE_SWE_V");
    strcpy(out_metadata[OUT_LAKE_SWE_V].long_name, "lwe_volume_of_snow_on_lake");
    strcpy(out_metadata[OUT_LAKE_SWE_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_SWE_V].description, "liquid water equivalent of snow on top of lake ice");

    /* lake volume [m3] */
    strcpy(out_metadata[OUT_LAKE_VOLUME].varname, "OUT_LAKE_VOLUME");
    strcpy(out_metadata[OUT_LAKE_VOLUME].long_name, "lake_volume");
    strcpy(out_metadata[OUT_LAKE_VOLUME].units, "m3");
    strcpy(out_metadata[OUT_LAKE_VOLUME].description, "lake volume");

    /* root zone soil moisture [mm] */
    strcpy(out_metadata[OUT_ROOTMOIST].varname, "OUT_ROOTMOIST");
    strcpy(out_metadata[OUT_ROOTMOIST].long_name, "soil_moisture_in_root_zone");
    strcpy(out_metadata[OUT_ROOTMOIST].units, "mm");
    strcpy(out_metadata[OUT_ROOTMOIST].description, "root zone soil moisture");

    /* fraction of soil moisture (by mass) that is ice, for each soil layer */
    strcpy(out_metadata[OUT_SMFROZFRAC].varname, "OUT_SMFROZFRAC");
    strcpy(out_metadata[OUT_SMFROZFRAC].long_name, "soil_moisture_ice_fraction");
    strcpy(out_metadata[OUT_SMFROZFRAC].units, "1");
    strcpy(out_metadata[OUT_SMFROZFRAC].description, "fraction of soil moisture (by mass) that is ice, for each soil layer");

    /* fraction of soil moisture (by mass) that is liquid, for each soil layer */
    strcpy(out_metadata[OUT_SMLIQFRAC].varname, "OUT_SMLIQFRAC");
    strcpy(out_metadata[OUT_SMLIQFRAC].long_name, "soil_moisture_liquid_fraction");
    strcpy(out_metadata[OUT_SMLIQFRAC].units, "1");
    strcpy(out_metadata[OUT_SMLIQFRAC].description, "fraction of soil moisture (by mass) that is liquid, for each soil layer");

    /* snow interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_SNOW_CANOPY].varname, "OUT_SNOW_CANOPY");
    strcpy(out_metadata[OUT_SNOW_CANOPY].long_name, "lwe_thickness_of_canopy_intercepted_snow");
    strcpy(out_metadata[OUT_SNOW_CANOPY].units, "mm");
    strcpy(out_metadata[OUT_SNOW_CANOPY].description, "snow interception storage in canopy");

    /* fractional area of snow cover [fraction] */
    strcpy(out_metadata[OUT_SNOW_COVER].varname, "OUT_SNOW_COVER");
    strcpy(out_metadata[OUT_SNOW_COVER].long_name, "snow_cover_area_fraction");
    strcpy(out_metadata[OUT_SNOW_COVER].units, "1");
    strcpy(out_metadata[OUT_SNOW_COVER].description, "fractional area of snow cover");

    /* depth of snow pack [cm] */
    strcpy(out_metadata[OUT_SNOW_DEPTH].varname, "OUT_SNOW_DEPTH");
    strcpy(out_metadata[OUT_SNOW_DEPTH].long_name, "thickness_of_snow");
    strcpy(out_metadata[OUT_SNOW_DEPTH].units, "cm");
    strcpy(out_metadata[OUT_SNOW_DEPTH].description, "depth of snow pack");

    /* soil ice content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_ICE].varname, "OUT_SOIL_ICE");
    strcpy(out_metadata[OUT_SOIL_ICE].long_name, "soil_moisture_ice_depth");
    strcpy(out_metadata[OUT_SOIL_ICE].units, "mm");
    strcpy(out_metadata[OUT_SOIL_ICE].description, "soil ice content for each soil layer");

    /* soil liquid moisture content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_LIQ].varname, "OUT_SOIL_LIQ");
    strcpy(out_metadata[OUT_SOIL_LIQ].long_name, "soil_moisture_liquid_depth");
    strcpy(out_metadata[OUT_SOIL_LIQ].units, "mm");
    strcpy(out_metadata[OUT_SOIL_LIQ].description, "soil liquid moisture content for each soil layer");

    /* soil ice content [1] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].varname, "OUT_SOIL_ICE_FRAC");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].long_name, "soil_moisture_ice_depth_fraction");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].units, "1");
    strcpy(out_metadata[OUT_SOIL_ICE_FRAC].description, "soil ice content fraction for each soil layer");

    /* soil liquid moisture content [1] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].varname, "OUT_SOIL_LIQ_FRAC");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].long_name, "soil_moisture_liquid_depth_fraction");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].units, "1");
    strcpy(out_metadata[OUT_SOIL_LIQ_FRAC].description, "soil liquid moisture content fraction for each soil layer");

    /* soil total moisture content [mm] for each soil layer */
    strcpy(out_metadata[OUT_SOIL_MOIST].varname, "OUT_SOIL_MOIST");
    strcpy(out_metadata[OUT_SOIL_MOIST].long_name, "unset");
    strcpy(out_metadata[OUT_SOIL_MOIST].units, "unset");
    strcpy(out_metadata[OUT_SOIL_MOIST].description, "something");

    /* vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
    strcpy(out_metadata[OUT_SOIL_WET].varname, "OUT_SOIL_WET");
    strcpy(out_metadata[OUT_SOIL_WET].long_name, "soil_moisture_wetness_fraction");
    strcpy(out_metadata[OUT_SOIL_WET].units, "1");
    strcpy(out_metadata[OUT_SOIL_WET].description, "vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point)");

    /* storage of liquid water on surface (ponding) [mm] */
    strcpy(out_metadata[OUT_SURFSTOR].varname, "OUT_SURFSTOR");
    strcpy(out_metadata[OUT_SURFSTOR].long_name, "surface_liquid_water_storage");
    strcpy(out_metadata[OUT_SURFSTOR].units, "mm");
    strcpy(out_metadata[OUT_SURFSTOR].description, "storage of liquid water on surface (ponding)");

    /* fraction of soil surface that is frozen [fraction] */
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].varname, "OUT_SURF_FROST_FRAC");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].long_name, "frozen_soil_surface_fraction");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].units, "1");
    strcpy(out_metadata[OUT_SURF_FROST_FRAC].description, "fraction of soil surface that is frozen");

    /* snow water equivalent in snow pack [mm] */
    strcpy(out_metadata[OUT_SWE].varname, "OUT_SWE");
    strcpy(out_metadata[OUT_SWE].long_name, "lwe_thickness_of_snow");
    strcpy(out_metadata[OUT_SWE].units, "mm");
    strcpy(out_metadata[OUT_SWE].description, "snow water equivalent in snow pack");

    /* total moisture interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_WDEW].varname, "OUT_WDEW");
    strcpy(out_metadata[OUT_WDEW].long_name, "soil_moisture_storage_depth_in_canopy");
    strcpy(out_metadata[OUT_WDEW].units, "mm");
    strcpy(out_metadata[OUT_WDEW].description, "total moisture interception storage in canopy");

    /* water table position [cm] (zwt within lowest unsaturated layer) */
    strcpy(out_metadata[OUT_ZWT].varname, "OUT_ZWT");
    strcpy(out_metadata[OUT_ZWT].long_name, "water_table_position_lowest_layer");
    strcpy(out_metadata[OUT_ZWT].units, "cm");
    strcpy(out_metadata[OUT_ZWT].description, "water table position (zwt within lowest unsaturated layer)");

    /* lumped water table position [cm] (zwt of total moisture across all layers, lumped together) */
    strcpy(out_metadata[OUT_ZWT_LUMPED].varname, "OUT_ZWT_LUMPED");
    strcpy(out_metadata[OUT_ZWT_LUMPED].long_name, "lumped_water_table_position");
    strcpy(out_metadata[OUT_ZWT_LUMPED].units, "cm");
    strcpy(out_metadata[OUT_ZWT_LUMPED].description, "lumped water table position (zwt of total moisture across all layers, lumped together)");

    // Water Balance Terms - fluxes
    /* baseflow out of the bottom layer [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_BASEFLOW].varname, "OUT_BASEFLOW");
    strcpy(out_metadata[OUT_BASEFLOW].long_name, "baseflow_amount");
    strcpy(out_metadata[OUT_BASEFLOW].units, "mm");
    strcpy(out_metadata[OUT_BASEFLOW].description, "baseflow out of the bottom layer");

    /* change in canopy interception storage [mm] */
    strcpy(out_metadata[OUT_DELINTERCEPT].varname, "OUT_DELINTERCEPT");
    strcpy(out_metadata[OUT_DELINTERCEPT].long_name, "change_in_canopy_interception_amount");
    strcpy(out_metadata[OUT_DELINTERCEPT].units, "mm");
    strcpy(out_metadata[OUT_DELINTERCEPT].description, "change in canopy interception storage");

    /* change in soil water content [mm] */
    strcpy(out_metadata[OUT_DELSOILMOIST].varname, "OUT_DELSOILMOIST");
    strcpy(out_metadata[OUT_DELSOILMOIST].long_name, "change_in_soil_water_amount");
    strcpy(out_metadata[OUT_DELSOILMOIST].units, "mm");
    strcpy(out_metadata[OUT_DELSOILMOIST].description, "change in soil water content");

    /* change in snow water equivalent [mm] */
    strcpy(out_metadata[OUT_DELSWE].varname, "OUT_DELSWE");
    strcpy(out_metadata[OUT_DELSWE].long_name, "change_in_snow_lwe_thickness");
    strcpy(out_metadata[OUT_DELSWE].units, "mm");
    strcpy(out_metadata[OUT_DELSWE].description, "change in snow water equivalent");

    /* change in surface liquid water storage  [mm] */
    strcpy(out_metadata[OUT_DELSURFSTOR].varname, "OUT_DELSURFSTOR");
    strcpy(out_metadata[OUT_DELSURFSTOR].long_name, "change_in_surface_liquid_water_storage");
    strcpy(out_metadata[OUT_DELSURFSTOR].units, "mm");
    strcpy(out_metadata[OUT_DELSURFSTOR].description, "change in surface liquid water storage");

    /* total net evaporation [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_EVAP].varname, "OUT_EVAP");
    strcpy(out_metadata[OUT_EVAP].long_name, "water_evaporation_flux_net");
    strcpy(out_metadata[OUT_EVAP].units, "mm");
    strcpy(out_metadata[OUT_EVAP].description, "total net evaporation");

    /* net evaporation from bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_EVAP_BARE].varname, "OUT_EVAP_BARE");
    strcpy(out_metadata[OUT_EVAP_BARE].long_name, "water_evaporation_from_soil");
    strcpy(out_metadata[OUT_EVAP_BARE].units, "mm");
    strcpy(out_metadata[OUT_EVAP_BARE].description, "net evaporation from bare soil");

    /* net evaporation from canopy interception [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_EVAP_CANOP].varname, "OUT_EVAP_CANOP");
    strcpy(out_metadata[OUT_EVAP_CANOP].long_name, "water_evaporation_from_canopy");
    strcpy(out_metadata[OUT_EVAP_CANOP].units, "mm");
    strcpy(out_metadata[OUT_EVAP_CANOP].description, "net evaporation from canopy interception");

    /* moisture that reaches top of soil column [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_INFLOW].varname, "OUT_INFLOW");
    strcpy(out_metadata[OUT_INFLOW].long_name, "unset");
    strcpy(out_metadata[OUT_INFLOW].units, "unset");
    strcpy(out_metadata[OUT_INFLOW].description, "something");

    /* incoming baseflow from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_BF_IN].varname, "OUT_LAKE_BF_IN");
    strcpy(out_metadata[OUT_LAKE_BF_IN].long_name, "infiltration_amount");
    strcpy(out_metadata[OUT_LAKE_BF_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_BF_IN].description, "oisture that reaches top of soil column");

    /* incoming volumetric baseflow from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].varname, "OUT_LAKE_BF_IN_V");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].long_name, "baseflow_volume_from_lake_catchment");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_BF_IN_V].description, "incoming volumetric baseflow from lake catchment");

    /* outgoing baseflow lake [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_BF_OUT].varname, "OUT_LAKE_BF_OUT");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].long_name, "baseflow_outgoing_lake");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].units, "mm");
    strcpy(out_metadata[OUT_LAKE_BF_OUT].description, "somethoutgoing baseflow lakeing");

    /* outgoing volumetric baseflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].varname, "OUT_LAKE_BF_OUT_V");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].long_name, "baseflow_outgoing_volume_lake");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_BF_OUT_V].description, "outgoing volumetric baseflow from lake");

    /* channel inflow into lake [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].varname, "OUT_LAKE_CHAN_IN");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].long_name, "channel_inflow_into_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN].description, "channel inflow into lake");

    /* volumetric channel inflow into lake [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].varname, "OUT_LAKE_CHAN_IN_V");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].long_name, "channel_inflow_volume_into_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_CHAN_IN_V].description, "volumetric channel inflow into lake");

    /* channel outflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].varname, "OUT_LAKE_CHAN_OUT");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].long_name, "channel_outflow_from_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].units, "mm");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT].description, "channel outflow from lake");

    /* volumetric channel outflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].varname, "OUT_LAKE_CHAN_OUT_V");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].long_name, "channel_outflow_volume_from_lake");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_CHAN_OUT_V].description, "volumetric channel outflow from lake");

    /* change in lake moisture storage (liquid plus ice cover) [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_DSTOR].varname, "OUT_LAKE_DSTOR");
    strcpy(out_metadata[OUT_LAKE_DSTOR].long_name, "change_in_lake_moisture_amount");
    strcpy(out_metadata[OUT_LAKE_DSTOR].units, "mm");
    strcpy(out_metadata[OUT_LAKE_DSTOR].description, "change in lake moisture storage (liquid plus ice cover)");

    /* volumetric change in lake moisture storage (liquid plus ice cover) [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].varname, "OUT_LAKE_DSTOR_V");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].long_name, "change_in_lake_moisture_volume");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_DSTOR_V].description, "volumetric change in lake moisture storage (liquid plus ice cover)");

    /* change in snowpack on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_DSWE].varname, "OUT_LAKE_DSWE");
    strcpy(out_metadata[OUT_LAKE_DSWE].long_name, "change_in_snow_lwe_thickness_on_lake_ice");
    strcpy(out_metadata[OUT_LAKE_DSWE].units, "mm");
    strcpy(out_metadata[OUT_LAKE_DSWE].description, "change in snowpack on top of lake ice");

    /* volumetric change in snowpack on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_DSWE_V].varname, "OUT_LAKE_DSWE_V");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].long_name, "change_in_snow_lwe_volume_on_lake_ice");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_DSWE_V].description, "volumetric change in snowpack on top of lake ice");

    /* net evaporation from lake surface [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_EVAP].varname, "OUT_LAKE_EVAP");
    strcpy(out_metadata[OUT_LAKE_EVAP].long_name, "water_evaporation_from_lake");
    strcpy(out_metadata[OUT_LAKE_EVAP].units, "mm");
    strcpy(out_metadata[OUT_LAKE_EVAP].description, "net evaporation from lake surface");

    /* net volumetric evaporation from lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_EVAP_V].varname, "OUT_LAKE_EVAP_V");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].long_name, "water_evaporation_volume_from_lake");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_EVAP_V].description, "net volumetric evaporation from lake surface");

    /* volumetric precipitation over lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_PREC_V].varname, "OUT_LAKE_PREC_V");
    strcpy(out_metadata[OUT_LAKE_PREC_V].long_name, "precipitation_over_lake_volume");
    strcpy(out_metadata[OUT_LAKE_PREC_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_PREC_V].description, "volumetric precipitation over lake surface");

    /* recharge from lake to surrounding wetland [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_RCHRG].varname, "OUT_LAKE_RCHRG");
    strcpy(out_metadata[OUT_LAKE_RCHRG].long_name, "recharge_from_lake_to_wetland");
    strcpy(out_metadata[OUT_LAKE_RCHRG].units, "mm");
    strcpy(out_metadata[OUT_LAKE_RCHRG].description, "recharge from lake to surrounding wetland");

    /* volumetric recharge from lake to surrounding wetland [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].varname, "OUT_LAKE_RCHRG_V");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].long_name, "");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].units, "unset");
    strcpy(out_metadata[OUT_LAKE_RCHRG_V].description, "volumetric recharge from lake to surrounding wetland");

    /* incoming runoff from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_RO_IN].varname, "OUT_LAKE_RO_IN");
    strcpy(out_metadata[OUT_LAKE_RO_IN].long_name, "recharge_volume_from_lake_to_wetland");
    strcpy(out_metadata[OUT_LAKE_RO_IN].units, "mm");
    strcpy(out_metadata[OUT_LAKE_RO_IN].description, "something");

    /* incoming volumetric runoff from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].varname, "OUT_LAKE_RO_IN_V");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].long_name, "runoff_volume_from_lake_catchment");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_RO_IN_V].description, "incoming volumetric runoff from lake catchment");

    /* sublimation from lake snow pack [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_LAKE_VAPFLX].varname, "OUT_LAKE_VAPFLX");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].long_name, "water_sublimation_flux_from_lake_snow");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].units, "mm");
    strcpy(out_metadata[OUT_LAKE_VAPFLX].description, "sublimation from lake snow pack");

    /* volumetric sublimation from lake snow pack [m3] (ALMA_OUTPUT: [m3/s]) */
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].varname, "OUT_LAKE_VAPFLX_V");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].long_name, "water_sublimation_flux_volume_from_lake_snow");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].units, "m3");
    strcpy(out_metadata[OUT_LAKE_VAPFLX_V].description, "volumetric sublimation from lake snow pack");

    /* Potential evapotranspiration (= area-weighted sum of potential transpiration and potential soil evaporation). [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_PET].varname, "OUT_PET");
    strcpy(out_metadata[OUT_PET].long_name, "water_potential_evaporation_amount");
    strcpy(out_metadata[OUT_PET].units, "mm");
    strcpy(out_metadata[OUT_PET].description, "Potential evapotranspiration (= area-weighted sum of potential transpiration and potential soil evaporation)");

    /* incoming precipitation [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_PREC].varname, "OUT_PREC");
    strcpy(out_metadata[OUT_PREC].long_name, "precipitation_amount");
    strcpy(out_metadata[OUT_PREC].units, "mm");
    strcpy(out_metadata[OUT_PREC].description, "incoming precipitation");

    /* rainfall [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_RAINF].varname, "OUT_RAINF");
    strcpy(out_metadata[OUT_RAINF].long_name, "rainfall_amount");
    strcpy(out_metadata[OUT_RAINF].units, "mm");
    strcpy(out_metadata[OUT_RAINF].description, "liquid rainfall amount");

    /* refreezing of water in the snow [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_REFREEZE].varname, "OUT_REFREEZE");
    strcpy(out_metadata[OUT_REFREEZE].long_name, "lwe_thickness_of_refreezing_water_in_snow");
    strcpy(out_metadata[OUT_REFREEZE].units, "mm");
    strcpy(out_metadata[OUT_REFREEZE].description, "refreezing of water in the snow");

    /* surface runoff [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_RUNOFF].varname, "OUT_RUNOFF");
    strcpy(out_metadata[OUT_RUNOFF].long_name, "runoff_amount");
    strcpy(out_metadata[OUT_RUNOFF].units, "mm");
    strcpy(out_metadata[OUT_RUNOFF].description, "surface runoff");

    /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SNOW_MELT].varname, "OUT_SNOW_MELT");
    strcpy(out_metadata[OUT_SNOW_MELT].long_name, "snow_melt_amount");
    strcpy(out_metadata[OUT_SNOW_MELT].units, "mm");
    strcpy(out_metadata[OUT_SNOW_MELT].description, "snow melt");

    /* snowfall [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SNOWF].varname, "OUT_SNOWF");
    strcpy(out_metadata[OUT_SNOWF].long_name, "snowfall_lwe_amount");
    strcpy(out_metadata[OUT_SNOWF].units, "mm");
    strcpy(out_metadata[OUT_SNOWF].description, "snowfall");

    /* net sublimation of blowing snow [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SUB_BLOWING].varname, "OUT_SUB_BLOWING");
    strcpy(out_metadata[OUT_SUB_BLOWING].long_name, "submlimation_amount_from_blowing_snow");
    strcpy(out_metadata[OUT_SUB_BLOWING].units, "mm");
    strcpy(out_metadata[OUT_SUB_BLOWING].description, "net sublimation of blowing snow");

    /* net sublimation from snow stored in canopy [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SUB_CANOP].varname, "OUT_SUB_CANOP");
    strcpy(out_metadata[OUT_SUB_CANOP].long_name, "sublimation_amount_from_canopy_snow");
    strcpy(out_metadata[OUT_SUB_CANOP].units, "mm");
    strcpy(out_metadata[OUT_SUB_CANOP].description, "net sublimation from snow stored in canopy");

    /* net sublimation from snow pack (surface and blowing) [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SUB_SNOW].varname, "OUT_SUB_SNOW");
    strcpy(out_metadata[OUT_SUB_SNOW].long_name, "sublimation_amount_from_snow_pack");
    strcpy(out_metadata[OUT_SUB_SNOW].units, "mm");
    strcpy(out_metadata[OUT_SUB_SNOW].description, "net sublimation from snow pack (surface and blowing)");

    /* net sublimation from snow pack surface [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SUB_SURFACE].varname, "OUT_SUB_SURFACE");
    strcpy(out_metadata[OUT_SUB_SURFACE].long_name, "sublimation_amount_from_snow_pack_surface");
    strcpy(out_metadata[OUT_SUB_SURFACE].units, "mm");
    strcpy(out_metadata[OUT_SUB_SURFACE].description, "net sublimation from snow pack surface");

    /* net transpiration from vegetation [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_TRANSP_VEG].varname, "OUT_TRANSP_VEG");
    strcpy(out_metadata[OUT_TRANSP_VEG].long_name, "transpiration_amount");
    strcpy(out_metadata[OUT_TRANSP_VEG].units, "mm");
    strcpy(out_metadata[OUT_TRANSP_VEG].description, "net transpiration from vegetation");

    // Energy Balance Terms - state variables
    /* albedo [fraction] */
    strcpy(out_metadata[OUT_ALBEDO].varname, "OUT_ALBEDO");
    strcpy(out_metadata[OUT_ALBEDO].long_name, "surface_albedo");
    strcpy(out_metadata[OUT_ALBEDO].units, "1");
    strcpy(out_metadata[OUT_ALBEDO].description, "albedo");

    /* bare soil surface temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_BARESOILT].varname, "OUT_BARESOILT");
    strcpy(out_metadata[OUT_BARESOILT].long_name, "surface_temperature_bare_soil");
    strcpy(out_metadata[OUT_BARESOILT].units, "C");
    strcpy(out_metadata[OUT_BARESOILT].description, "bare soil surface temperature");

    /* depth of freezing fronts [cm] (ALMA_OUTPUT: [m]) for each freezing front */
    strcpy(out_metadata[OUT_FDEPTH].varname, "OUT_FDEPTH");
    strcpy(out_metadata[OUT_FDEPTH].long_name, "freezing_fronts_depth");
    strcpy(out_metadata[OUT_FDEPTH].units, "cm");
    strcpy(out_metadata[OUT_FDEPTH].description, "depth of freezing fronts for each freezing front");

    /* lake ice temperature [K] */
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].varname, "OUT_LAKE_ICE_TEMP");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].long_name, "lake_ice_temperature");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].units, "K");
    strcpy(out_metadata[OUT_LAKE_ICE_TEMP].description, "lake ice temperature");

    /* lake surface temperature [K] */
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].varname, "OUT_LAKE_SURF_TEMP");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].long_name, "lake_surface_temperature");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].units, "K");
    strcpy(out_metadata[OUT_LAKE_SURF_TEMP].description, "lake surface temperature");

    /* average radiative surface temperature [K] */
    strcpy(out_metadata[OUT_RAD_TEMP].varname, "OUT_RAD_TEMP");
    strcpy(out_metadata[OUT_RAD_TEMP].long_name, "surface_radiative_temperature");
    strcpy(out_metadata[OUT_RAD_TEMP].units, "K");
    strcpy(out_metadata[OUT_RAD_TEMP].description, "average radiative surface temperature");

    /* snow albedo [fraction] */
    strcpy(out_metadata[OUT_SALBEDO].varname, "OUT_SALBEDO");
    strcpy(out_metadata[OUT_SALBEDO].long_name, "snow_albedo");
    strcpy(out_metadata[OUT_SALBEDO].units, "1");
    strcpy(out_metadata[OUT_SALBEDO].description, "snow albedo");

    /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].varname, "OUT_SNOW_PACK_TEMP");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].long_name, "snow_pack_temperature");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].units, "C");
    strcpy(out_metadata[OUT_SNOW_PACK_TEMP].description, "snow pack temperature");

    /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].varname, "OUT_SNOW_SURF_TEMP");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].long_name, "snow_surface_temperature");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].units, "C");
    strcpy(out_metadata[OUT_SNOW_SURF_TEMP].description, "snow surface temperature");

    /* snow surface temperature flag */
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].varname, "OUT_SNOWT_FBFLAG");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].long_name, "snow_surface_temperature_flag");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SNOWT_FBFLAG].description, "snow surface temperature flag");

    /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer */
    strcpy(out_metadata[OUT_SOIL_TEMP].varname, "OUT_SOIL_TEMP");
    strcpy(out_metadata[OUT_SOIL_TEMP].long_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TEMP].units, "C");
    strcpy(out_metadata[OUT_SOIL_TEMP].description, "soil temperature for each soil layer");

    /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node */
    strcpy(out_metadata[OUT_SOIL_TNODE].varname, "OUT_SOIL_TNODE");
    strcpy(out_metadata[OUT_SOIL_TNODE].long_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TNODE].units, "C");
    strcpy(out_metadata[OUT_SOIL_TNODE].description, "soil temperature for each soil thermal node");

    /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node in the wetland */
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].varname, "OUT_SOIL_TNODE_WL");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].long_name, "soil_temperature");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].units, "C");
    strcpy(out_metadata[OUT_SOIL_TNODE_WL].description, "soil temperature for each soil thermal node in the wetland");

    /* soil temperature flag for each soil thermal node */
    strcpy(out_metadata[OUT_SOILT_FBFLAG].varname, "OUT_SOILT_FBFLAG");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].long_name, "soil_temperature_flag");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SOILT_FBFLAG].description, "soil temperature flag for each soil thermal node");

    /* average surface temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_SURF_TEMP].varname, "OUT_SURF_TEMP");
    strcpy(out_metadata[OUT_SURF_TEMP].long_name, "surface_temperature");
    strcpy(out_metadata[OUT_SURF_TEMP].units, "C");
    strcpy(out_metadata[OUT_SURF_TEMP].description, "average surface temperature");

    /* surface temperature flag */
    strcpy(out_metadata[OUT_SURFT_FBFLAG].varname, "OUT_SURFT_FBFLAG");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].long_name, "surface_temperature_flag");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_SURFT_FBFLAG].description, "surface temperature flag");

    /* Tcanopy flag */
    strcpy(out_metadata[OUT_TCAN_FBFLAG].varname, "OUT_TCAN_FBFLAG");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].long_name, "canopy_temperature_flag");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TCAN_FBFLAG].description, "Tcanopy flag");

    /* depth of thawing fronts [cm] (ALMA_OUTPUT: [m]) for each thawing front */
    strcpy(out_metadata[OUT_TDEPTH].varname, "OUT_TDEPTH");
    strcpy(out_metadata[OUT_TDEPTH].long_name, "depth_of_thawing_fronts");
    strcpy(out_metadata[OUT_TDEPTH].units, "cm");
    strcpy(out_metadata[OUT_TDEPTH].description, "depth of thawing fronts for each thawing front");

    /* Tfoliage flag */
    strcpy(out_metadata[OUT_TFOL_FBFLAG].varname, "OUT_TFOL_FBFLAG");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].long_name, "foliage_temperature_flag");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].units, "1");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].description, "depth of thawing fronts for each thawing front");

    /* average vegetation canopy temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_VEGT].varname, "OUT_VEGT");
    strcpy(out_metadata[OUT_VEGT].long_name, "canopy_temperature");
    strcpy(out_metadata[OUT_VEGT].units, "C");
    strcpy(out_metadata[OUT_TFOL_FBFLAG].description, "average vegetation canopy temperature");

    // Energy Balance Terms - fluxes
    /* net sensible heat advected to snow pack [W/m2] */
    strcpy(out_metadata[OUT_ADV_SENS].varname, "OUT_ADV_SENS");
    strcpy(out_metadata[OUT_ADV_SENS].long_name, "net_sensible_heat_flux_to_snow_pack");
    strcpy(out_metadata[OUT_ADV_SENS].units, "W m-2");
    strcpy(out_metadata[OUT_ADV_SENS].description, "net sensible heat advected to snow pack");

    /* advected energy [W/m2] */
    strcpy(out_metadata[OUT_ADVECTION].varname, "OUT_ADVECTION");
    strcpy(out_metadata[OUT_ADVECTION].long_name, "advected_energy");
    strcpy(out_metadata[OUT_ADVECTION].units, "W m-2");
    strcpy(out_metadata[OUT_ADVECTION].description, "advected energy ");

    /* rate of change in cold content in snow pack [W/m2] */
    strcpy(out_metadata[OUT_DELTACC].varname, "OUT_DELTACC");
    strcpy(out_metadata[OUT_DELTACC].long_name, "rate_change_in_snow_pack_cold_content");
    strcpy(out_metadata[OUT_DELTACC].units, "W m-2");
    strcpy(out_metadata[OUT_DELTACC].description, "rate of change in cold content in snow pack");

    /* rate of change in heat storage [W/m2] */
    strcpy(out_metadata[OUT_DELTAH].varname, "OUT_DELTAH");
    strcpy(out_metadata[OUT_DELTAH].long_name, "rate_change_in_heat_storage");
    strcpy(out_metadata[OUT_DELTAH].units, "W m-2");
    strcpy(out_metadata[OUT_DELTAH].description, "rate of change in heat storage");

    /* energy budget error [W/m2] */
    strcpy(out_metadata[OUT_ENERGY_ERROR].varname, "OUT_ENERGY_ERROR");
    strcpy(out_metadata[OUT_ENERGY_ERROR].long_name, "energy_budget_error");
    strcpy(out_metadata[OUT_ENERGY_ERROR].units, "W m-2");
    strcpy(out_metadata[OUT_ENERGY_ERROR].description, "energy budget error");

    /* water budget error [mm] */
    strcpy(out_metadata[OUT_WATER_ERROR].varname, "OUT_WATER_ERROR");
    strcpy(out_metadata[OUT_WATER_ERROR].long_name, "water_budget_error");
    strcpy(out_metadata[OUT_WATER_ERROR].units, "mm");
    strcpy(out_metadata[OUT_WATER_ERROR].description, "water budget error");

    /* net energy used to melt/freeze soil moisture [W/m2] */
    strcpy(out_metadata[OUT_FUSION].varname, "OUT_FUSION");
    strcpy(out_metadata[OUT_FUSION].long_name, "energy_of_fusion_in_soil_moisture");
    strcpy(out_metadata[OUT_FUSION].units, "W m-2");
    strcpy(out_metadata[OUT_FUSION].description, "net energy used to melt/freeze soil moisture");

    /* net heat flux into ground [W/m2] */
    strcpy(out_metadata[OUT_GRND_FLUX].varname, "OUT_GRND_FLUX");
    strcpy(out_metadata[OUT_GRND_FLUX].long_name, "downward_heat_flux_in_soil");
    strcpy(out_metadata[OUT_GRND_FLUX].units, "W m-2");
    strcpy(out_metadata[OUT_GRND_FLUX].description, "net heat flux into ground");

    /* incoming longwave flux at surface (under veg) [W/m2] */
    strcpy(out_metadata[OUT_IN_LONG].varname, "OUT_IN_LONG");
    strcpy(out_metadata[OUT_IN_LONG].long_name, "downwelling_longwave_flux_at_ground_surface");
    strcpy(out_metadata[OUT_IN_LONG].units, "W m-2");
    strcpy(out_metadata[OUT_IN_LONG].description, "incoming longwave flux at surface (under veg)");

    /* net upward latent heat flux [W/m2] */
    strcpy(out_metadata[OUT_LATENT].varname, "OUT_LATENT");
    strcpy(out_metadata[OUT_LATENT].long_name, "surface_upward_latent_heat_flux");
    strcpy(out_metadata[OUT_LATENT].units, "W m-2");
    strcpy(out_metadata[OUT_LATENT].description, "net upward latent heat flux");

    /* net upward latent heat flux from sublimation [W/m2] */
    strcpy(out_metadata[OUT_LATENT_SUB].varname, "OUT_LATENT_SUB");
    strcpy(out_metadata[OUT_LATENT_SUB].long_name, "surface_net_latent_heat_flux_from_sublimation");
    strcpy(out_metadata[OUT_LATENT_SUB].units, "W m-2");
    strcpy(out_metadata[OUT_LATENT_SUB].description, "net upward latent heat flux from sublimation");

    /* energy of fusion (melting) [W/m2] */
    strcpy(out_metadata[OUT_MELT_ENERGY].varname, "OUT_MELT_ENERGY");
    strcpy(out_metadata[OUT_MELT_ENERGY].long_name, "energy_of_fusion");
    strcpy(out_metadata[OUT_MELT_ENERGY].units, "W m-2");
    strcpy(out_metadata[OUT_MELT_ENERGY].description, "energy of fusion (melting)");

    /* net downward longwave flux [W/m2] */
    strcpy(out_metadata[OUT_LWNET].varname, "OUT_LWNET");
    strcpy(out_metadata[OUT_LWNET].long_name, "net_downward_longwave_flux_at_surface");
    strcpy(out_metadata[OUT_LWNET].units, "W m-2");
    strcpy(out_metadata[OUT_LWNET].description, "net downward longwave flux");

    /* net downward shortwave flux [W/m2] */
    strcpy(out_metadata[OUT_SWNET].varname, "OUT_SWNET");
    strcpy(out_metadata[OUT_SWNET].long_name, "net_downward_shortwave_flux_at_surface");
    strcpy(out_metadata[OUT_SWNET].units, "W m-2");
    strcpy(out_metadata[OUT_SWNET].description, "net downward shortwave flux");

    /* net downward radiation flux [W/m2] */
    strcpy(out_metadata[OUT_R_NET].varname, "OUT_R_NET");
    strcpy(out_metadata[OUT_R_NET].long_name, "net_downward_radiation_flux_at_surface");
    strcpy(out_metadata[OUT_R_NET].units, "W m-2");
    strcpy(out_metadata[OUT_R_NET].description, "net downward radiation flux");

    /* net energy used to refreeze liquid water in snowpack [W/m2] */
    strcpy(out_metadata[OUT_RFRZ_ENERGY].varname, "OUT_RFRZ_ENERGY");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].long_name, "net_energy_used_to_refreeze_water_in_snowpack");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].units, "W m-2");
    strcpy(out_metadata[OUT_RFRZ_ENERGY].description, "net energy used to refreeze liquid water in snowpack");

    /* net upward sensible heat flux [W/m2] */
    strcpy(out_metadata[OUT_SENSIBLE].varname, "OUT_SENSIBLE");
    strcpy(out_metadata[OUT_SENSIBLE].long_name, "surface_upward_net_sensible_heat_flux");
    strcpy(out_metadata[OUT_SENSIBLE].units, "W m-2");
    strcpy(out_metadata[OUT_SENSIBLE].description, "net upward sensible heat flux");

    /* energy flux through snow pack [W/m2] */
    strcpy(out_metadata[OUT_SNOW_FLUX].varname, "OUT_SNOW_FLUX");
    strcpy(out_metadata[OUT_SNOW_FLUX].long_name, "energy_flux_through_snow_pack");
    strcpy(out_metadata[OUT_SNOW_FLUX].units, "W m-2");

    // Miscellaneous Terms
    /* "scene" aerodynamic conductance [m/s] (tiles with overstory contribute overstory conductance; others contribue surface conductance) */
    strcpy(out_metadata[OUT_AERO_COND].varname, "OUT_AERO_COND");
    strcpy(out_metadata[OUT_AERO_COND].long_name, "aerodynamic_conductance");
    strcpy(out_metadata[OUT_AERO_COND].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND].description, "scene aerodynamic conductance (tiles with overstory contribute overstory conductance; others contribue surface conductance)");

    /* surface aerodynamic conductance [m/s] */
    strcpy(out_metadata[OUT_AERO_COND1].varname, "OUT_AERO_COND1");
    strcpy(out_metadata[OUT_AERO_COND1].long_name, "aerodynamic_conductance_surface");
    strcpy(out_metadata[OUT_AERO_COND1].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND1].description, "surface aerodynamic conductance");

    /* overstory aerodynamic conductance [m/s] */
    strcpy(out_metadata[OUT_AERO_COND2].varname, "OUT_AERO_COND2");
    strcpy(out_metadata[OUT_AERO_COND2].long_name, "aerodynamic_conductance_overstory");
    strcpy(out_metadata[OUT_AERO_COND2].units, "m/s");
    strcpy(out_metadata[OUT_AERO_COND2].description, "overstory aerodynamic conductance");

    /* "scene" aerodynamic resistance [s/m] (tiles with overstory contribute overstory resistance; others contribue surface resistance)*/
    strcpy(out_metadata[OUT_AERO_RESIST].varname, "OUT_AERO_RESIST");
    strcpy(out_metadata[OUT_AERO_RESIST].long_name, "aerodynamic_resistance");
    strcpy(out_metadata[OUT_AERO_RESIST].units, "s/m");
    strcpy(out_metadata[OUT_AERO_RESIST].description, "scene aerodynamic resistance (tiles with overstory contribute overstory resistance; others contribue surface resistance)");

    /* surface aerodynamic resistance [m/s] */
    strcpy(out_metadata[OUT_AERO_RESIST1].varname, "OUT_AERO_RESIST1");
    strcpy(out_metadata[OUT_AERO_RESIST1].long_name, "aerodynamic_resistance_surface");
    strcpy(out_metadata[OUT_AERO_RESIST1].units, "s/m");
    strcpy(out_metadata[OUT_AERO_RESIST1].description, "surface aerodynamic resistance");

    /* overstory aerodynamic resistance [m/s] */
    strcpy(out_metadata[OUT_AERO_RESIST2].varname, "OUT_AERO_RESIST2");
    strcpy(out_metadata[OUT_AERO_RESIST2].long_name, "aerodynamic_resistance_overstory");
    strcpy(out_metadata[OUT_AERO_RESIST2].units, "s/m");
    strcpy(out_metadata[OUT_AERO_RESIST2].description, "overstory aerodynamic resistance");

    /* air temperature [C] */
    strcpy(out_metadata[OUT_AIR_TEMP].varname, "OUT_AIR_TEMP");
    strcpy(out_metadata[OUT_AIR_TEMP].long_name, "air_temperature");
    strcpy(out_metadata[OUT_AIR_TEMP].units, "C");
    strcpy(out_metadata[OUT_AIR_TEMP].description, "air temperature");

    /* atmospheric CO2 concentration [ppm] */
    strcpy(out_metadata[OUT_CATM].varname, "OUT_CATM");
    strcpy(out_metadata[OUT_CATM].long_name, "concentration_of_carbon_dioxide_in_air");
    strcpy(out_metadata[OUT_CATM].units, "ppm");
    strcpy(out_metadata[OUT_CATM].description, "atmospheric CO2 concentration");

    /* near-surface atmospheric density [kg/m3] */
    strcpy(out_metadata[OUT_DENSITY].varname, "OUT_DENSITY");
    strcpy(out_metadata[OUT_DENSITY].long_name, "air_density");
    strcpy(out_metadata[OUT_DENSITY].units, "kg/m3");
    strcpy(out_metadata[OUT_DENSITY].description, "near-surface atmospheric density");

    /* fractional area covered by plant canopy [fraction] */
    strcpy(out_metadata[OUT_FCANOPY].varname, "OUT_FCANOPY");
    strcpy(out_metadata[OUT_FCANOPY].long_name, "canopy_cover_area_fraction");
    strcpy(out_metadata[OUT_FCANOPY].units, "1");
    strcpy(out_metadata[OUT_FCANOPY].description, "fractional area covered by plant canopy");

    /* fraction of incoming shortwave that is direct [fraction] */
    strcpy(out_metadata[OUT_FDIR].varname, "OUT_FDIR");
    strcpy(out_metadata[OUT_FDIR].long_name, "fraction_of_incoming_shorwave_radiation_that_is_direct");
    strcpy(out_metadata[OUT_FDIR].units, "1");
    strcpy(out_metadata[OUT_FDIR].description, "fraction of incoming shortwave that is direct");

    /* leaf area index [m2/m2] */
    strcpy(out_metadata[OUT_LAI].varname, "OUT_LAI");
    strcpy(out_metadata[OUT_LAI].long_name, "leaf_area_index");
    strcpy(out_metadata[OUT_LAI].units, "1");
    strcpy(out_metadata[OUT_LAI].description, "leaf area index");

    /* incoming longwave [W/m2] */
    strcpy(out_metadata[OUT_LWDOWN].varname, "OUT_LWDOWN");
    strcpy(out_metadata[OUT_LWDOWN].long_name, "downwelling_longwave_flux_in_air");
    strcpy(out_metadata[OUT_LWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_LWDOWN].description, "incoming longwave");

    /* incoming photosynthetically active radiation [W/m2] */
    strcpy(out_metadata[OUT_PAR].varname, "OUT_PAR");
    strcpy(out_metadata[OUT_PAR].long_name, "surface_downwelling_photosynthetic_radiative_flux_in_air");
    strcpy(out_metadata[OUT_PAR].units, "W m-2");
    strcpy(out_metadata[OUT_PAR].description, "incoming photosynthetically active radiation");

    /* near surface atmospheric pressure [kPa] */
    strcpy(out_metadata[OUT_PRESSURE].varname, "OUT_PRESSURE");
    strcpy(out_metadata[OUT_PRESSURE].long_name, "surface_air_pressure");
    strcpy(out_metadata[OUT_PRESSURE].units, "kPa");
    strcpy(out_metadata[OUT_PRESSURE].description, "near surface atmospheric pressure");

    /* specific humidity [kg/kg] */
    strcpy(out_metadata[OUT_QAIR].varname, "OUT_QAIR");
    strcpy(out_metadata[OUT_QAIR].long_name, "specific_humidity");
    strcpy(out_metadata[OUT_QAIR].units, "1");
    strcpy(out_metadata[OUT_QAIR].description, "specific humidity");

    /* relative humidity [fraction]*/
    strcpy(out_metadata[OUT_REL_HUMID].varname, "OUT_REL_HUMID");
    strcpy(out_metadata[OUT_REL_HUMID].long_name, "relative_humidity");
    strcpy(out_metadata[OUT_REL_HUMID].units, "1");
    strcpy(out_metadata[OUT_REL_HUMID].description, "relative humidity");

    /* incoming shortwave [W/m2] */
    strcpy(out_metadata[OUT_SWDOWN].varname, "OUT_SWDOWN");
    strcpy(out_metadata[OUT_SWDOWN].long_name, "incoming shortwave");
    strcpy(out_metadata[OUT_SWDOWN].units, "W m-2");
    strcpy(out_metadata[OUT_SWDOWN].description, "incoming shortwave");

    /* surface conductance [m/s] */
    strcpy(out_metadata[OUT_SURF_COND].varname, "OUT_SURF_COND");
    strcpy(out_metadata[OUT_SURF_COND].long_name, "downwelling_shortwave_flux_in_air");
    strcpy(out_metadata[OUT_SURF_COND].units, "unset");
    strcpy(out_metadata[OUT_SURF_COND].description, "something");

    /* near surface vapor pressure [kPa] */
    strcpy(out_metadata[OUT_VP].varname, "OUT_VP");
    strcpy(out_metadata[OUT_VP].long_name, "water_vapor_pressure");
    strcpy(out_metadata[OUT_VP].units, "kPa");
    strcpy(out_metadata[OUT_VP].description, "near surface vapor pressure");

    /* near surface vapor pressure deficit [kPa] */
    strcpy(out_metadata[OUT_VPD].varname, "OUT_VPD");
    strcpy(out_metadata[OUT_VPD].long_name, "water_vapor_saturation_deficit");
    strcpy(out_metadata[OUT_VPD].units, "kPa");
    strcpy(out_metadata[OUT_VPD].description, "near surface vapor pressure deficit");

    /* near surface wind speed [m/s] */
    strcpy(out_metadata[OUT_WIND].varname, "OUT_WIND");
    strcpy(out_metadata[OUT_WIND].long_name, "wind_speed");
    strcpy(out_metadata[OUT_WIND].units, "unset");
    strcpy(out_metadata[OUT_WIND].description, "near surface wind speed");

    // Carbon-cycling Terms
    /* absorbed PAR [W/m2] */
    strcpy(out_metadata[OUT_APAR].varname, "OUT_APAR");
    strcpy(out_metadata[OUT_APAR].long_name, "absorbed_surface_diffuse_downwelling_photosynthetic_radiative_flux");
    strcpy(out_metadata[OUT_APAR].units, "W m-2");
    strcpy(out_metadata[OUT_APAR].description, "absorbed PAR");

    /* gross primary productivity [g C/m2s] */
    strcpy(out_metadata[OUT_GPP].varname, "OUT_GPP");
    strcpy(out_metadata[OUT_GPP].long_name, "gross_primary_productivity_of_biomass_expressed_as_carbon");
    strcpy(out_metadata[OUT_GPP].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_GPP].description, "gross primary productivity");

    /* autotrophic respiration [g C/m2s] */
    strcpy(out_metadata[OUT_RAUT].varname, "OUT_RAUT");
    strcpy(out_metadata[OUT_RAUT].long_name, "autotrophic_respiration_carbon_flux");
    strcpy(out_metadata[OUT_RAUT].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_RAUT].description, "autotrophic respiration");

    /* net primary productivity [g C/m2s] */
    strcpy(out_metadata[OUT_NPP].varname, "OUT_NPP");
    strcpy(out_metadata[OUT_NPP].long_name, "net_primary_productivity_of_biomass_expressed_as_carbon");
    strcpy(out_metadata[OUT_NPP].units, "g m-2 s-1");
    strcpy(out_metadata[OUT_NPP].description, "et primary productivity");

    /* flux of carbon from living biomass into soil [g C/m2d] */
    strcpy(out_metadata[OUT_LITTERFALL].varname, "OUT_LITTERFALL");
    strcpy(out_metadata[OUT_LITTERFALL].long_name, "carbon_mass_flux_into_soil_from_litter");
    strcpy(out_metadata[OUT_LITTERFALL].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_LITTERFALL].description, "flux of carbon from living biomass into soil");

    /* heterotrophic respiration [g C/m2d] */
    strcpy(out_metadata[OUT_RHET].varname, "OUT_RHET");
    strcpy(out_metadata[OUT_RHET].long_name, "unset");
    strcpy(out_metadata[OUT_RHET].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_RHET].description, "something");

    /* net ecosystem exchange [g C/m2d] */
    strcpy(out_metadata[OUT_NEE].varname, "OUT_NEE");
    strcpy(out_metadata[OUT_NEE].long_name, "net_ecosystem_exhanged_expressed_as_carbon");
    strcpy(out_metadata[OUT_NEE].units, "g m-2 d-1");
    strcpy(out_metadata[OUT_NEE].description, "net ecosystem exchange");

    /* litter pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CLITTER].varname, "OUT_CLITTER");
    strcpy(out_metadata[OUT_CLITTER].long_name, "litter_carbon_content");
    strcpy(out_metadata[OUT_CLITTER].units, "g m-2");
    strcpy(out_metadata[OUT_CLITTER].description, "litter pool carbon density");

    /* intermediate pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CINTER].varname, "OUT_CINTER");
    strcpy(out_metadata[OUT_CINTER].long_name, "intermediate_pool_carbon_content");
    strcpy(out_metadata[OUT_CINTER].units, "g m-2");
    strcpy(out_metadata[OUT_CINTER].description, "intermediate pool carbon density");

    /* slow pool carbon density [g C/m2] */
    strcpy(out_metadata[OUT_CSLOW].varname, "OUT_CSLOW");
    strcpy(out_metadata[OUT_CSLOW].long_name, "slow_pool_carbon_content");
    strcpy(out_metadata[OUT_CSLOW].units, "g m-2");
    strcpy(out_metadata[OUT_CSLOW].description, "slow pool carbon density");

    // Band-specific quantities
    /* net sensible heat flux advected to snow pack [W/m2] */
    strcpy(out_metadata[OUT_ADV_SENS_BAND].varname, "OUT_ADV_SENS_BAND");
    strcpy(out_metadata[OUT_ADV_SENS_BAND].long_name, out_metadata[OUT_ADV_SENS].long_name);
    strcpy(out_metadata[OUT_ADV_SENS_BAND].units, out_metadata[OUT_ADV_SENS].units);
    strcpy(out_metadata[OUT_ADV_SENS_BAND].description, out_metadata[OUT_ADV_SENS].description);

    strcpy(out_metadata[OUT_ADV_SENS].varname, "OUT_ADV_SENS");
    strcpy(out_metadata[OUT_ADV_SENS].long_name, "net_sensible_heat_flux_to_snow_pack");
    strcpy(out_metadata[OUT_ADV_SENS].units, "W m-2");
    strcpy(out_metadata[OUT_ADV_SENS].description, "net sensible heat advected to snow pack");

    /* advected energy [W/m2] */
    strcpy(out_metadata[OUT_ADVECTION_BAND].varname, "OUT_ADVECTION_BAND");
    strcpy(out_metadata[OUT_ADVECTION_BAND].long_name, out_metadata[OUT_ADVECTION].long_name);
    strcpy(out_metadata[OUT_ADVECTION_BAND].units, out_metadata[OUT_ADVECTION].units);
    strcpy(out_metadata[OUT_ADVECTION_BAND].description, out_metadata[OUT_ADVECTION].description);

    /* albedo [fraction] */
    strcpy(out_metadata[OUT_ALBEDO_BAND].varname, "OUT_ALBEDO_BAND");
    strcpy(out_metadata[OUT_ALBEDO_BAND].long_name, out_metadata[OUT_ALBEDO].long_name);
    strcpy(out_metadata[OUT_ALBEDO_BAND].units, out_metadata[OUT_ALBEDO].units);
    strcpy(out_metadata[OUT_ALBEDO_BAND].description, out_metadata[OUT_ALBEDO].description);

    /* change in cold content in snow pack [W/m2] */
    strcpy(out_metadata[OUT_DELTACC_BAND].varname, "OUT_DELTACC_BAND");
    strcpy(out_metadata[OUT_DELTACC_BAND].long_name, out_metadata[OUT_DELTACC].long_name);
    strcpy(out_metadata[OUT_DELTACC_BAND].units, out_metadata[OUT_DELTACC].units);
    strcpy(out_metadata[OUT_DELTACC_BAND].description, out_metadata[OUT_DELTACC].description);

    /* net heat flux into ground [W/m2] */
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].varname, "OUT_GRND_FLUX_BAND");
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].long_name, out_metadata[OUT_GRND_FLUX].long_name);
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].units, out_metadata[OUT_GRND_FLUX].units);
    strcpy(out_metadata[OUT_GRND_FLUX_BAND].description, out_metadata[OUT_GRND_FLUX].description);

    /* incoming longwave flux at surface (under veg) [W/m2] */
    strcpy(out_metadata[OUT_IN_LONG_BAND].varname, "OUT_IN_LONG_BAND");
    strcpy(out_metadata[OUT_IN_LONG_BAND].long_name, out_metadata[OUT_IN_LONG].long_name);
    strcpy(out_metadata[OUT_IN_LONG_BAND].units, out_metadata[OUT_IN_LONG].units);
    strcpy(out_metadata[OUT_IN_LONG_BAND].description, out_metadata[OUT_IN_LONG].description);

    /* net upward latent heat flux [W/m2] */
    strcpy(out_metadata[OUT_LATENT_BAND].varname, "OUT_LATENT_BAND");
    strcpy(out_metadata[OUT_LATENT_BAND].long_name, out_metadata[OUT_LATENT].long_name);
    strcpy(out_metadata[OUT_LATENT_BAND].units, out_metadata[OUT_LATENT].units);
    strcpy(out_metadata[OUT_LATENT_BAND].description, out_metadata[OUT_LATENT].description);

    /* net upward latent heat flux from sublimation [W/m2] */
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].varname, "OUT_LATENT_SUB_BAND");
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].long_name, out_metadata[OUT_LATENT_SUB].long_name);
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].units, out_metadata[OUT_LATENT_SUB].units);
    strcpy(out_metadata[OUT_LATENT_SUB_BAND].description, out_metadata[OUT_LATENT_SUB].description);

    /* energy of fusion (melting) [W/m2] */
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].varname, "OUT_MELT_ENERGY_BAND");
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].long_name, out_metadata[OUT_MELT_ENERGY].long_name);
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].units, out_metadata[OUT_MELT_ENERGY].units);
    strcpy(out_metadata[OUT_MELT_ENERGY_BAND].description, out_metadata[OUT_MELT_ENERGY].description);

    /* net downward longwave flux [W/m2] */
    strcpy(out_metadata[OUT_LWNET_BAND].varname, "OUT_LWNET_BAND");
    strcpy(out_metadata[OUT_LWNET_BAND].long_name, out_metadata[OUT_LWNET].long_name);
    strcpy(out_metadata[OUT_LWNET_BAND].units, out_metadata[OUT_LWNET].units);
    strcpy(out_metadata[OUT_LWNET_BAND].description, out_metadata[OUT_LWNET].description);

    /* net downward shortwave flux [W/m2] */
    strcpy(out_metadata[OUT_SWNET_BAND].varname, "OUT_SWNET_BAND");
    strcpy(out_metadata[OUT_SWNET_BAND].long_name, out_metadata[OUT_SWNET].long_name);
    strcpy(out_metadata[OUT_SWNET_BAND].units, out_metadata[OUT_SWNET].units);
    strcpy(out_metadata[OUT_SWNET_BAND].description, out_metadata[OUT_SWNET].description);

    /* net energy used to refreeze liquid water in snowpack [W/m2] */
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].varname, "OUT_RFRZ_ENERGY_BAND");
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].long_name, out_metadata[OUT_RFRZ_ENERGY].long_name);
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].units, out_metadata[OUT_RFRZ_ENERGY].units);
    strcpy(out_metadata[OUT_RFRZ_ENERGY_BAND].description, out_metadata[OUT_RFRZ_ENERGY].description);

    /* net upward sensible heat flux [W/m2] */
    strcpy(out_metadata[OUT_SENSIBLE_BAND].varname, "OUT_SENSIBLE_BAND");
    strcpy(out_metadata[OUT_SENSIBLE_BAND].long_name, out_metadata[OUT_SENSIBLE].long_name);
    strcpy(out_metadata[OUT_SENSIBLE_BAND].units, out_metadata[OUT_SENSIBLE].units);
    strcpy(out_metadata[OUT_SENSIBLE_BAND].description, out_metadata[OUT_SENSIBLE].description);

    /* snow interception storage in canopy [mm] */
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].varname, "OUT_SNOW_CANOPY_BAND");
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].long_name, out_metadata[OUT_SNOW_CANOPY].long_name);
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].units, out_metadata[OUT_SNOW_CANOPY].units);
    strcpy(out_metadata[OUT_SNOW_CANOPY_BAND].description, out_metadata[OUT_SNOW_CANOPY].description);

    /* fractional area of snow cover [fraction] */
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].varname, "OUT_SNOW_COVER_BAND");
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].long_name, out_metadata[OUT_SNOW_COVER].long_name);
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].units, out_metadata[OUT_SNOW_COVER].units);
    strcpy(out_metadata[OUT_SNOW_COVER_BAND].description, out_metadata[OUT_SNOW_COVER].description);

    /* depth of snow pack [cm] */
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].varname, "OUT_SNOW_DEPTH_BAND");
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].long_name, out_metadata[OUT_SNOW_DEPTH].long_name);
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].units, out_metadata[OUT_SNOW_DEPTH].units);
    strcpy(out_metadata[OUT_SNOW_DEPTH_BAND].description, out_metadata[OUT_SNOW_DEPTH].description);

    /* energy flux through snow pack [W/m2] */
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].varname, "OUT_SNOW_FLUX_BAND");
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].long_name, out_metadata[OUT_SNOW_FLUX].long_name);
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].units, out_metadata[OUT_SNOW_FLUX].units);
    strcpy(out_metadata[OUT_SNOW_FLUX_BAND].description, out_metadata[OUT_SNOW_FLUX].description);

    /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].varname, "OUT_SNOW_MELT_BAND");
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].long_name, out_metadata[OUT_SNOW_MELT].long_name);
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].units, out_metadata[OUT_SNOW_MELT].units);
    strcpy(out_metadata[OUT_SNOW_MELT_BAND].description, out_metadata[OUT_SNOW_MELT].description);

    /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].varname, "OUT_SNOW_PACKT_BAND");
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].long_name, out_metadata[OUT_SNOW_PACK_TEMP].long_name);
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].units, out_metadata[OUT_SNOW_PACK_TEMP].units);
    strcpy(out_metadata[OUT_SNOW_PACKT_BAND].description, out_metadata[OUT_SNOW_PACK_TEMP].description);

    /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].varname, "OUT_SNOW_SURFT_BAND");
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].long_name, out_metadata[OUT_SNOW_SURF_TEMP].long_name);
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].units, out_metadata[OUT_SNOW_SURF_TEMP].units);
    strcpy(out_metadata[OUT_SNOW_SURFT_BAND].description, out_metadata[OUT_SNOW_SURF_TEMP].description);

    /* snow water equivalent in snow pack [mm] */
    strcpy(out_metadata[OUT_SWE_BAND].varname, "OUT_SWE_BAND");
    strcpy(out_metadata[OUT_SWE_BAND].long_name, out_metadata[OUT_SWE].long_name);
    strcpy(out_metadata[OUT_SWE_BAND].units, out_metadata[OUT_SWE].units);
    strcpy(out_metadata[OUT_SWE_BAND].description, out_metadata[OUT_SWE].description);

    if (options.FROZEN_SOIL) {
        out_metadata[OUT_FDEPTH].nelem = MAX_FRONTS;
        out_metadata[OUT_TDEPTH].nelem = MAX_FRONTS;
    }

    out_metadata[OUT_SMLIQFRAC].nelem = options.Nlayer;
    out_metadata[OUT_SMFROZFRAC].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_ICE].nelem = options.Nlayer;
    out_metadata[OUT_SOIL_LIQ].nelem = options.Nlayer;
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

/******************************************************************************
 * @brief    This routine creates the list of output data.
 *****************************************************************************/
double **
create_outdata()
{
    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t               i;
    size_t               j;
    double             **out_data;

    out_data = calloc(N_OUTVAR_TYPES, sizeof(*out_data));

    // Allocate space for data
    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        out_data[i] = calloc(out_metadata[i].nelem, sizeof(*(out_data[i])));

        // initialize data member
        for (j = 0; j < out_metadata[i].nelem; j++) {
            out_data[i][j] = 0;
        }
    }

    return out_data;
}

/******************************************************************************
 * @brief    This routine creates the list of output streams.
 *****************************************************************************/
void
setup_stream(stream_struct      *stream,
             stream_file_struct *stream_file,
             size_t              nvars,
             unsigned int        nextagg)
{

    size_t i;

    stream->nvars = nvars;
    stream->nextagg = nextagg;
    stream->counter = 0;

    stream_file->file_format = UNSET_FILE_FORMAT;
    stream_file->compress = false;
    stream_file->output_steps_per_day = 1;
    stream_file->output_steps_per_day = 0;

    // Allocate stream members of shape [nvars]
    stream->varid = calloc(nvars, sizeof(*(stream->varid)));
    if (stream->varid == NULL) {
        log_err("Memory allocation error in setup_stream().");
    }
    stream->aggtype = calloc(nvars, sizeof(*(stream->aggtype)));
    if (stream->aggtype == NULL) {
        log_err("Memory allocation error in setup_stream().");
    }
    stream->aggdata = calloc(nvars, sizeof(*(stream->aggdata)));
    if (stream->aggdata == NULL) {
        log_err("Memory allocation error in setup_stream().");
    }

    // Allocate stream_file members of shape [nvars]
    stream_file->type = calloc(nvars, sizeof(*(stream_file->type)));
    if (stream_file->type == NULL) {
        log_err("Memory allocation error in setup_stream_file().");
    }
    stream_file->mult = calloc(nvars, sizeof(*(stream_file->mult)));
    if (stream_file->mult == NULL) {
        log_err("Memory allocation error in setup_stream_file().");
    }
    // Question: do we have to dynamically allocate the length of each string
    stream_file->format = calloc(nvars, sizeof(*(stream_file->format)));
    if (stream_file->format == NULL) {
        log_err("Memory allocation error in setup_stream().");
    }
    for (i = 0; i < nvars; i++) {
        stream_file->format[i] = calloc(MAXSTRING, sizeof(*(stream_file->format[i])));
        if (stream_file->format[i] == NULL) {
            log_err("Memory allocation error in setup_stream().");
        }
    }
}

/******************************************************************************
 * @brief   This routine resets an output stream
 *****************************************************************************/
void
reset_stream(stream_struct *stream) {

    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t  i;
    size_t  j;
    size_t  v;

    stream->counter = 0;

    for (i = 0; i < stream->nvars; i++) {
        v = stream->varid[i];
        for (j = 0; j < out_metadata[v].nelem; j++) {
            stream->aggdata[i][j][0] = 0.;
        }
    }
}


/******************************************************************************
 * @brief   This routine sets the default aggregation type for variables in an
            output stream
 *****************************************************************************/
unsigned int
get_default_outvar_aggtype(unsigned int varid)
{

    unsigned int agg_type;

    switch (varid) {
        // AGG_TYPE_END
        case OUT_ASAT:
        case OUT_LAKE_AREA_FRAC:
        case OUT_LAKE_DEPTH:
        case OUT_LAKE_ICE:
        case OUT_LAKE_ICE_FRACT:
        case OUT_LAKE_ICE_HEIGHT:
        case OUT_LAKE_MOIST:
        case OUT_LAKE_SURF_AREA:
        case OUT_LAKE_SWE:
        case OUT_LAKE_SWE_V:
        case OUT_LAKE_VOLUME:
        case OUT_ROOTMOIST:
        case OUT_SMFROZFRAC:
        case OUT_SMLIQFRAC:
        case OUT_SNOW_CANOPY:
        case OUT_SNOW_COVER:
        case OUT_SNOW_DEPTH:
        case OUT_SOIL_ICE:
        case OUT_SOIL_LIQ:
        case OUT_SOIL_MOIST:
        case OUT_SOIL_WET:
        case OUT_SURFSTOR:
        case OUT_SURF_FROST_FRAC:
        case OUT_SWE:
        case OUT_WDEW:
        case OUT_ZWT:
        case OUT_ZWT_LUMPED:
        case OUT_SNOW_CANOPY_BAND:
        case OUT_SNOW_COVER_BAND:
        case OUT_SNOW_DEPTH_BAND:
        case OUT_SWE_BAND:
            agg_type = AGG_TYPE_END;
            break;
        // AGG_TYPE_SUM
        case OUT_BASEFLOW:
        case OUT_DELINTERCEPT:
        case OUT_DELSOILMOIST:
        case OUT_DELSWE:
        case OUT_DELSURFSTOR:
        case OUT_EVAP:
        case OUT_EVAP_BARE:
        case OUT_EVAP_CANOP:
        case OUT_INFLOW:
        case OUT_LAKE_BF_IN:
        case OUT_LAKE_BF_IN_V:
        case OUT_LAKE_BF_OUT:
        case OUT_LAKE_BF_OUT_V:
        case OUT_LAKE_CHAN_IN:
        case OUT_LAKE_CHAN_IN_V:
        case OUT_LAKE_CHAN_OUT:
        case OUT_LAKE_CHAN_OUT_V:
        case OUT_LAKE_DSTOR:
        case OUT_LAKE_DSTOR_V:
        case OUT_LAKE_DSWE:
        case OUT_LAKE_DSWE_V:
        case OUT_LAKE_EVAP:
        case OUT_LAKE_EVAP_V:
        case OUT_LAKE_PREC_V:
        case OUT_LAKE_RCHRG:
        case OUT_LAKE_RCHRG_V:
        case OUT_LAKE_RO_IN:
        case OUT_LAKE_RO_IN_V:
        case OUT_LAKE_VAPFLX:
        case OUT_LAKE_VAPFLX_V:
        case OUT_PET:
        case OUT_PREC:
        case OUT_RAINF:
        case OUT_REFREEZE:
        case OUT_RUNOFF:
        case OUT_SNOWF:
        case OUT_SUB_BLOWING:
        case OUT_SUB_CANOP:
        case OUT_SUB_SNOW:
        case OUT_SUB_SURFACE:
        case OUT_TRANSP_VEG:
        case OUT_DELTACC_BAND:
        case OUT_SNOW_MELT:
        case OUT_SNOWT_FBFLAG:
        case OUT_SOILT_FBFLAG:
        case OUT_SURFT_FBFLAG:
        case OUT_TCAN_FBFLAG:
        case OUT_TFOL_FBFLAG:
            agg_type = AGG_TYPE_SUM;
            break;
        default:
            agg_type = AGG_TYPE_AVG;
    }
    return agg_type;
}

/******************************************************************************
 * @brief    This routine updates the output information for a given output
 *           variable.
 *****************************************************************************/
void
set_output_var(stream_struct      *stream,
               stream_file_struct *out_data_file,
               char               *varname,
               size_t              varnum,
               char               *format,
               unsigned short int  type,
               double              mult)
{

    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    int varid;
    int found = false;

    for (varid = 0; varid < N_OUTVAR_TYPES; varid++) {
        if (strcmp(out_metadata[varid].varname, varname) == 0) {
            found = true;
            break;
        }
    }
    if (!found) {
        log_err("set_output_var: \"%s\" was not found in the list of "
                "supported output variable names.  Please use the exact name "
                "listed in vic_driver_shared.h.", varname);
    }

    // Set stream_file_struct members
    if (strcmp(format, "*") != 0) {
        strcpy(out_data_file->format[varnum], format);
    }
    else {
        strcpy(out_data_file->format[varnum], "%.4f");
    }
    if (type != 0) {
        out_data_file->type[varnum] = type;
    }
    else {
        out_data_file->type[varnum] = OUT_TYPE_FLOAT;
    }
    if (mult != 0) {
        out_data_file->mult[varnum] = mult;
    }
    else {
        out_data_file->mult[varnum] = 1.;
    }

    // Set stream members
    stream->varid[varnum] = varid;

    // TODO: set aggtype
    stream->aggtype[varnum] = get_default_outvar_aggtype(varid);

}

/******************************************************************************
 * @brief    This routine frees the memory in the out_data_files array.
 *****************************************************************************/
void
free_out_data_files(stream_file_struct **out_data_files)
{
    extern option_struct options;
    size_t               i;

    // free output
    for (i = 0; i < options.Noutstreams; i++) {
        free((char*) out_data_files[i]->type);
        free((char*) out_data_files[i]->mult);
        free((char*) out_data_files[i]->format);
    }
    free((char*) out_data_files);
}

/******************************************************************************
 * @brief    This routine frees the memory in the streams array.
 *****************************************************************************/
void
free_out_data_streams(stream_struct **streams)
{
    extern option_struct options;
    size_t               i;

    // free output
    for (i = 0; i < options.Noutstreams; i++) {
        free((char*) streams[i]->varid);
        free((char*) streams[i]->aggtype);
        // TODO: Free all 3 dimensions
        free((char*) streams[i]->aggdata);
    }
    free((char*) streams);
}

/******************************************************************************
 * @brief    This routine frees the memory in the out_data array.
 *****************************************************************************/
void
free_out_data(double **out_data)
{
    extern out_metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t  varid;

    if (out_data == NULL) {
        return;
    }

    for (varid = 0; varid < N_OUTVAR_TYPES; varid++) {
        free(out_data[varid]);
    }

    free((char*) (*out_data));
}

/******************************************************************************
 * @brief    Validate the streams settings.
 *****************************************************************************/
void
validate_stream_settings(stream_file_struct **out_data_file)
{

    extern global_param_struct global_param;
    extern option_struct options;

    size_t i;

    for (i = 0; i < options.Noutstreams; i++) {

        // Validate the output step
        if (out_data_file[i]->output_steps_per_day == 0) {
            out_data_file[i]->output_steps_per_day = global_param.model_steps_per_day;
        }
        if (out_data_file[i]->output_steps_per_day > global_param.model_steps_per_day) {
            log_err("[stream=%zu] Invalid value for OUTPUT_STEPS_PER_DAY (%zu). "
                    "OUTPUT_STEPS_PER_DAY must be <= MODEL_STEPS_PER_DAY (%zu)",
                    i, out_data_file[i]->output_steps_per_day,
                    global_param.model_steps_per_day);
        }
        else if (global_param.model_steps_per_day % out_data_file[i]->output_steps_per_day != 0) {
            log_err("[stream=%zu] Invalid value for OUTPUT_STEPS_PER_DAY (%zu). "
                    "MODEL_STEPS_PER_DAY (%zu) must be a multiple of "
                    "OUTPUT_STEPS_PER_DAY.", i,
                    out_data_file[i]->output_steps_per_day,
                    global_param.model_steps_per_day);
        }
        else if (out_data_file[i]->output_steps_per_day != 1 &&
                 out_data_file[i]->output_steps_per_day < MIN_SUBDAILY_STEPS_PER_DAY) {
            log_err("[stream=%zu] The specified number of output steps per day "
                    "(%zu) > 1 and < the minimum number of subdaily steps per "
                    "day (%d). Make sure that the global file defines "
                    "OUTPUT_STEPS_PER_DAY of at least (%d).", i,
                    global_param.model_steps_per_day,
                    MIN_SUBDAILY_STEPS_PER_DAY, MIN_SUBDAILY_STEPS_PER_DAY);
        }
        else if (out_data_file[i]->output_steps_per_day > MAX_SUBDAILY_STEPS_PER_DAY) {
            log_err("[stream=%zu] The specified number of model steps per day "
                    "(%zu) > the the maximum number of subdaily steps per day "
                    "(%d). Make sure that the global file defines "
                    "MODEL_STEPS_PER_DAY of at most (%d).", i,
                    global_param.model_steps_per_day,
                    MAX_SUBDAILY_STEPS_PER_DAY, MAX_SUBDAILY_STEPS_PER_DAY);
        }
        else {
            out_data_file[i]->out_dt = SEC_PER_DAY / (double) out_data_file[i]->output_steps_per_day;
        }
    }
}
