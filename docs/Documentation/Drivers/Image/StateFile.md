# VIC Model State File - Image Driver

VIC can save the hydrologic state from any point in the simulation (usually the final state) to a file for the purpose of re-starting the simulation later (as an initial state file). This is useful for simulations that require lengthy spin-up periods or ensemble methods. The initial state file is not required; if it is not specified, VIC will use a default starting condition, and will take initial soil moisture contents from the values specified in the [parameters file](Params.md). An initial state file may be prepared simply by running VIC with the necessary state file options in the [global parameter file](GlobalParam.md#DefineStateFiles).

The model state file contains all information needed by the VIC model to "warm"-start a simulation (i.e. start from "realistic" conditions, or re-start a simulation exactly where the model previously stopped). To read an initial state file, or to save a "final" state file, the appropriate options should be set in the [global parameter file](GlobalParam.md#DefineStateFiles).

The timestamp of the state file represents the instantaneous time for which the values of the state variables are valid. This corresponds to the end of the time interval after which the state file is written out and the beginning of the time interval for which the model is started. For example, if the MODEL_STEPS_PER_DAY is 24 (hourly) and the last time step for which the model is run is 1999-09-20 23:00:00, the the state file will be stamped 1999-09-21 00:00:00. This state file can then be used to restart a model run whose starting time will be 1999-09-21 00:00:00.

The state file in image driver contains the same variables as in [classic driver state file](../Classic/StateFile.md), but is in netCDF format. The following is a detailed description of the dimensions and variables of the netCDF state file.

* * *

## netCDF State File Dimensions

The netCDF state file dimensions include:

| Dimension name | Description                                      |
|----------------|--------------------------------------------------|
| lon            | Number of longitudes                             |
| lat            | Number of latitudes                              |
| nlayer         | Number of soil layers                            |
| soil_node      | Number of soil thermal nodes                     |
| veg_class      | Number of vegetation types (including bare soil) |
| snow_band      | Number of snow bands                             |
| frost_area     | Number of frost areas                            |

* * *

## netCDF State File Variables - Dimensions, soil layers and thermal nodes

The following variables define the basic model information, including grid cell lat and lon, vegetation classes, snow bands, soil layers and thermal nodes and frost area.

| Variable   name | Dimension  | Type   | Description                                                                                                                 |
|-----------------|------------|--------|-----------------------------------------------------------------------------------------------------------------------------|
| lat             | lat        | double | List of latitudes                                                                                                           |
| lon             | lon        | double | List of longitudes                                                                                                          |
| veg_class       | veg_class  | int    | Vegetation class indices                                                                                                    |
| snow_band       | snow_band  | int    | Snow band indices                                                                                                           |
| layer           | nlayer     | int    | Soil layer indices                                                                                                          |
| frost_area      | frost_area | int    | Frost area indices                                                                                                          |
| dz_node         | [soil_node, lat, lon]  | double | Distances between soil thermal   nodes [m]                                                                                  |
| node_depth      | [soil_node, lat, lon]  | double | Depth   from surface of each soil thermal node (first node should have a depth of 0m   indicating it is at the surface) [m] |

* * *

## netCDF State File Variables - state variables

The following variables contain information about the storage of moisture and thermal state for each vegetation and snow band tile within each grid cell. If the model is being run with distributed precipitation, the wet and dry fractions are averaged before the model state is stored and the model is always initialized with a mu value of 1.



| State   variable name       | Dimension                                             | Type            | Description                                                                                                                                          |
|-----------------------------|-------------------------------------------------------|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| STATE_SOIL_MOISTURE         | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Soil   total moisture contents including ice for each soil layer [mm]                                                                                |
| STATE_SOIL_ICE              | [veg_class, snow_band, nlayer,   frost_area, at, lon] | double          | Soil ice   content for each soil layer [mm]                                                                                                          |
| STATE_CANOPY_WATER          | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Amount   of water stored in the vegetation canopy [mm]. Not defined for bare soil                                                                    |
| STATE_SNOW_AGE              | [veg_class, snow_band, nlayer,   lat, lon]            | model_time_step | Number of   model time steps since the last new snow                                                                                                 |
| STATE_SNOW_MELT_STATE       | [veg_class, snow_band, nlayer,   lat, lon]            | int             | flag to   indicate whether snowpack is in accumulation or melting phase (1 melting, 0   not melting)                                                 |
| STATE_SNOW_COVERAGE         | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Fraction   of grid cell area covered by snow [1]                                                                                                     |
| STATE_SNOW_WATER_EQUIVALENT | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snow   water equivalent [m]                                                                                                                          |
| STATE_SNOW_SURF_TEMP        | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snow   surface layer temperature [C]                                                                                                                 |
| STATE_SNOW_SURF_WATER       | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Liquid   water content of the snow surface layer [m]                                                                                                 |
| STATE_SNOW_PACK_TEMP        | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snow pack   layer temperature [C]                                                                                                                    |
| STATE_SNOW_PACK_WATER       | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Liquid   water content of the snow pack layer [m]                                                                                                    |
| STATE_SNOW_DENSITY          | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snowpack   density [kg/m3]                                                                                                                           |
| STATE_SNOW_COLD_CONTENT     | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snowpack   cold content [J/m2]                                                                                                                       |
| STATE_SNOW_CANOPY           | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Snow   interception storage in canopy [m]                                                                                                            |
| STATE_SOIL_NODE_TEMP        | [veg_class, snow_band,   soil_node, nlayer, lat, lon] | double          | Soil   temperature of each soil thermal node [C]                                                                                                     |
| STATE_FOLIAGE_TEMPERATURE   | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Overstory vegetaion temperature   [C]                                                                                                                |
| STATE_ENERGY_LONGUNDEROUT   | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Outgoing longwave flux from   understory vegetation [W/m2] (Note: this is a flux variable temporarily saved   in state file to ensure exact restart) |
| STATE_ENERGY_SNOW_FLUX      | [veg_class, snow_band, nlayer,   lat, lon]            | double          | Thermal flux through snowpack   [W/m2] (Note: this is a flux variable temporarily saved in state file to   ensure exact restart)                     |

## netCDF State File Variables - RVIC-Routing (only when RVIC extension is activated in the makefile before compiling)

| State   variable name       | Dimension                                             | Type            | Description                                                                                                                                          |
|-----------------------------|-------------------------------------------------------|-----------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| STATE_ROUT_RING             | [routing_timestep, outlet]                            | double          | Routing ring   unit hydrographs in the routing ring                                                                                                  |

* * *

## Carbon Information

If `CARBON=TRUE` in the [global parameter file](GlobalParam.md), the following variables will appear in the state file.

| State   variable name | Type   | Description                                   |
|-----------------------|--------|-----------------------------------------------|
| STATE_ANNUALNPP       | double | Running   total annual NPP [gC/m2]            |
| STATE_ANNUALNPPPREV   | double | Total   annual NPP from previous year [gC/m2] |
| STATE_CLITTER         | double | Carbon   storage in litter pool [gC/m2]       |
| STATE_CINTER          | double | Carbon   storage in intermediate pool [gC/m2] |
| STATE_CSLOW           | double | Carbon   storage in slow pool [gC/m2]         |

* * *

## Lake Information (only when LAKES are turned on in the [global parameter file](GlobalParam.md#DefineStateFiles))

| State   variable name            | Type              | Description                                                                                                      |
|----------------------------------|-------------------|------------------------------------------------------------------------------------------------------------------|
| STATE_LAKE_SOIL_MOISTURE         | double            | Soil   moisture below lake [mm]                                                                                  |
| STATE_LAKE_SOIL_ICE              | double            | Soil ice   content below lake [mm]                                                                               |
| STATE_LAKE_CLITTER               | double            | Carbon   storage in litter pool below lake [gC/m2] (Note: this variable only appears   if `CARBON=TRUE`)         |
| STATE_LAKE_CINTER                | double            | Carbon   storage in intermediate pool below lake [gC/m2] (Note: this variable only   appears if `CARBON=TRUE`)   |
| STATE_LAKE_CSLOW                 | double            | Carbon   storage in slow pool below lake [gC/m2] (Note: this variable only appears if   `CARBON=TRUE`)           |
| STATE_LAKE_SNOW_AGE              | model time   step | Number of   model time steps since the last new snow on lake ice                                                 |
| STATE_LAKE_SNOW_MELT_STATE       | int               | Flag to   indicate whether snowpack is in accumulation or melting phase on lake ice (1   melting, 0 not melting) |
| STATE_LAKE_SNOW_COVERAGE         | 1                 | Fraction   of grid cell area covered by snow on lake ice                                                         |
| STATE_LAKE_SNOW_WATER_EQUIVALENT | double            | Lake snow   water equivalent on lake ice [m]                                                                     |
| STATE_LAKE_SNOW_SURF_TEMP        | double            | Snow   surface layer temperature on lake ice [C]                                                                 |
| STATE_LAKE_SNOW_SURF_WATER       | double            | Liquid   water content of the snow surface layer on lake ice [m]                                                 |
| STATE_LAKE_SNOW_PACK_TEMP        | double            | Snow pack   layer temperature on lake ice [C]                                                                    |
| STATE_LAKE_SNOW_PACK_WATER       | double            | Liquid   water content of the snow surface layer on lake ice [m]                                                 |
| STATE_LAKE_SNOW_DENSITY          | double            | Snowpack   density on lake ice [kg/m3]                                                                           |
| STATE_LAKE_SNOW_COLD_CONTENT     | double            | Snowpack   cold content on lake ice [J/m2]                                                                       |
| STATE_LAKE_SNOW_CANOPY           | double            | Snow   interception storage in canopy on lake ice [m]                                                            |
| STATE_LAKE_SOIL_NODE_TEMP        | double            | Soil   temperature of each soil thermal node below lake[C]                                                       |
| STATE_LAKE_ACTIVE_LAYERS         | int               | Number of nodes whose corresponding layers   currently contain water                                             |
| STATE_LAKE_LAYER_DZ              | double            | Vertical   thickness of all horizontal lake water layers below the surface layer [,]                             |
| STATE_LAKE_SURF_LAYER_DZ         | double            | Vertical   thickness of surface water layer in lake [m]                                                          |
| STATE_LAKE_DEPTH                 | double            | distance   from surface to deepest point in lake [m]                                                             |
| STATE_LAKE_LAYER_SURF_AREA       | double            | Surface   area of liquid water in lake at each node [m2]                                                         |
| STATE_LAKE_SURF_AREA             | double            | Surface   area of liquid plus ice water on lake surface [m2]                                                     |
| STATE_LAKE_VOLUME                | double            | Lake   total volume including liquid water equivalent of lake ice [m3]                                           |
| STATE_LAKE_LAYER_TEMP            | double            | Lake   water temperature at each node [C]                                                                        |
| STATE_LAKE_AVERAGE_TEMP          | double            | Average   liquid water temperature of entire lake [C]                                                            |
| STATE_LAKE_ICE_AREA              | double            | Area of   ice coverage on lake at beginning of time step [m2]                                                    |
| STATE_LAKE_ICE_AREA_NEW          | double            | Area of   ice coverage on lake at end of time step [m2]                                                          |
| STATE_LAKE_ICE_WATER_EQUIVALENT  | double            | Liquid   water equivalent volume of lake ice [m3]                                                                |
| STATE_LAKE_ICE_HEIGHT            | double            | Lake ice   height at ghickest point [m]                                                                          |
| STATE_LAKE_ICE_TEMP              | double            | Lake ice   temperature [C]                                                                                       |
| STATE_LAKE_ICE_SWE               | double            | liquid   water equivalent depth of lake snow [m]                                                                 |
| STATE_LAKE_ICE_SNOW_SURF_TEMP    | double            | Temperature   of snow surface layer of lake snow [C]                                                             |
| STATE_LAKE_ICE_SNOW_PACK_TEMP    | double            | Temperature   of snow pack layer of lake snow [C]                                                                |
| STATE_LAKE_ICE_SNOW_COLD_CONTENT | double            | Snowpack   cold content of snow lake [J/m2]                                                                      |
| STATE_LAKE_ICE_SNOW_SURF_WATER   | double            | Liquid   water content of surface snow layer of lake snow [m]                                                    |
| STATE_LAKE_ICE_SNOW_PACK_WATER   | double            | Liquid   water content of pack snow layer of lake snow [m]                                                       |
| STATE_LAKE_ICE_SNOW_ALBEDO       | double            | Albedo of   lake snow [1]                                                                                        |
| STATE_LAKE_ICE_SNOW_DEPTH        | double            | Depth of   snow on lake ice [m]                                                                                  |
