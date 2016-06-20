# Snow Output File - Image Driver

The snow file contains information about the snowpack, averaged across all elevation bands and vegetation tiles. The set of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the global parameter file. Snow files are always written, regardless of the mode of operation. Snow file name is `snow.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

The default snow output file is in netCDF format, with the same data structure described in the [image driver output format](OutputFormatting.md). the list of default flux variables are:

| Variable   Name                                                                                                | Dimension          | Units    | Description                                              |
|----------------------------------------------------------------------------------------------------------------|--------------------|----------|----------------------------------------------------------|
| OUT_SWE                                                                                                        | [time,   lat, lon] | mm       | Snow water   equivalent in snow pack                     |
| OUT_SNOW_DEPTH                                                                                                 | [time,   lat, lon] | cm       | Depth of   snow pack                                     |
| OUT_SNOW_CANOPY                                                                                                | [time,   lat, lon] | mm       | Snow   interception storage in canopy                    |
| OUT_SNOW_COVER                                                                                                 | [time,   lat, lon] | fraction | Fractional   area of snow cover                          |
| The following variables are only output when   FULL_ENERGY or FROZEN_SOIL is TRUE in the global parameter file |                    |          |                                                          |
| OUT_ADVECTION                                                                                                  | [time,   lat, lon] | W/m2     | Advected   energy                                        |
| OUT_DELTACC                                                                                                    | [time,   lat, lon] | W/m2     | Rate of   change in cold content in snow pack            |
| OUT_SNOW_FLUX                                                                                                  | [time,   lat, lon] | W/m2     | Energy   flux though snow pack                           |
| OUT_RFRZ_ENERGY                                                                                                | [time,   lat, lon] | W/m2     | Net energy   used to refreeze liquid water in snowpack   |
| OUT_MELT_ENERGY                                                                                                | [time,   lat, lon] | W/m2     | Energy of   fusion (melting) in snowpack                 |
| OUT_ADV_SENS                                                                                                   | [time,   lat, lon] | W/m2     | Net   sensible heat flux advected to snow pack           |
| OUT_LATENT_SUB                                                                                                 | [time,   lat, lon] | W/m2     | Net upward   latent heat flux due to sublimation         |
| OUT_SNOW_SURF_TEMP                                                                                             | [time,   lat, lon] | C        | Snow   surface temperature                               |
| OUT_SNOW_PACK_TEMP                                                                                             | [time,   lat, lon] | C        | Snow pack   temperature                                  |
| OUT_SNOW_MELT                                                                                                  | [time,   lat, lon] | mm       | Snow melt                                                |
| The following variables are only output when BLOWING   is TRUE in the global parameter file                    |                    |          |                                                          |
| OUT_SUB_BLOWING                                                                                                | [time,   lat, lon] | mm       | Net   sublimation of blowing snow                        |
| OUT_SUB_SURFACE                                                                                                | [time,   lat, lon] | mm       | Net   sublimation from snow pack surface                 |
| OUT_SUB_SNOW                                                                                                   | [time,   lat, lon] | mm       | Total   sublimation from snow pack (surface and blowing) |