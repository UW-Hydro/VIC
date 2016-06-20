# Default Flux Output Files - Image Driver

The primary output file type for the VIC Classic Driver is the flux file, which contains information about moisture and energy fluxes for each time step. The number of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the [global parameter file](GlobalParam.md); when either FULL_ENERGY or FROZEN_SOIL are true, the file will contain the same number of variables as found in the PILPS-2C files for comparison purposes. Flux files are always written, regardless of the mode of operation. Flux file name is `fluxes.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

The default flux output file is in netCDF format, with the same data structure described in the [image driver output format](OutputFormatting.md). the list of default flux variables are:

| Variable   Name | Dimension                  | Units    | Description                                                                                                                       |
|-----------------|----------------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------|
| OUT_PREC        | [time,   lat, lon]         | mm       | Precipitation   for current record                                                                                                |
| OUT_EVAP        | [time,   lat, lon]         | mm       | Evaporation   for current record                                                                                                  |
| OUT_RUNOFF      | [time,   lat, lon]         | mm       | Runoff for   current record                                                                                                       |
| OUT_BASEFLOW    | [time,   lat, lon]         | mm       | Baseflow   for current record                                                                                                     |
| OUT_WDEW        | [time,   lat, lon]         | mm       | Canopy   interception of liquid water                                                                                             |
| OUT_SOIL_LIQ    | [time,   nlayer, lat, lon] | mm       | Moisture   content of each soil layer                                                                                             |
| OUT_RAD_TEMP    | [time,   lat, lon]         | K        | Radiative   temperature of the surface (NOTE: only output when FULL_ENERGY or FROZEN_SOIL   is TRUE in the global parameter file) |
| OUT_SWNET       | [time,   lat, lon]         | W/m2     | Net   shortwave radiation at the surface                                                                                          |
| OUT_R_NET       | [time,   lat, lon]         | W/m2     | Net   radiation at the surface, includes long and shortwave radiation                                                             |
| OUT_LATENT      | [time,   lat, lon]         | W/m2     | Latent   heat from the surface (NOTE: only output when FULL_ENERGY or FROZEN_SOIL is   TRUE in the global parameter file)         |
| OUT_EVAP_CANOP  | [time,   lat, lon]         | mm       | Evaporation   from canopy storage                                                                                                 |
| OUT_TRANSP_VEG  | [time,   lat, lon]         | mm       | Transpiration   from the vegetation                                                                                               |
| OUT_EVAP_BARE   | [time,   lat, lon]         | mm       | Evaporation   from bare soil                                                                                                      |
| OUT_SUB_CANOP   | [time,   lat, lon]         | mm       | Sublimation   from canopy interception                                                                                            |
| OUT_SUB_SNOW    | [time,   lat, lon]         | mm       | Sublimation   from ground snow pack                                                                                               |
| OUT_SENSIBLE    | [time,   lat, lon]         | W/m2     | Sensible   heat flux from the surface (NOTE: only output when FULL_ENERGY or FROZEN_SOIL   is TRUE in the global parameter file)  |
| OUT_GRND_FLUX   | [time,   lat, lon]         | W/m2     | Ground   heat flux plus heat storage in the top soil layer                                                                        |
| OUT_DELTAH      | [time,   lat, lon]         | W/m2     | Rate of   change in heat storage                                                                                                  |
| OUT_FUSION      | [time,   lat, lon]         | W/m2     | Net energy   used to melt/freeze soil moisture                                                                                    |
| OUT_AERO_RESIST | [time,   lat, lon]         | s/m      | Aerodynamic   resistance                                                                                                          |
| OUT_SURF_TEMP   | [time,   lat, lon]         | C        | Surface   temperature                                                                                                             |
| OUT_ALBEDO      | [time,   lat, lon]         | fraction | Albedo of   surface cover                                                                                                         |
| OUT_REL_HUMID   | [time,   lat, lon]         | fraction | Relative   humidity                                                                                                               |
| OUT_IN_LONG     | [time,   lat, lon]         | W/m2     | Incoming   longwave at ground surface (under vegetation)                                                                          |
| OUT_AIR_TEMP    | [time,   lat, lon]         | C        | Air   temperature                                                                                                                 |
| OUT_WIND        | [time,   lat, lon]         | m/s      | Near   surface wind speed                                                                                                         |