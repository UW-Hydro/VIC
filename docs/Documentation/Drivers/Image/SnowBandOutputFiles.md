# Snow Band Output Files

When the model is run with snow elevation bands, and PRT_SNOW_BAND is turned on in the [global parameter file](GlobalParam.md), snow pack information for each elevation band will be output to files in the results directory with file name `snow.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md). Energy fluxes are output only for the full energy balance model, so there are file descriptions for both full energy and water balance model output files.

The default flux output file is in netCDF format, with the same data structure described in the [image driver output format](OutputFormatting.md). the list of default flux variables are:

| Variable   Name      | Dimension                   | Units    | Description                                                                                                                              |
|----------------------|-----------------------------|----------|------------------------------------------------------------------------------------------------------------------------------------------|
| OUT_SWE_BAND         | [time, snow_band, lat, lon] | mm       | Snow   water equivalence (the amount of water stored in the snow pack)                                                                   |
| OUT_SNOW_DEPTH_BAND  | [time, snow_band, lat, lon] | cm       | Snow pack   depth                                                                                                                        |
| OUT_SNOW_CANOPY_BAND | [time, snow_band, lat, lon] | mm       | Amount of   snow stored in the canopy is present                                                                                         |
| OUT_ADVECTION_BAND   | [time, snow_band, lat, lon] | W/m2     | Advective   flux into the canopy by rain (NOTE: only output when FULL_ENERGY=TRUE in the   global parameter file)                        |
| OUT_DELTACC_BAND     | [time, snow_band, lat, lon] | W/m2     | Change in   the cold content or energy storage of the snow pack (NOTE: only output when   FULL_ENERGY=TRUE in the global parameter file) |
| OUT_SNOW_FLUX_BAND   | [time, snow_band, lat, lon] | W/m2     | Thermal   flux through the snow pack (NOTE: only output when FULL_ENERGY=TRUE in the   global parameter file)                            |
| OUT_RFRZ_ENERGY_BAND | [time, snow_band, lat, lon] | W/m2     | Energy   used to refreeze the snow pack (NOTE: only output when FULL_ENERGY=TRUE in   the global parameter file)                         |
| OUT_SWNET_BAND       | [time, snow_band, lat, lon] | W/m2     | Net   downward shortwave flux                                                                                                            |
| OUT_LWNET_BAND       | [time, snow_band, lat, lon] | W/m2     | Net   downward longwave flux                                                                                                             |
| OUT_ALBEDO_BAND      | [time, snow_band, lat, lon] | fraction | Average   surface albedo                                                                                                                 |
| OUT_LATENT_BAND      | [time, snow_band, lat, lon] | W/m2     | Net   upward latent heat flux                                                                                                            |
| OUT_SENSIBLE_BAND    | [time, snow_band, lat, lon] | W/m2     | Net   upward sensible heat flux                                                                                                          |
| OUT_GRND_FLUX_BAND   | [time, snow_band, lat, lon] | W/m2     | Net heat   flux into ground                                                                                                              |
