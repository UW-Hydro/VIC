# Snow Band Output Files

When the model is run with snow elevation bands, and PRT_SNOW_BAND is turned on in the [global parameter file](GlobalParam.md), snow pack information for each elevation band will be output to files in the results directory with the prefix "snow_band_". Energy fluxes are output only for the full energy balance model, so there are file descriptions for both full energy and water balance model output files.

Each output file contains model output from one grid cell. The files use the naming convention "snow_band_xx.xx_yy.yy", where xx.xx is the latitude, and yy.yy is the longitude. The number of decimal places in the output filename is determined by GRID_DECIMAL in the [global parameter file](GlobalParam.md), while the latitude and longitude values come from the [soil parameter file](SoilParam.md).

| Column                                                                                             | Variable Name        | Units    | Description                                                            |
|----------------------------------------------------------------------------------------------------|----------------------|----------|------------------------------------------------------------------------|
| All the following columns repeat for each snow band                                                |                      |          |                                                                        |
| 5+(SNOW_BAND-1)*7                                                                                  | OUT_SWE_BAND         | mm       | Snow   water equivalence (the amount of water stored in the snow pack) |
| 6+(SNOW_BAND-1)*7                                                                                  | OUT_SNOW_DEPTH_BAND  | cm       | Snow pack   depth                                                      |
| 7+(SNOW_BAND-1)*7                                                                                  | OUT_SNOW_CANOPY_BAND | mm       | Amount of   snow stored in the canopy is present                       |
| (The   following 4 variables are only output when FULL_ENERGY=TRUE in the global   parameter file) |                      |          |                                                                        |
| 9+(SNOW_BAND-1)*7                                                                                  | OUT_ADVECTION_BAND   | W/m2     | Advective   flux into the canopy by rain                               |
| 10+(SNOW_BAND-1)*7                                                                                 | OUT_DELTACC_BAND     | W/m2     | Change in   the cold content or energy storage of the snow pack        |
| 11+(SNOW_BAND-1)*7                                                                                 | OUT_SNOW_FLUX_BAND   | W/m2     | Thermal   flux through the snow pack                                   |
| 12+(SNOW_BAND-1)*7                                                                                 | OUT_RFRZ_ENERGY_BAND | W/m2     | Energy   used to refreeze the snow pack                                |
|                                                                                                    |                      |          |                                                                        |
| 19+(SNOW_BAND-1)*7                                                                                 | OUT_SWNET_BAND       | W/m2     | Net   downward shortwave flux                                          |
| 20+(SNOW_BAND-1)*7                                                                                 | OUT_LWNET_BAND       | W/m2     | Net   downward longwave flux                                           |
| 21+(SNOW_BAND-1)*7                                                                                 | OUT_ALBEDO_BAND      | fraction | Average   surface albedo                                               |
| 22+(SNOW_BAND-1)*7                                                                                 | OUT_LATENT_BAND      | W/m2     | Net   upward latent heat flux                                          |
| 23+(SNOW_BAND-1)*7                                                                                 | OUT_SENSIBLE_BAND    | W/m2     | Net   upward sensible heat flux                                        |
| 24+(SNOW_BAND-1)*7                                                                                 | OUT_GRND_FLUX_BAND   | W/m2     | Net heat   flux into ground                                            |
