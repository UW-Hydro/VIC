# VIC Elevation Band Parameter File

By default, VIC assumes each grid cell is flat. This assumption can lead to inaccuracies in snow pack estimates in mountainous grid cells. For this reason, the option exists to have VIC break each grid cell up into a number of _elevation bands_ (also called _snow bands_) and simulate them separately. Each band's mean elevation is used to lapse the grid cell average temperature, pressure, and precipitation to a more accurate local estimate.

This file contains information needed to define the properties of each elevation band used by the snow model. Snow elevation bands are used to improve the model's performance in areas with pronounced topography, especially mountainous regions, where the effects of elevation on snow pack accumulation and ablation might be lost in a large grid cell.

The number of snow elevation bands (_option.SNOW_BAND_) to be used with the model is defined in the [global parameter file](GlobalParam.md). The elevation band (or snow band) file is only read if the number of snow elevation bands is greater than 1.

It is not necessary that all grid cells in a basin have the same number of elevation bands. _SNOW_BAND_ is simply the maximum number of elevation bands anywhere in the basin. For relatively flat grid cells, some of the elevation bands will have _AreaFract_ values of 0\. For these zero-area bands, a value of 0 may be supplied for _elevation_ and _Pfactor_.

## Elevation Band File Format

| Column                                | Variable Name     | Units     | Number of Values  | Description                                                                                                                                                                               |
|-----------------------------------    |---------------    |---------- |------------------ |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  |
| 1                                     | cellnum           | N/A       | 1                 | Grid cell number (should match numbers assigned in soil parameter file)                                                                                                                   |
| 2 : (SNOW_BAND+1)                     | AreaFract         | fraction  | SNOW_BAND         | Fraction of grid cell covered by each elevation band. Sum of the fractions must equal 1.                                                                                                  |
| (SNOW_BAND+2) : (2\*SNOW_BAND+1)      | elevation         | m         | SNOW_BAND         | Mean (or median) elevation of elevation band. This is used to compute the change in air temperature from the grid cell mean elevation.                                                    |
| (2\*SNOW_BAND+2) : (3\*SNOW_BAND+1)   | Pfactor           | fraction  | SNOW_BAND         | Fraction of cell precipitation that falls on each elevation band. Total must equal 1. To ignore effects of elevation on precipitation, set these fractions equal to the area fractions.   |
