# VIC Forcings Files

The VIC Classic Driver requires subdaily forcings (meteorological or other).  The required forcing variables vary depending options set in the global parameter file.

#### Meteorological Forcings, Required in all simulations:

| Variable   | Description                         | Units           |
|------------|-------------------------------------|---------------- |
| AIR_TEMP   | Average air temperature             | C               |
| PREC       | Total precipitation (rain and snow) | mm              |
| PRESSURE   | Atmospheric pressure                | kPa             |
| SWDOWN     | Incoming shortwave radiation        | W/m<sup>2</sup> |
| LWDOWN     | Incoming longwave radiation         | W/m<sup>2</sup> |
| VP         | Vapor pressure                      | kPa             |
| WIND       | Wind speed                          | m/s             |

#### Vegetation Timeseries Forcings (Optional):

| Variable   | Description                                              | Units                       |
|------------|----------------------------------------------------------|---------------------------- |
| ALBEDO     | Surface Albedo                                           | fraction (between 0 and 1)  |
| LAI_IN     | Leaf Area Index                                          | m<sup>2</sup>/m<sup>2</sup> |
| VEGCOVER   | Partial veg cover fraction ( = 1 - canopy gap fraction ) | fraction (between 0 and 1)  |

#### Lake Forcings, Required when LAKES is TRUE:

| Variable   | Description                                              | Units           |
|------------|----------------------------------------------------------|---------------- |
| CHANNEL_IN | Incoming channel flow (total volume over the time step)  | m<sup>3</sup>   |

#### Carbon Cycle Forcings, Required when CARBON is TRUE:

| Variable   | Description                                   | Units           |
|------------|-----------------------------------------------|---------------- |
| CATM       | Atmospheric CO2 mixing ratio                  | ppm             |
| FDIR       | Fraction of incoming shortwave that is direct | fraction        |
| PAR        | Photosynthetically active radiation           | W/m<sup>2</sup> |

As of October 2015, work is currently underway to develop a next generation meteorological forcing generator. This work, in part, will support the development of VIC forcings.  Follow or contribute to the development of these tools by visiting https://github.com/jhamman/mtclim5.
