# Snow Output File

The snow file contains information about the snowpack, averaged across all elevation bands and vegetation tiles. The set of variables in this file depends on the values of `FULL_ENERGY` and `FROZEN_SOIL` in the global parameter file. Snow files are always written, regardless of the mode of operation. Snow files begin with the prefix `snow_`.

Each output file contains model output from one grid cell. The files use the naming convention `fluxes_xx.xx_yy.yy`, where `xx.xx` is the latitude, and `yy.yy` is the longitude. The number of decimal places in the output filename is determined by `GRID_DECIMAL` in the [global parameter file](GlobalParam.md), while the latitude and longitude values come from the [soil parameter file](SoilParam.md).

| Column                                                                                                       	| Variable Name      	| Units             | Description                                            	|
|--------------------------------------------------------------------------------------------------------------	|--------------------	|------------------ |--------------------------------------------------------	|
| 1                                                                                                            	| OUT_SWE            	| mm                | Snow water equivalent in snow pack                     	|
| 2                                                                                                            	| OUT_SNOW_DEPTH     	| cm               	| Depth of snow pack                                     	|
| 3                                                                                                            	| OUT_SNOW_CANOPY    	| mm               	| Snow interception storage in canopy                    	|
| 4                                                                                                            	| OUT_SNOW_COVER     	| fraction         	| Fractional area of snow cover                          	|
| The following variables are only output when FULL_ENERGY or FROZEN_SOIL is TRUE in the global parameter file 	|                    	|                  	|                                                        	|
| 5                                                                                                            	| OUT_ADVECTION      	| W/m<sup>2</sup>   | Advected energy                                        	|
| 6                                                                                                            	| OUT_DELTACC        	| W/m<sup>2</sup>  	| Rate of change in cold content in snow pack            	|
| 7                                                                                                            	| OUT_SNOW_FLUX      	| W/m<sup>2</sup>   | Energy flux though snow pack                           	|
| 8                                                                                                            	| OUT_RFRZ_ENERGY    	| W/m<sup>2</sup>   | Net energy used to refreeze liquid water in snowpack   	|
| 9                                                                                                            	| OUT_MELT_ENERGY    	| W/m<sup>2</sup>   | Energy of fusion (melting) in snowpack                 	|
| 10                                                                                                           	| OUT_ADV_SENS       	| W/m<sup>2</sup>   | Net sensible heat flux advected to snow pack           	|
| 11                                                                                                           	| OUT_LATENT_SUB     	| W/m<sup>2</sup>   | Net upward latent heat flux due to sublimation         	|
| 12                                                                                                           	| OUT_SNOW_SURF_TEMP 	| C               	| Snow surface temperature                               	|
| 13                                                                                                           	| OUT_SNOW_PACK_TEMP 	| C               	| Snow pack temperature                                  	|
| 14                                                                                                           	| OUT_SNOW_MELT      	| mm              	| Snow melt                                              	|
| The following variables are only output when BLOWING is TRUE in the global parameter file                    	|                    	|                   |                                                        	|
| 15                                                                                                           	| OUT_SUB_BLOWING    	| mm              	| Net sublimation of blowing snow                        	|
| 16                                                                                                           	| OUT_SUB_SURFACE    	| mm              	| Net sublimation from snow pack surface                 	|
| 17                                                                                                           	| OUT_SUB_SNOW       	| mm              	| Total sublimation from snow pack (surface and blowing) 	|
