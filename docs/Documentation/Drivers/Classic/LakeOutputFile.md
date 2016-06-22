# Lake Output File

The lake file contains information about the lake fraction of the grid cell. Lake files are only written when LAKES is set equal to a valid lake parameter file, in the [global parameter file](GlobalParam.md). Lake files begin with the prefix `lake_`.

Each output file contains model output from one grid cell. The files use the naming convention `fluxes_xx.xx_yy.yy`, where `xx.xx` is the latitude, and `yy.yy` is the longitude. The number of decimal places in the output filename is determined by the `GRID_DECIMAL` in the [global parameter file](GlobalParam.md), while the latitude and longitude values come from the [soil parameter file](SoilParam.md).

| Column 	| Variable Name       	| Units        	| Description                                             	|
|--------	|---------------------	|--------------	|---------------------------------------------------------	|
| 1      	| OUT_LAKE_ICE_TEMP   	| C             | Temperature of lake ice                                 	|
| 2      	| OUT_LAKE_ICE_HEIGHT 	| cm            | Thickness of lake ice                                   	|
| 3      	| OUT_LAKE_ICE_FRACT  	| fraction      | Fractional coverage of lake ice                         	|
| 4      	| OUT_LAKE_DEPTH      	| m             | Lake depth (distance between surface and deepest point) 	|
| 5      	| OUT_LAKE_SURF_AREA  	| m<sup>2</sup>	| Lake surface area                                       	|
| 6      	| OUT_LAKE_VOLUME     	| m<sup>3</sup>	| Lake volume                                             	|
| 7      	| OUT_LAKE_SURF_TEMP  	| C             | Lake surface temperature                                	|
| 8      	| OUT_EVAP_LAKE       	| mm            | Net evaporation from lake surface                       	|
