# Frozen Soil Output File

When the model is run with the frozen soil algorithms, a third output file is produced which contains soil thermal output parameters. Since the frozen soils algorithm only works in with the full energy balance, there is only one format for the output file.

Each output file contains model output from one grid cell. The files use the naming convention `fdepth_xx.xx_yy.yy`, where `xx.xx` is the latitude, and `yy.yy` is the longitude. The number of decimal places in the output filename is determined by GRID_DECIMAL in the [global parameter file](GlobalParam.md), while the latitude and longitude values come from the [soil parameter file](SoilParam.md).

* * *

## Frozen Soil Output File

| Column                              	| Variable Name       	| Units    	| Description                                                                     	|
|-------------------------------------	|---------------------	|----------	|---------------------------------------------------------------------------------	|
| 1                                   	| year                	| year     	| Year of current record                                                          	|
| 2                                   	| month               	| month    	| Month of current record                                                         	|
| 3                                   	| day                 	| day      	| Day of current record                                                           	|
| 4                                   	| hour                	| hour     	| Hour of current record (present only if model is run at a sub-daily time step). 	|
| Repeats for all frost fronts        	|                     	|          	|                                                                                 	|
| 5:(Nfronts\*2)+4:2                   	| OUT_FDEPTH          	| cm       	| Depth of the freezing front                                                     	|
| 6:(Nfronts\*2)+5:2                   	| OUT_TDEPTH          	| cm       	| Depth of the thawing front                                                      	|
| End of file does not repeat.        	|                     	|          	|                                                                                 	|
| (Nfronts\*2)+6:(Nfronts\*2)+Nlayers+5 | OUT_SOIL_MOIST      	| mm       	| Total soil moisture in each layer (liquid water plus ice)                       	|
| (Nfronts\*2)+7:(Nfronts\*2)+Nlayers+5 | OUT_SOIL_MOIST_FRAC 	| fraction 	| Fraction of total soil moisture in each layer                                   	|
