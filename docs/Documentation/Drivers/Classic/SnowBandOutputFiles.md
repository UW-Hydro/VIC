# Snow Band Output Files

When the model is run with snow elevation bands, and default outputs are specified via [global parameter file](GlobalParam.md), snow pack information for each elevation band will be output to files in the results directory with the prefix `snow_band_`. Energy fluxes are output only for the full energy balance model, so there are file descriptions for both full energy and water balance model output files.

Each output file contains model output from one grid cell. The files use the naming convention `snow_band_xx.xx_yy.yy`, where `xx.xx` is the latitude, and `yy.yy` is the longitude. The number of decimal places in the output filename is determined by GRID_DECIMAL in the [global parameter file](GlobalParam.md), while the latitude and longitude values come from the [soil parameter file](SoilParam.md).

## Snow Band File Format

| Column                                          	| Variable Name        	| Units              	| Description                                                                     	|
|-------------------------------------------------	|----------------------	|--------------------	|---------------------------------------------------------------------------------	|
| 1                                               	| year                 	| year               	| Year of current record                                                          	|
| 2                                               	| month                	| month              	| Month of current record                                                         	|
| 3                                               	| day                  	| day                	| Day of current record                                                           	|
| 4                                               	| hour                 	| hour               	| Hour of current record (present only if model is run at a sub-daily time step). 	|
| The following columns repeat for each snow band 	|                      	|                    	|                                                                                 	|
| 5+(SNOW_BAND-1)*7                               	| OUT_SWE_BAND         	| mm                 	| Snow water equivalence (the amount of water stored in the snow pack)            	|
| 6+(SNOW_BAND-1)*7                               	| OUT_SNOW_DEPTH_BAND  	| cm                 	| Snow pack depth                                                                 	|
| 7+(SNOW_BAND-1)*7                               	| OUT_SNOW_CANOPY_BAND 	| mm                 	| Amount of snow stored in the canopy is present                                  	|
| 8+(SNOW_BAND-1)*7                               	| OUT_SNOW_COVER_BAND  	| fraction           	| Fractional area of snow cover                                                   	|
| 9+(SNOW_BAND-1)*7                               	| OUT_ADVECTION_BAND   	| W/m2               	| Advective flux into the canopy by rain                                          	|
| 10+(SNOW_BAND-1)*7                              	| OUT_DELTACC_BAND     	| W/m2               	| Change in the cold content or energy storage of the snow pack                   	|
| 11+(SNOW_BAND-1)*7                              	| OUT_SNOW_FLUX_BAND   	| W/m2               	| Thermal flux through the snow pack                                              	|
| 12+(SNOW_BAND-1)*7                              	| OUT_RFRZ_ENERGY_BAND 	| W/m2               	| Energy used to refreeze the snow pack                                           	|
| 13+(SNOW_BAND-1)*7                              	| OUT_MELT_ENERGY_BAND 	| W/m2               	| Energy of fusion (melting) in snow pack                                         	|
| 14+(SNOW_BAND-1)*7                              	| OUT_ADV_SENS_BAND    	| W/m2               	| Net sensible heat flux advected to snow pack                                    	|
| 15+(SNOW_BAND-1)*7                              	| OUT_LATENT_SUB_BAND  	| W/m2               	| Net upward latent heat flux due to sublimation                                  	|
| 16+(SNOW_BAND-1)*7                              	| OUT_SNOW_SURFT_BAND  	| C                 	| Snow surface temperature                                                        	|
| 17+(SNOW_BAND-1)*7                              	| OUT_SNOW_PACKT_BAND  	| C                 	| Snow pack temperature                                                           	|
| 18+(SNOW_BAND-1)*7                              	| OUT_MELT_BAND        	| mm                 	| Amount of snow melt                                                             	|
| 19+(SNOW_BAND-1)*7                              	| OUT_NET_SHORT_BAND   	| W/m2               	| Net downward shortwave flux                                                     	|
| 20+(SNOW_BAND-1)*7                              	| OUT_NET_LONG_BAND    	| W/m2               	| Net downward longwave flux                                                      	|
| 21+(SNOW_BAND-1)*7                              	| OUT_ALBEDO_BAND      	| fraction           	| Average surface albedo                                                          	|
| 22+(SNOW_BAND-1)*7                              	| OUT_LATENT_BAND      	| W/m2               	| Net upward latent heat flux                                                     	|
| 23+(SNOW_BAND-1)*7                              	| OUT_SENSIBLE_BAND    	| W/m2               	| Net upward sensible heat flux                                                   	|
| 24+(SNOW_BAND-1)*7                              	| OUT_GRND_FLUX_BAND   	| W/m2               	| Net heat flux into ground                                                       	|
