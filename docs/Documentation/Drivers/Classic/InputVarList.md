# Input Variables

Below is a list of meteorological input variables. These input variables may always be found in the vicNl_def.h file.

NOTE: except where noted, variables are interpreted according to the following convention:

*   For moisture fluxes in units of mm, the variable is the total input over the time step.
*   For moisture fluxes in units of mm/s, the variable is the average input rate over the time step.
*   For **TMAX** and **TMIN**, the variable is the maximum or minimum over the time step.
*   For most other variables the variable is the average over the time step.

NOTE: For sub-daily time steps, VIC assigns the **beginning** time to the record. For example, if the [global parameter file](GlobalParam.md) variables FORCEYEAR, FORCEMONTH, FORCEDAY, and FORCEHOUR are 1948, 01, 01, and 00, respectively, and the variable FORCE_DT is set to 1 hour, then VIC assumes the first forcing record corresponds to the interval spanning 1948-01-01 00:00:00 and 1948-01-01 00:59:59\. VIC's output for the first time step will also be labelled as 1948 01 01 00\. The final time step of each day will be hour 23.

| Variable   	| Description                                             	| Units                 	|
|------------	|---------------------------------------------------------	|-----------------------	|
| AIR_TEMP   	| Average air temperature                                 	| C (ALMA_INPUT: K)     	|
| ALBEDO     	| Surface albedo                                          	| fraction              	|
| CHANNEL_IN 	| Incoming channel flow (total volume over the time step) 	| m<sup>3</sup>          	|
| CATM       	| Atmospheric CO2 mixing ratio                            	| ppm                   	|
| CRAINF     	| Convective rainfall                                     	| mm (ALMA_INPUT: mm/s) 	|
| CSNOWF     	| Convective snowfall                                     	| mm (ALMA_INPUT: mm/s) 	|
| DENSITY    	| Atmospheric density                                     	| kg/m<sup>3</sup>       	|
| FDIR       	| Fraction of incoming shortwave that is direct           	| fraction              	|
| LONGWAVE   	| Incoming longwave radiation                             	| W/m<sup>2</sup>        	|
| LSRAINF    	| Large-scale rainfall                                    	| mm (ALMA_INPUT: mm/s) 	|
| LSSNOWF    	| Large-scale snowfall                                    	| mm (ALMA_INPUT: mm/s) 	|
| PREC       	| Total precipitation (rain and snow)                     	| mm (ALMA_INPUT: mm/s) 	|
| PRESSURE   	| Atmospheric pressure                                    	| kPa (ALMA_INPUT: Pa)  	|
| QAIR       	| Specific humidity                                       	| kg/kg                 	|
| RAINF      	| Rainfall (convective and large-scale)                   	| mm (ALMA_INPUT: mm/s) 	|
| REL_HUMID  	| Relative humidity                                       	| fraction              	|
| SHORTWAVE  	| Incoming shortwave                                      	| W/m<sup>2</sup>        	|
| SNOWF      	| Snowfall (convective and large-scale)                   	| mm (ALMA_INPUT: mm/s) 	|
| TMAX       	| Maximum daily temperature                               	| C (ALMA_INPUT: K)     	|
| TMIN       	| Minimum daily temperature                               	| C (ALMA_INPUT: K)     	|
| TSKC       	| Cloud cover                                             	| fraction              	|
| VEGCOVER   	| partial veg cover fraction                              	| fraction              	|
| VP         	| Vapor pressure                                          	| kPa (ALMA_INPUT: Pa)  	|
| WIND       	| Wind speed                                              	| m/s                   	|
| WIND_E     	| Zonal component of wind speed                           	| m/s                   	|
| WIND_N     	| Meridional component of wind speed                      	| m/s                   	|
| SKIP       	| Place holder for unused data columns                    	| Non-Data              	|
