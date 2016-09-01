# Output Variables

Using options within the *global parameter file*, any combination of the variables listed below may be output by VIC.

## Water Balance Terms (state variables)
| Variable            	| Description                                                                                 | Units                 |
|---------------------	|-------------------------------------------------------------------------------------------- |---------------------- |
| OUT_ASAT            	| Saturated Area Fraction (of exposed land, i.e. non-lake; the total fraction of the grid cell that is completely full of water would be OUT_ASAT plus OUT_LAKE_AREA_FRAC) 	| fraction |
| OUT_LAKE_AREA_FRACT 	| Lake surface area as fraction of grid cell area                                             | fraction              |
| OUT_LAKE_DEPTH      	| Lake depth (distance between surface and deepest point)                                     | m                     |
| OUT_LAKE_ICE        	| Moisture stored as lake ice                                                                 | mm over lake ice area |
| OUT_LAKE_ICE_FRACT  	| Fractional coverage of lake ice                                                             | fraction              |
| OUT_LAKE_ICE_HEIGHT 	| Thickness of lake ice                                                                       | cm  |
| OUT_LAKE_MOIST      	| Liquid water and ice stored in lake                                                         | mm over grid cell     |
| OUT_LAKE_SURF_AREA  	| Lake surface area                                                                           | m<sup>2</sup>         |
| OUT_LAKE_SWE        	| Liquid water equivalent of snow on top of lake ice                                          | m over lake ice area  |
| OUT_LAKE_SWE_V      	| Volumetric liquid water equivalent of snow on top of lake ice                               | m<sup>3</sup>         |
| OUT_LAKE_VOLUME     	| Lake volume                                                                                 | m<sup>3</sup>         |
| OUT_ROOTMOIST       	| Total soil moisture in layers that contain roots.                                           | mm                    |
| OUT_SMFROZFRAC      	| Fraction of soil moisture (by mass) that is ice, for each soil layer                        | fraction              |
| OUT_SMLIQFRAC       	| Fraction of soil moisture (by mass) that is liquid, for each soil layer                     | fraction              |
| OUT_SNOW_CANOPY     	| Snow interception storage in canopy                                                         | mm                    |
| OUT_SNOW_COVER      	| Fractional area of snow cover                                                               | fraction              |
| OUT_SNOW_DEPTH      	| Depth of snow pack                                                                          | cm                    |
| OUT_SOIL_ICE        	| Soil ice content for each soil layer                                                        | mm                    |
| OUT_SOIL_LIQ        	| Soil liquid content for each soil layer                                                     | mm                    |
| OUT_SOIL_ICE_FRAC   	| Fractional soil ice content for each soil layer                                             | fraction              |
| OUT_SOIL_LIQ_FRAC   	| Fractional soil liquid content for each soil layer                                          | fraction              |
| OUT_SOIL_MOIST      	| Total soil moisture content for each soil layer                                             | mm                    |
| OUT_SOIL_WET        	| Vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) | mm/mm                 |
| OUT_SURFSTOR        	| Storage of liquid water and ice (not snow) on surface (ponding)                             | mm                    |
| OUT_SURF_FROST_FRAC 	| Fraction of soil surface that is frozen                                                     | fraction              |
| OUT_SWE             	| Snow water equivalent in snow pack (including vegetation-intercepted snow)                  | mm                    |
| OUT_WDEW            	| Total moisture interception storage in canopy                                               | mm                    |
| OUT_ZWT             	| Water table position, using lowest unsaturated soil layer                                   | cm (positive upwards, i.e. negative values indicate below soil surface; 0 = at soil surface) 	|
| OUT_ZWT_LUMPED      	| Water table position, lumping all layers' moistures together                                | cm (positive upwards, i.e. negative values indicate below soil surface; 0 = at soil surface) 	|

## Water Balance Terms (fluxes)
| Variable               	| Description                                                                               	| Units        	|
|------------------------	|-------------------------------------------------------------------------------------------	|-------------	|
| OUT_BASEFLOW           	| Baseflow out of the bottom layer                                                          	| mm  	        |
| OUT_DELINTERCEPT       	| Change in canopy interception storage                                                     	| mm            |
| OUT_DELSOILMOIST       	| Change in soil water content                                                              	| mm            |
| OUT_DELSURFSTOR        	| Change in surface liquid water storage                                                    	| mm            |
| OUT_DELSWE             	| Change in snow water equivalent                                                           	| mm            |
| OUT_EVAP               	| Total net evaporation                                                                     	| mm  	        |
| OUT_EVAP_BARE          	| Net evaporation from bare soil                                                            	| mm  	        |
| OUT_EVAP_CANOP         	| Net evaporation from canopy interception                                                  	| mm  	        |
| OUT_INFLOW             	| Moisture that reaches top of soil column                                                  	| mm  	        |
| OUT_LAKE_BF_IN         	| Incoming baseflow from lake catchment                            	                          | mm  	        |
| OUT_LAKE_BF_IN_V       	| Incoming volumetric baseflow from lake catchment                 	                          | m<sup>3</sup> |
| OUT_LAKE_BF_OUT        	| Outgoing baseflow from lake to channel network                  	                          | mm  	        |
| OUT_LAKE_BF_OUT_V      	| Outgoing volumetric baseflow from lake to channel network        	                          | m<sup>3</sup> |
| OUT_LAKE_CHANNEL_IN    	| Channel inflow from upstream                                    	                          | mm  	        |
| OUT_LAKE_CHANNEL_IN_V  	| Volumetric channel inflow from upstream                          	                          | m<sup>3</sup> |
| OUT_LAKE_CHANNEL_OUT   	| Channel outflow from lake to channel network                     	                          | mm  	        |
| OUT_LAKE_CHANNEL_OUT_V 	| Volumetric channel outflow from lake to channel network          	                          | m<sup>3</sup> |
| OUT_LAKE_DSTOR         	| Change in lake moisture storage (liquid plus ice cover)         	                          | mm  	        |
| OUT_LAKE_DSTOR_V       	| Volumetric change in lake moisture storage (liquid plus ice cover)                          | m<sup>3</sup> |
| OUT_LAKE_DSWE          	| Change in swe on top of lake ice                                	                          | mm  	        |
| OUT_LAKE_DSWE_V        	| Volumetric change in swe on top of lake ice                      	                          | m<sup>3</sup> |
| OUT_LAKE_EVAP          	| Net evaporation from lake surface                                	                          | mm  	        |
| OUT_LAKE_EVAP_V        	| Net volumetric evaporation from lake surface                     	                          | m<sup>3</sup> |
| OUT_LAKE_PREC_V        	| Volumetric precipitation over lake surface                       	                          | m<sup>3</sup> |
| OUT_LAKE_RCHRG         	| Recharge from lake to surrounding wetland                        	                          | mm  	        |
| OUT_LAKE_RCHRG_V       	| Volumetric recharge from lake to surrounding wetland             	                          | m<sup>3</sup> |
| OUT_LAKE_RO_IN         	| Incoming runoff from lake catchment                             	                          | mm  	        |
| OUT_LAKE_RO_IN_V       	| Incoming volumetric runoff from lake catchment                   	                          | m<sup>3</sup> |
| OUT_LAKE_VAPFLX        	| Outgoing sublimation from snow on top of lake ice                	                          | mm  	        |
| OUT_LAKE_VAPFLX_V      	| Outgoing volumetric sublimation from snow on top of lake ice     	                          | m<sup>3</sup> |
| OUT_PET                	| Potential evapotranspiration (= area-weighted sum of potential transpiration and potential soil evaporation).  Potential transpiration is computed using the Penman-Monteith eqn with architectural resistance and LAI of the current veg cover. | mm  	        |
| OUT_PREC               	| Incoming precipitation                                                                    	| mm  	        |
| OUT_RAINF              	| Rainfall                                                                                  	| mm  	        |
| OUT_REFREEZE           	| Refreezing of water in the snow                                                           	| mm  	        |
| OUT_RUNOFF             	| Surface runoff                                                                            	| mm  	        |
| OUT_SNOW_MELT          	| Snow melt                                                                                 	| mm  	        |
| OUT_SNOWF              	| Snowfall                                                                                  	| mm  	        |
| OUT_SUB_BLOWING        	| Net sublimation of blowing snow                                                           	| mm  	        |
| OUT_SUB_CANOP          	| Net sublimation from snow stored in canopy                                                	| mm  	        |
| OUT_SUB_SNOW           	| Total net sublimation from snow pack (surface and blowing)                                	| mm  	        |
| OUT_SUB_SURFACE        	| Net sublimation from snow pack surface                                                    	| mm  	        |
| OUT_TRANSP_VEG         	| Net transpiration from vegetation                                                         	| mm  	        |
| OUT_WATER_ERROR        	| Water budget error                                                                        	| mm          	|

## Energy Balance Terms (state variables)
| Variable           	| Description                                                	| Units    |
|--------------------	|------------------------------------------------------------	|--------- |
| OUT_ALBEDO         	| Average surface albedo                                     	| fraction |
| OUT_BARESOILT      	| Bare soil surface temperature                              	| C        |
| OUT_FDEPTH         	| Depth of freezing fronts for each freezing front           	| cm       |
| OUT_LAKE_ICE_TEMP  	| Temperature of lake ice                                    	| K        |
| OUT_LAKE_SURF_TEMP 	| Lake surface temperature                                   	| K        |
| OUT_RAD_TEMP       	| Average radiative surface temperature                      	| K        |
| OUT_SALBEDO        	| Snow pack albedo                                           	| fraction |
| OUT_SNOW_PACK_TEMP 	| Snow pack temperature                                      	| C        |
| OUT_SNOW_SURF_TEMP 	| Snow surface temperature                                   	| C        |
| OUT_SNOWT_FBFLAG   	| Snow surface temperature fallback flag                     	| 0 or 1   |
| OUT_SOIL_TEMP      	| Soil temperature for each soil layer                       	| C        |
| OUT_SOIL_TNODE     	| Soil temperature for each soil thermal node                	| C        |
| OUT_SOIL_TNODE_WL  	| Soil temperature for each soil thermal node in the wetland 	| C        |
| OUT_SOILT_FBFLAG   	| Soil temperature flag for each soil thermal node           	| 0 or 1   |
| OUT_SURF_TEMP      	| Average surface temperature                                	| C        |
| OUT_SURFT_FBFLAG   	| Surface temperature fallback flag                          	| 0 or 1   |
| OUT_TCAN_FBFLAG    	| Tcanopy fallback flag                                      	| 0 or 1   |
| OUT_TDEPTH         	| Depth of thawing fronts for each thawing front             	| cm       |
| OUT_TFOL_FBFLAG    	| Tfoliage fallback flag                                     	| 0 or 1   |
| OUT_VEGT           	| Average vegetation canopy temperature                      	| C        |

## Energy Balance Terms (fluxes)
| Variable         	| Description                                          	| Units           |
|------------------	|------------------------------------------------------	|---------------- |
| OUT_ADV_SENS     	| Net sensible flux advected to snow pack              	| W/m<sup>2</sup> |
| OUT_ADVECTION    	| Advected energy                                      	| W/m<sup>2</sup> |
| OUT_DELTACC      	| Rate of change in cold content in snow pack          	| W/m<sup>2</sup> |
| OUT_DELTAH       	| Rate of change in heat storage                       	| W/m<sup>2</sup> |
| OUT_ENERGY_ERROR 	| Energy budget error                                  	| W/m<sup>2</sup> |
| OUT_FUSION       	| Net energy used to melt/freeze soil moisture         	| W/m<sup>2</sup> |
| OUT_GRND_FLUX    	| Net heat flux into ground                            	| W/m<sup>2</sup> |
| OUT_IN_LONG      	| Incoming longwave at ground surface (under veg)      	| W/m<sup>2</sup> |
| OUT_LATENT       	| Net upward latent heat flux                          	| W/m<sup>2</sup> |
| OUT_LATENT_SUB   	| Net upward latent heat flux from sublimation         	| W/m<sup>2</sup> |
| OUT_MELT_ENERGY  	| Energy of fusion (melting) in snowpack               	| W/m<sup>2</sup> |
| OUT_LWNET        	| Net downward longwave flux                           	| W/m<sup>2</sup> |
| OUT_SWNET        	| Net downward shortwave flux                          	| W/m<sup>2</sup> |
| OUT_R_NET        	| Net downward radiation flux                          	| W/m<sup>2</sup> |
| OUT_RFRZ_ENERGY  	| Net energy used to refreeze liquid water in snowpack 	| W/m<sup>2</sup> |
| OUT_SENSIBLE     	| Net upward sensible heat flux                        	| W/m<sup>2</sup> |
| OUT_SNOW_FLUX    	| Energy flux through snow pack                        	| W/m<sup>2</sup> |

## Miscellaneous Terms
| OUT_AERO_COND    	| Scene aerodynamic conductance (tiles with overstory contribute overstory conductance; others contribute surface conductance)     	| m/s              |
|------------------	|----------------------------------------------------------------------------------------------------------------------------------	|----------------- |
| OUT_AERO_COND   	| "scene" aerodynamic conductance (tiles with overstory contribute overstory conductance; others contribute surface conductance)  	| m/s              |
| OUT_AERO_COND1   	| Surface aerodynamic conductance                                                                                                  	| m/s              |
| OUT_AERO_COND2   	| Overstory aerodynamic conductance                                                                                                	| m/s              |
| OUT_AERO_RESIST  	| Scenecanopy aerodynamic resistance (tiles with overstory contribute over story resistance; others contribute surface resistance) 	| s/m              |
| OUT_AERO_RESIST1 	| Surface aerodynamic resistance                                                                                                   	| s/m              |
| OUT_AERO_RESIST2 	| Overstory aerodynamic resistance                                                                                                 	| s/m              |
| OUT_AIR_TEMP     	| Air temperature                                                                                                                  	| C     	         |
| OUT_CATM       	  | atmospheric CO2 concentrtaion                                                                                                    	| ppm    	         |
| OUT_DENSITY      	| Near-surface atmospheric density                                                                                                 	| kg/m<sup>3</sup> |
| OUT_FCANOPY     	| Vegetation canopy cover fraction                                                                                                	| fraction         |
| OUT_FDIR         	| fraction of incoming shortwave that is direct                                                                                    	| fraction         |
| OUT_LAI          	| Leaf Area Index                                                                                                                  	| fraction         |
| OUT_LWDOWN       	| Incoming longwave                                                                                                              	  | W/m<sup>2</sup>  |
| OUT_PAR          	| Incoming photosynthetically active radiation                                                                                     	| W/m<sup>2</sup>  |
| OUT_PRESSURE     	| Near surface atmospheric pressure                                                                                                	| kPa        	     |
| OUT_QAIR         	| Specific humidity                                                                                                                	| kg/kg            |
| OUT_REL_HUMID    	| Relative humidity                                                                                                                	| fraction         |
| OUT_SWDOWN      	| Incoming shortwave                                                                                                               	| W/m<sup>2</sup>  |
| OUT_SURF_COND    	| Surface conductance                                                                                                              	| m/s              |
| OUT_VP           	| Near surface vapor pressure                                                                                                      	| kPa  	           |
| OUT_VPD          	| Near surface vapor pressure deficit                                                                                              	| kPa  	           |
| OUT_WIND         	| Near surface wind speed                                                                                                          	| m/s              |

## Band-specific quantities
| Variable             	| Description                                          	| Units           |
|----------------------	|------------------------------------------------------	|---------------- |
| OUT_ADV_SENS_BAND    	| Net sensible heat flux advected to snow pack         	| W/m<sup>2</sup> |
| OUT_ADVECTION_BAND   	| Advected energy                                      	| W/m<sup>2</sup> |
| OUT_ALBEDO_BAND      	| Average surface albedo                               	| fraction        |
| OUT_DELTACC_BAND     	| Change in cold content in snow pack                  	| W/m<sup>2</sup> |
| OUT_GRND_FLUX_BAND   	| Net heat flux into ground                            	| W/m<sup>2</sup> |
| OUT_IN_LONG_BAND     	| Incoming longwave at ground surface (under veg)      	| W/m<sup>2</sup> |
| OUT_LATENT_BAND      	| Net upward latent heat flux                          	| W/m<sup>2</sup> |
| OUT_LATENT_SUB_BAND  	| Net upward latent heat flux due to sublimation       	| W/m<sup>2</sup> |
| OUT_MELT_ENERGY_BAND 	| Energy of fusion (melting) in snowpack               	| W/m<sup>2</sup> |
| OUT_LWNET_BAND      	| Net downward longwave flux                           	| W/m<sup>2</sup> |
| OUT_SWNET_BAND       	| Net downward shortwave flux                      	    | W/m<sup>2</sup> |
| OUT_RFRZ_ENERGY_BAND 	| Net energy used to refreeze liquid water in snowpack 	| W/m<sup>2</sup> |
| OUT_SENSIBLE_BAND    	| Net upward sensible heat flux                        	| W/m<sup>2</sup> |
| OUT_SNOW_CANOPY_BAND 	| Snow interception storage in canopy                  	| mm              |
| OUT_SNOW_COVER_BAND  	| Fractional area of snow cover                        	| fraction        |
| OUT_SNOW_DEPTH_BAND  	| Depth of snow pack                                   	| cm              |
| OUT_SNOW_FLUX_BAND   	| Energy flux through snow pack                        	| W/m<sup>2</sup> |
| OUT_SNOW_MELT_BAND   	| Snow melt                                            	| mm              |
| OUT_SNOW_PACKT_BAND  	| Snow pack temperature                                	| C  	            |
| OUT_SNOW_SURFT_BAND  	| Snow surface temperature                             	| C               |
| OUT_SWE_BAND         	| Snow water equivalent in snow pack                   	| mm              |

## Carbon Cycle Terms
| Variable       	| Description                                  	| Units   	|
|----------------	|----------------------------------------------	|---------	|
| OUT_APAR       	| absorbed PAR                                 	| W/m<sup>2</sup>    	|
| OUT_GPP        	| gross primary productivity                   	| g C/m<sup>2</sup>d 	|
| OUT_RAUT       	| autotrophic respiration                      	| g C/m<sup>2</sup>d 	|
| OUT_NPP        	| net primary productivity                     	| g C/m<sup>2</sup>d 	|
| OUT_LITTERFALL 	| flux of carbon from living biomass into soil 	| g C/m<sup>2</sup>d 	|
| OUT_RHET       	| soil (heterotrophic) respiration             	| g C/m<sup>2</sup>d 	|
| OUT_NEE        	| net ecosystem exchange (=NPP-RHET)           	| g C/m<sup>2</sup>d 	|
| OUT_CLITTER    	| carbon density in litter pool                	| g C/m<sup>2</sup>  	|
| OUT_CINTER     	| carbon density in intermediate pool          	| g C/m<sup>2</sup>  	|
| OUT_CSLOW      	| carbon density in slow pool                  	| g C/m<sup>2</sup>  	|

## Profiling/Timing Terms
| Variable             | Description                    | Units   |
|--------------------- |------------------------------- |-------- |
| OUT_TIME_VICRUN_WALL | Wall time spent inside vic_run | seconds |
| OUT_TIME_VICRUN_CPU  | CPU time spent inside vic_run  | seconds |
