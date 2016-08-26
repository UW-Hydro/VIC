# Snow Model

VIC uses an energy balance approach to represent snow accumulation and ablation on the ground ([Andreadis et al., 2009](../Documentation/References.md#primary-historical-reference)). Although low-lying vegetation is assumed to be completely covered by snowfall, and hence does not affect the ground snow pack energy balance, the model does contain an explicit canopy snow interception scheme that accounts for sublimation, drip and release of intercepted snow. Additionally, each grid cell is subdivided into elevation bands were the user defines the precipitation fraction in each band.

## Ground Snowpack

Ground snow accumulation and melt are simulated using a two-layer energy-balance model at the snow surface, similar to that described by [Anderson (1968)](http://dx.doi.org/10.1029/WR004i001p00019). The snowpack is divided into two layers (a thin surface layer and the pack layer) and all the important heat and energy fluxes are considered (longwave and shortwave radiation, sensible and latent heat, convective energy). Internal energy of the snowpack is also considered. The ground heat flux is ignored (unless the frozen soil model is used). Water can be added to the snowpack as rain, snow, or drip/throughfall from the canopy. If snow is present it is assumed to completely cover the ground, thereby affecting radiation transfer and the wind profiles via increased albedo and decreased surface roughness (the snow surface roughness is used).

In every time step the model calculates the rain or snow fraction that is added to the snowpack. Then all the energy fluxes are calculated and if the energy balance is positive melt occurs. If the liquid water holding capacity of both the surface and pack layers are exceeded then the excess liquid water is immediately released as snowpack outflow. If the energy balance is negative, then the energy balance is solved by iterating on the snow surface temperature.

## Intercepted snow

Snow can be intercepted and stored by the canopy. The amount of snow that is intercepted is calculated as a function of leaf area index (LAI). Maximum water storage in the canopy is a function of wind and air temperature as well as LAI.

To reduce computational expense, snowmelt from the canopy is calculated using a simplified version of the ground snowpack energy balance. If the air temperature is below freezing, then the snow surface temperature in the canopy is set to the ambient air temperature, otherwise the snow surface temperature is set to zero.

## Calibrating the snow model

During calibration three main parameters are adjusted for grid cells without overlying vegetation.

1.  Maximum air temperature at which snowfall occurs
2.  Minimum air temperature at which rainfall occurs
3.  Snow surface roughness

Usually the two first parameters are set to 1.5 and -0.5, respectively. Our experience suggests that the snow roughness parameter should be in the range from 0.001 m to 0.03 m.

For grid cells with an overstory the model calculates a windspeed profile through the canopy. The wind profile through the canopy is controlled by the exponential decay through the overstory and the canopy displacement and roughness heights. In general these parameters are fixed for a given vegetation type.

Shortwave radiation attenuation trough the canopy is calculated as a function of LAI and an exponential decay coefficient. A value of 0.5 is set as default but values in the range 0.1 to 0.6 have been used in our simulations. Care must be taken so that manipulation of the vegetation attenuation parameters for shortwave radiation and wind reflect real conditions and values. These parameters should be treated as calibration parameters only as needed to correctly describe under-canopy radiation and wind.
