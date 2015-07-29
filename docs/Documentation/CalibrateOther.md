# VIC Model Calibration Using Parameters Other Than Soil or Snow

For grid cells with an overstory the model calculates a windspeed profile through the canopy. The wind profile through the canopy is controlled by

*   the exponential decay through the overstory and the canopy displacement
*   and roughness heights. In general these parameters are fixed for a given
*   vegetation type

Shortwave radiation attenuation trough the canopy is calculated as a function of LAI and an exponential decay coefficient. A value of 0.5 is set as default but values in the range 0.1 to 0.6 have been used in our simulations. Care must be taken so that manipulation of the vegetation attenuation parameters for shortwave radiation and wind reflect real conditions and values. These parameters should be treated as calibration parameters only as needed to correctly describe under-canopy radiation and wind. Both the wind and the shortwave radiation under the canopy should be in the range from 0.1 to 0.5 of the above canopy values.
