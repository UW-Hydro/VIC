# Snow Bands Formulation

Snow elevation bands represent the effect of sub-grid topography on snow accumulation and melt through orographic controls on precipitation and temperature. By including snow bands in areas of high topographic relief the user can obtain a more attenuated hydrograph during spring melt by representing the offset in melt timing between high and low elevations.

The snow band algorithm works as follows:

*   The user specifies a variable number of snow bands with a fractional area and elevation associated with each band. Fractional area must sum to one.
*   Mean pixel temperature from the meteorological forcings file is lapsed to each elevation band using the dry adiabatic lapse rate.
*   The increase in precipitation with elevation is specified through precipitation fractions that sum to one. If the precipitation fraction is equal to the area fraction, precipitation is uniformly distributed across the pixel. If the precipitation fraction exceeds the area fraction, than that fractional area receives a proportionally larger quantity of precipitation.
*   Precipitation falls as snow or rain depending on the lapsed temperature.
*   Specified vegetation fractions are repeated for each snow band. Therefore, the snow model must be run independently for each vegetation type for each snow band.

!!!Note "Regarding Evapotranspiration and Time Step"
    Model run time will increase accordingly

Although the specified number of elevation bands is constant for all pixels in a model run, it is possible to specify bands with zero area fractions, in which case the elevation band will not be run.
