# VIC Model Parameters

The Image Driver uses the [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) file format for its input model parameters. It is possible to convert the VIC ASCII style parameters to this format. We have put together an example ([Tutorial](Ascii_to_NetCDF_params.md) and [Ipython Notebook](https://github.com/UW-Hydro/VIC/blob/develop/samples/notebooks/example_reformat_vic4_parameters_to_vic5image.ipynb)) that provide examples of how to do this conversion. Our example uses the `tonic` [Python](https://www.python.org/) Package.

!!! Note
	It is the user's responsibility to ensure that parameter files are formatted appropriately. Notably, the variables `AreaFract`, `Pfactor`, `zone_fract`, and `Cv` must sum exactly to 1.0. If using the `SNOW_BAND` option, the area weighted `elevation` must match the mean grid cell elevation (`elev`). VIC will print *** warnings *** if any of these criteria are violated.

# Soil Parameters

The Soil Parameters serve three main purposes:

*   Define the cell ID number of each grid cell. This ID number is essentially a database key that links a grid cell to its parameters in the various parameter files.
*   Define the grid cell soil parameters
*   Define initial soil moisture conditions, to be used in the absence of an initial state file.

The soil parameters are supplied to VIC in a NetCDF file, with a separate variable for each soil parameter.

Below is a list of soil parameters.

| Variable Name            | Dimension          | Units    | Type   | Number of Values | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|--------------------------|--------------------|----------|--------|------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| run_cell                 | [lat, lon]         | N/A      | int    | 1                | 1 = Run Grid Cell, 0 = Do Not Run. Must be zero for all grid cells outside of the mask defined in the domain netCDF file.                              |
| gridcel                  | [lat, lon]         | N/A      | int    | 1                | Grid cell number                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| lat                      | [lat, lon]         | degrees  | double | 1                | Latitude of grid cell                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| lon                      | [lat, lon]         | degrees  | double | 1                | Longitude of grid cell                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| infilt                   | [lat, lon]         | N/A      | double | 1                | Variable infiltration curve parameter (binfilt)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| Ds                       | [lat, lon]         | fraction | double | 1                | Fraction of Dsmax where non-linear baseflow begins                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| Dsmax                    | [lat, lon]         | mm/day   | double | 1                | Maximum velocity of baseflow                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| Ws                       | [lat, lon]         | fraction | double | 1                | Fraction of maximum soil moisture where non-linear baseflow occurs                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| c                        | c                  | N/A      | double | 1                | Exponent used in baseflow curve, normally set to 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| expt                     | [nlayer, lat, lon] | N/A      | double | Nlayer           | Exponent n (=3+2/lambda) in Campbell's eqn for hydraulic conductivity, HBH 5.6 (where lambda = soil pore size distribution parameter). Values should be > 3.0.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Ksat                     | [nlayer, lat, lon] | mm/day   | double | Nlayer           | Saturated hydrologic conductivity                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| phi_s                    | [nlayer, lat, lon] | mm/mm    | double | Nlayer           | Soil moisture diffusion parameter                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| init_moist               | [nlayer, lat, lon] | mm       | double | Nlayer           | Initial layer moisture content                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| elev                     | [lat, lon]         | m        | double | 1                | Average elevation of grid cell                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| depth                    | [nlayer, lat, lon] | m        | double | Nlayer           | Thickness of each soil moisture layer                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| avg_T                    | [lat, lon]         | C        | double | 1                | Average soil temperature, used as the bottom boundary for soil heat flux solutions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| dp                       | [lat, lon]         | m        | double | 1                | Soil thermal damping depth (depth at which soil temperature remains constant through the year, ~4 m)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| bubble                   | [nlayer, lat, lon] | cm       | double | 1                | Bubbling pressure of soil. Values should be >0.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| quartz                   | [nlayer, lat, lon] | N/A      | double | Nlayer           | Quartz content of soil                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| bulk_density             | [nlayer, lat, lon] | kg/m3    | double | Nlayer           | Bulk density of soil layer                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| soil_density             | [nlayer, lat, lon] | kg/m3    | double | Nlayer           | Soil particle density, normally 2685 kg/m3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| organic (optional)       | [nlayer, lat, lon] | fraction | double | Nlayer           | Fraction of soil layer that is organic. If ORGANIC_FRACT is TRUE in the global parameter file, this variable must be included in the soil parameter file. If ORGANIC_FRACT is FALSE then this variable must not appear in the soil parameter file. (release 4.1.2 and later)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| bulk_dens_org (optional) | [nlayer, lat, lon] | kg/m3    | double | Nlayer           | Bulk density of organic portion of soil. If ORGANIC_FRACT is TRUE in the global parameter file, this variable must be included in the soil parameter file. If ORGANIC_FRACT is FALSE then this variable must not appear in the soil parameter file. (release 4.1.2 and later)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| soil_dens_org (optional) | [nlayer, lat, lon] | kg/m3    | double | Nlayer           | Soil particle density of organic portion of soil, normally 1300 kg/m3. If ORGANIC_FRACT is TRUE in the global parameter file, this variable must be included in the soil parameter file. If ORGANIC_FRACT is FALSE then this variable must not appear in the soil parameter file. (release 4.1.2 and later)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| off_gmt                  | [lat, lon]         | hours    | double | 1                | Time zone offset from GMT. This parameter determines how VIC interprets sub-daily time steps relative to the model start date and time. We have adopted the following convention for off_gmt:An off_gmt value of 0 indicates that the model start date/time is relative to Greenwich Mean Time (GMT).An off_gmt value of (grid_cell_longitude*24/360) indicates that the model start date/time is relative to local time.When outputting sub-daily results, VIC's output files are referenced to the model start date/time; therefore they are controlled by off_gmt (off_gmt=0 means VIC results are referenced to GMT; off_gmt=(grid_cell_longitude*24/360) means VIC results are referenced to local time). Daily supplied forcings are assumed to start/end at midnight in local time; the forcing start date/time is thus in local time. When VIC disaggregates daily forcings into sub-daily forcings, off_gmt will be used to determine the time lag between the start of the forcing's diurnal cycle and the start of the VIC simulation.Sub-daily supplied forcings are assumed to occur relative to the time zone indicated by off_gmt.  Therefore, if VIC outputs these sub-daily forcings, they will occur at the exact same time of day as in the input files. Therefore, if mixing daily and sub-daily forcing inputs, it is important that any sub-daily forcing inputs be shifted as necessary to be in the time zone indicated by off_gmt. |
| Wcr_fract                | [nlayer, lat, lon] | fraction | double | Nlayer           | Fractional soil moisture content at the critical point (~70% of field capacity) (fraction of maximum moisture)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Wpwp_FRACT               | [nlayer, lat, lon] | fraction | double | Nlayer           | Fractional soil moisture content at the wilting point (fraction of maximum moisture)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| rough                    | [lat, lon]         | m        | double | 1                | Surface roughness of bare soil                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| snow_rough               | [lat, lon]         | m        | double | 1                | Surface roughness of snowpack                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| annual_prec              | [lat, lon]         | mm       | double | 1                | Average annual precipitation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| resid_moist              | [nlayer, lat, lon] | fraction | double | Nlayer           | Soil moisture layer residual moisture                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| fs_active                | [lat, lon]         | 1 or 0   | int    | 1                | If set to 1, then frozen soil algorithm is activated for the grid cell. A 0 indicates that frozen soils are not computed even if soil temperatures fall below 0C.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| frost_slope              | [lat, lon]         | C        | double | 1                | Slope of uniform distribution of soil temperature (if SPATIAL_FROST == TRUE in the global parameter file).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| max_snow_distrib_slope   | [lat, lon]         | m        | double | 1                | Maximum slope of the snow depth distribution. This is only used if SPATIAL_SNOW == TRUE in the global parameter file. This parameter should be set to twice the spatial average snow depth at which coverage == 1.0. In other words, if we define depth_thresh to be the minimum spatial average snow depth below which coverage < 1.0, then max_snow_distrib_slope = 2depth_thresh. NOTE*: Partial snow coverage is only computed when the snow pack has started melting and thespatial average snow pack depth <= max_snow_distrib_slope/2. During the accumulation season, coverage is 1.0. Even after the pack has started melting anddepth <= max_snow_distrib_slope/2, new snowfall resets coverage to 1.0, and the previous partial coverage is stored. Coverage remains at 1.0 until the new snow has melted away, at which point the previous partial coverage is recovered.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| July_Tavg                | [lat, lon]         | C        | double | 1                | Average July air temperature, used for treeline computations (required if COMPUTE_TREELINE == TRUE in the global parameter file).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |

Note: if BASEFLOW is set to NIJSSEN2001 in the [global parameter file](GlobalParam.md), VIC will interpret the baseflow parameters Ds, Dsmax, Ws, and c to be the alternative baseflow parameters D1, D2, D3, and D4.</font="red">

# Vegetation Parameter File

Vegetation parameters needed for the different VIC model set-ups are listed below. The number of vegetation tiles and fraction of grid cell covered are defined for each grid cell.

| Variable Name     | Units     | Description                                   |
|---------------    |-------    |---------------------------------------------  |
| gridcel           | N/A       | Grid cell number                              |
| Nveg              | N/A       | Number of vegetation tiles in the grid cell   |

Repeats for each vegetation tile in the grid cell:

| Variable Name   | Units     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|---------------  |---------- |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  |
| veg_class       | N/A       | Vegetation class identification number (reference index to [vegetation library](../Classic/VegLib.md)) <br><br>*NOTE 1*: it is common practice to define only one tile for each vegetation class in the grid cell. But this is not strictly necessary. It is OK to define multiple tiles having the same vegetation class. <br><br>*NOTE 2*: As of VIC 4.1.1, if you are simulating lakes, you MUST designate one of the tiles from each grid cell as the tile that contains the lake(s). This designation happens in the [lake parameter file](../Classic/LakeParam.md). You can either choose an existing tile to host the lakes, or insert a new tile (just make sure that the sum of the tile areas in the grid cell = 1.0). This extra lake/wetland tile may have the same vegetation class as one of the other existing tiles (see NOTE 1).   |
| Cv              | fraction  | Fraction of grid cell covered by vegetation tile                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |           

For each vegetation tile, repeats for each defined root zone:

| Variable Name   | Units     | Description                                                             |   
|---------------  |---------- |------------------------------------------------------------------------ |
| root_depth      | m         | Root zone thickness (sum of depths is total depth of root penetration)  |
| root_fract      | fraction  | Fraction of root in the current root zone.                              |   

OPTIONAL - If BLOWING_SNOW is TRUE in global parameter file, then for each vegetation tile, the following parameters will also be included:

| Variable Name   | Units   | Description                                                                           |   
|---------------  |-------  |-------------------------------------------------------------------------------------  |
| sigma_slope     | N/A     | Standard deviation of terrain slopes within vegetation tile                           |   
| lag_one         | N/A     | Lag-one autocorrelation of terrain slope within vegetation tile                       |   
| fetch           | m       | Average fetch (distance the wind blows without obstructions) within vegetation tile   |   

OPTIONAL - If VEGPARAM_LAI is TRUE in global parameter file, then for each vegetation tile, the file will also contain the following parameters:

| Variable Name   | Units     | Description                     |   
|---------------  |---------- |-------------------------------- |
| LAI             | fraction  | Leaf Area Index, one per month  |

OPTIONAL - If VEGPARAM_FCAN is TRUE in global parameter file, then for each vegetation tile, the following parameters will also be included:

| Variable Name   | Units     | Description                                       |   
|---------------  |---------- |-------------------------------------------------- |
| FCANOPY         | fraction  | Partial vegetation cover fraction, one per month  |

OPTIONAL - If VEGPARAM_ALBEDO is TRUE in global parameter file, then for each vegetation tile, there following parameters will also be included in the NetCDF parameter file:

| Variable Name | Units    | Description           |
|---------------|----------|-----------------------|
| ALBEDO        | fraction | Albedo, one per month |


# Vegetation Library Parameters

Vegetation parameters are given for different vegetation types. Below are a list of vegetation parameters:

| Variable Name                                                     | Units                 | Number of Values  | Description                                                                                                                                                           |
|-----------------------------------------------------------------  |---------------------  |------------------ |---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| veg_class                                                         | N/A                   | 1                 | Vegetation class identification number (reference index for library table)                                                                                            |
| overstory                                                         | N/A                   | 1                 | Flag to indicate whether or not the current vegetation type has an overstory (TRUE for overstory present [*e.g. trees*], FALSE for overstory not present [*e.g. grass*])  |
| rarc                                                              | s/m                   | 1                 | Architectural resistance of vegetation type (~2 s/m)                                                                                                                  |       
| rmin                                                              | s/m                   | 1                 | Minimum stomatal resistance of vegetation type (~100 s/m)                                                                                                             |       
| LAI                                                               | fraction              | 12                | Leaf-area index of vegetation type                                                                                                                                    |       
| FCANOPY  (Only present if VEGLIB_FCAN=TRUE in [global parameter file](GlobalParam.md) | fraction              | 12                | Partial vegetation cover fraction                                                                                                                                     |       
| albedo                                                            | fraction              | 12                | Shortwave albedo for vegetation type                                                                                                                                  |       
| rough                                                             | m                     | 12                | Vegetation roughness length (typically `0.123 * vegetation height`)                                                                                                     |       
| displacement                                                      | m                     | 12                | Vegetation displacement height (typically `0.67 * vegetation height`)                                                                                                   |       
| wind_h                                                            | m                     | 1                 | Height at which wind speed is measured. **If using snow interception routines please read the [documentation on wind_h](Definitions/#wind_h)**.                                                          |       
| RGL                                                               | W/m<sup>2</sup>                 | 1                 | Minimum incoming shortwave radiation at which there will be transpiration. For trees this is about 30 W/m<sup>2</sup>, for crops about 100 W/m<sup>2</sup>.                               |       
| rad_atten                                                         | fract                 | 1                 | Radiation attenuation factor. Normally set to 0.5, though may need to be adjusted for high latitudes.                                                                 |       
| wind_atten                                                        | fract                 | 1                 | Wind speed attenuation through the overstory. The default value has been 0.5.                                                                                         |       
| trunk_ratio                                                       | fract                 | 1                 | Ratio of total tree height that is trunk (no branches). The default value has been 0.2.                                                                               |       
| Ctype (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | 0 or 1 (0 = C3, 1 = C4) | 1                 | Photosynthetic pathway (NOTE: the previous input values of "C3" and "C4" are still accepted but are deprecated.                                                       |       
| MaxCarboxRate (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | mol CO2/m2s           | 1                 | Maximum carboxlyation rate at 25 deg C                                                                                                                                |       
| MaxETransport (C3) or CO2Specificity (C4) (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | mol CO2/m2s           | 1                 | Maximum electron transport rate at 25 deg C (C3) or CO2 specificity at 25 deg C (C4)                                                                                  |       
| LightUseEff (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | mol CO2/mol photons   | 1                 | Light-use efficiency                                                                                                                                                  |       
| NscaleFlag (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | 0 or 1                | 1                 | 1 = nitrogen-scaling factors are applicable to this veg class                                                                                                         |       
| Wnpp_inhib (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | fract                 | 1                 | Fraction of maximum moisture storage in top soil layer above which photosynthesis begins to be inhibited by wet conditions                                            |       
| NPPfactor_sat (Only present if VEGLIB_PHOTO=TRUE in [global parameter file](GlobalParam.md) | fract                 | 1                 | Ratio of NPP under saturated conditions (moisture fraction = 1.0) to NPP at moisture fraction = Wnpp_inhib                                                            |
| comment                                                           | N/A                   | 1                 | Comment block for vegetation type. Model skips end of line so spaces are valid entries.                                                                                |


# Elevation Bands (Optional)

By default, VIC assumes each grid cell is flat. This assumption can lead to inaccuracies in snow pack estimates in mountainous grid cells. For this reason, the option exists to have VIC break each grid cell up into a number of _elevation bands_ (also called _snow bands_) and simulate them separately. Each band's mean elevation is used to lapse the grid cell average temperature, pressure, and precipitation to a more accurate local estimate.

The NetCDF Parameter file contains information needed to define the properties of each elevation band used by the snow model. Snow elevation bands are used to improve the model's performance in areas with pronounced topography, especially mountainous regions, where the effects of elevation on snow pack accumulation and ablation might be lost in a large grid cell.

The number of snow elevation bands (_option.SNOW_BAND_) to be used with the model is defined in the [global parameter file](GlobalParam.md). The elevation band information is only read if the number of snow elevation bands is greater than 1.

It is not necessary that all grid cells have the same number of elevation bands. _SNOW_BAND_ is simply the maximum number of elevation bands specified anywhere in the domain given by the domain file. For relatively flat grid cells, some of the elevation bands will have _AreaFract_ values of 0\. For these zero-area bands, a value of 0 may be supplied for _elevation_ and _Pfactor_.

Below is a list of elevation band properties:

| Variable Name | Dimension             | Units    | Type   | Number of Values | Description                                                                                                                                                                             |
|---------------|-----------------------|----------|--------|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| cellnum       | [lat, lon]            | N/A      | double | 1                | Grid cell number (should match numbers assigned in soil parameter file)                                                                                                                 |
| AreaFract     | [snow_band, lat, lon] | fraction | double | SNOW_BAND        | Fraction of grid cell covered by each elevation band. Sum of the fractions must equal 1.                                                                                                |
| elevation     | [snow_band, lat, lon] | m        | double | SNOW_BAND        | Mean (or median) elevation of elevation band. This is used to compute the change in air temperature from the grid cell mean elevation.                                                  |
| Pfactor       | [snow_band, lat, lon] | fraction | double | SNOW_BAND        | Fraction of cell precipitation that falls on each elevation band. Total must equal 1. To ignore effects of elevation on precipitation, set these fractions equal to the area fractions. |

# Example netCDF format VIC 5 image driver parameters

```shell
ncdump -h Stehekin_test_params_20160327.nc
netcdf Stehekin_test_params_20160327 {
dimensions:
	lon = 5 ;
	lat = 4 ;
	nlayer = 3 ;
	snow_band = 5 ;
	veg_class = 12 ;
	root_zone = 3 ;
	month = 12 ;
variables:
	double lat(lat) ;
		lat:units = "degrees_north" ;
		lat:long_name = "latitude of grid cell center" ;
	double lon(lon) ;
		lon:units = "degrees_east" ;
		lon:long_name = "longitude of grid cell center" ;
	int layer(nlayer) ;
		layer:long_name = "soil layer" ;
	int run_cell(lat, lon) ;
		run_cell:_FillValue = -2147483647 ;
		run_cell:units = "N/A" ;
		run_cell:description = "1 = Run Grid Cell, 0 = Do Not Run" ;
		run_cell:long_name = "run_cell" ;
	int gridcell(lat, lon) ;
		gridcell:_FillValue = -2147483647 ;
		gridcell:units = "N/A" ;
		gridcell:description = "Grid cell number" ;
		gridcell:long_name = "gridcell" ;
  double lats(lat, lon) ;
      lats:_FillValue = 9.96920996838687e+36 ;
      lats:units = "degrees" ;
      lats:description = "Latitude of grid cell" ;
      lats:long_name = "lats" ;
    double lons(lat, lon) ;
  		lons:_FillValue = 9.96920996838687e+36 ;
  		lons:units = "degrees" ;
  		lons:description = "Longitude of grid cell" ;
  		lons:long_name = "lons" ;
  	double infilt(lat, lon) ;
  		infilt:_FillValue = 9.96920996838687e+36 ;
  		infilt:units = "mm/day" ;
  		infilt:description = "Variable infiltration curve parameter (binfilt)" ;
  		infilt:long_name = "infilt" ;
  	double Ds(lat, lon) ;
  		Ds:_FillValue = 9.96920996838687e+36 ;
  		Ds:units = "fraction" ;
  		Ds:description = "Fraction of Dsmax where non-linear baseflow begins" ;
  		Ds:long_name = "Ds" ;
  	double Dsmax(lat, lon) ;
  		Dsmax:_FillValue = 9.96920996838687e+36 ;
  		Dsmax:units = "mm/day" ;
  		Dsmax:description = "Maximum velocity of baseflow" ;
  		Dsmax:long_name = "Dsmax" ;
  	double Ws(lat, lon) ;
  		Ws:_FillValue = 9.96920996838687e+36 ;
  		Ws:units = "fraction" ;
  		Ws:description = "Fraction of maximum soil moisture where non-linear baseflow occurs" ;
  		Ws:long_name = "Ws" ;
  	double c(lat, lon) ;
  		c:_FillValue = 9.96920996838687e+36 ;
  		c:units = "N/A" ;
  		c:description = "Exponent used in baseflow curve, normally set to 2" ;
  		c:long_name = "c" ;
    double expt(nlayer, lat, lon) ;
    	expt:_FillValue = 9.96920996838687e+36 ;
    	expt:units = "N/A" ;
    	expt:description = "Exponent n (=3+2/lambda) in Campbells eqn for hydraulic conductivity, HBH 5.6 (where lambda = soil pore size distribution parameter).  Values should be > 3.0." ;
    	expt:long_name = "expt" ;
    double Ksat(nlayer, lat, lon) ;
    	Ksat:_FillValue = 9.96920996838687e+36 ;
    	Ksat:units = "mm/day" ;
    	Ksat:description = "Saturated hydrologic conductivity" ;
    	Ksat:long_name = "Ksat" ;
    double phi_s(nlayer, lat, lon) ;
    	phi_s:_FillValue = 9.96920996838687e+36 ;
    	phi_s:units = "mm/mm" ;
    	phi_s:description = "Soil moisture diffusion parameter" ;
    	phi_s:long_name = "phi_s" ;
    double init_moist(nlayer, lat, lon) ;
  		init_moist:_FillValue = 9.96920996838687e+36 ;
  		init_moist:units = "mm" ;
  		init_moist:description = "Initial layer moisture content" ;
  		init_moist:long_name = "init_moist" ;
  	double elev(lat, lon) ;
  		elev:_FillValue = 9.96920996838687e+36 ;
  		elev:units = "m" ;
  		elev:description = "Average elevation of grid cell" ;
  		elev:long_name = "elev" ;
  	double depth(nlayer, lat, lon) ;
  		depth:_FillValue = 9.96920996838687e+36 ;
  		depth:units = "m" ;
  		depth:description = "Thickness of each soil moisture layer" ;
  		depth:long_name = "depth" ;
  	double avg_T(lat, lon) ;
  		avg_T:_FillValue = 9.96920996838687e+36 ;
  		avg_T:units = "C" ;
  		avg_T:description = "Average soil temperature, used as the bottom boundary for soil heat flux solutions" ;
  		avg_T:long_name = "avg_T" ;
    double dp(lat, lon) ;
  		dp:_FillValue = 9.96920996838687e+36 ;
  		dp:units = "m" ;
  		dp:description = "Soil thermal damping depth (depth at which soil temperature remains constant through the year, ~4 m)" ;
  		dp:long_name = "dp" ;
  	double bubble(nlayer, lat, lon) ;
  		bubble:_FillValue = 9.96920996838687e+36 ;
  		bubble:units = "cm" ;
  		bubble:description = "Bubbling pressure of soil. Values should be > 0.0" ;
  		bubble:long_name = "bubble" ;
  	double quartz(nlayer, lat, lon) ;
  		quartz:_FillValue = 9.96920996838687e+36 ;
  		quartz:units = "fraction" ;
  		quartz:description = "Quartz content of soil" ;
  		quartz:long_name = "quartz" ;
  	double bulk_density(nlayer, lat, lon) ;
  		bulk_density:_FillValue = 9.96920996838687e+36 ;
  		bulk_density:units = "kg/m3" ;
  		bulk_density:description = "Bulk density of soil layer" ;
  		bulk_density:long_name = "bulk_density" ;
  	double soil_density(nlayer, lat, lon) ;
  		soil_density:_FillValue = 9.96920996838687e+36 ;
  		soil_density:units = "kg/m3" ;
  		soil_density:description = "Soil particle density, normally 2685 kg/m3" ;
  		soil_density:long_name = "soil_density" ;
  	double off_gmt(lat, lon) ;
  		off_gmt:_FillValue = 9.96920996838687e+36 ;
  		off_gmt:units = "hours" ;
  		off_gmt:description = "Time zone offset from GMT. This parameter determines how VIC interprets sub-daily time steps relative to the model start date and time." ;
  		off_gmt:long_name = "off_gmt" ;
  	double Wcr_FRACT(nlayer, lat, lon) ;
  		Wcr_FRACT:_FillValue = 9.96920996838687e+36 ;
  		Wcr_FRACT:units = "fraction" ;
  		Wcr_FRACT:description = "Fractional soil moisture content at the critical point (~70%% of field capacity) (fraction of maximum moisture)" ;
  		Wcr_FRACT:long_name = "Wcr_FRACT" ;
    double Wpwp_FRACT(nlayer, lat, lon) ;
  		Wpwp_FRACT:_FillValue = 9.96920996838687e+36 ;
  		Wpwp_FRACT:units = "fraction" ;
  		Wpwp_FRACT:description = "Fractional soil moisture content at the wilting point (fraction of maximum moisture)" ;
  		Wpwp_FRACT:long_name = "Wpwp_FRACT" ;
  	double rough(lat, lon) ;
  		rough:_FillValue = 9.96920996838687e+36 ;
  		rough:units = "m" ;
  		rough:description = "Surface roughness of bare soil" ;
  		rough:long_name = "rough" ;
  	double snow_rough(lat, lon) ;
  		snow_rough:_FillValue = 9.96920996838687e+36 ;
  		snow_rough:units = "m" ;
  		snow_rough:description = "Surface roughness of snowpack" ;
  		snow_rough:long_name = "snow_rough" ;
  	double annual_prec(lat, lon) ;
  		annual_prec:_FillValue = 9.96920996838687e+36 ;
  		annual_prec:units = "mm" ;
  		annual_prec:description = "Average annual precipitation." ;
  		annual_prec:long_name = "annual_prec" ;
  	double resid_moist(nlayer, lat, lon) ;
  		resid_moist:_FillValue = 9.96920996838687e+36 ;
  		resid_moist:units = "fraction" ;
  		resid_moist:description = "Soil moisture layer residual moisture." ;
  		resid_moist:long_name = "resid_moist" ;
  	int fs_active(lat, lon) ;
  		fs_active:_FillValue = -2147483647 ;
  		fs_active:units = "binary" ;
  		fs_active:description = "If set to 1, then frozen soil algorithm is activated for the grid cell. A 0 indicates that frozen soils are not computed even if soil temperatures fall below 0C." ;
  		fs_active:long_name = "fs_active" ;
  	int snow_band(snow_band) ;
  		snow_band:long_name = "snow band" ;
        	double cellnum(lat, lon) ;
        		cellnum:_FillValue = 9.96920996838687e+36 ;
        		cellnum:units = "N/A" ;
        		cellnum:description = "Grid cell number (should match numbers assigned in soil parameter file)" ;
    double AreaFract(snow_band, lat, lon) ;
  		AreaFract:_FillValue = 9.96920996838687e+36 ;
  		AreaFract:units = "fraction" ;
  		AreaFract:description = "Fraction of grid cell covered by each elevation band. Sum of the fractions must equal 1." ;
  	double elevation(snow_band, lat, lon) ;
  		elevation:_FillValue = 9.96920996838687e+36 ;
  		elevation:units = "m" ;
  		elevation:description = "Mean (or median) elevation of elevation band. This is used to compute the change in air temperature from the grid cell mean elevation." ;
  	double Pfactor(snow_band, lat, lon) ;
  		Pfactor:_FillValue = 9.96920996838687e+36 ;
  		Pfactor:units = "fraction" ;
  		Pfactor:description = "Fraction of cell precipitation thatfalls on each elevation band. Total must equal 1. To ignore effects of elevation on precipitation, set these fractions equal to the area fractions." ;
  	int veg_class(veg_class) ;
  		veg_class:long_name = "Vegetation Class" ;
  	string veg_descr(veg_class) ;
  		veg_descr:long_name = "Vegetation Class Description" ;
  	int root_zone(root_zone) ;
  		root_zone:long_name = "root zone" ;
  	int month(month) ;
  		month:long_name = "month of year" ;
  	int Nveg(lat, lon) ;
  		Nveg:_FillValue = -2147483647 ;
  		Nveg:long_name = "Nveg" ;
  		Nveg:units = "N/A" ;
  		Nveg:description = "Number of vegetation tiles in the grid cell" ;
  	double Cv(veg_class, lat, lon) ;
  		Cv:_FillValue = 9.96920996838687e+36 ;
  		Cv:long_name = "Cv" ;
  		Cv:units = "fraction" ;
  		Cv:description = "Fraction of grid cell covered by vegetation tile" ;
  	double root_depth(veg_class, root_zone, lat, lon) ;
  		root_depth:_FillValue = 9.96920996838687e+36 ;
  		root_depth:long_name = "root_depth" ;
  		root_depth:units = "m" ;
  		root_depth:description = "Root zone thickness (sum of depths is total depth of root penetration)" ;
    double root_fract(veg_class, root_zone, lat, lon) ;
  		root_fract:_FillValue = 9.96920996838687e+36 ;
  		root_fract:long_name = "root_fract" ;
  		root_fract:units = "fraction" ;
  		root_fract:description = "Fraction of root in the current root zone" ;
  	double LAI(veg_class, month, lat, lon) ;
  		LAI:_FillValue = 9.96920996838687e+36 ;
  		LAI:long_name = "LAI" ;
  		LAI:units = "m2/m2" ;
  		LAI:description = "Leaf Area Index, one per month" ;
  	int overstory(veg_class, lat, lon) ;
  		overstory:_FillValue = -2147483647 ;
  		overstory:long_name = "overstory" ;
  		overstory:units = "N/A" ;
  		overstory:description = "Flag to indicate whether or not the current vegetation type has an overstory (1 for overstory present [e.g. trees], 0 for overstory not present [e.g. grass])" ;
  	double rarc(veg_class, lat, lon) ;
  		rarc:_FillValue = 9.96920996838687e+36 ;
  		rarc:long_name = "rarc" ;
  		rarc:units = "s/m" ;
  		rarc:description = "Architectural resistance of vegetation type (~2 s/m)" ;
  	double rmin(veg_class, lat, lon) ;
  		rmin:_FillValue = 9.96920996838687e+36 ;
  		rmin:long_name = "rmin" ;
  		rmin:units = "s/m" ;
  		rmin:description = "Minimum stomatal resistance of vegetation type (~100 s/m)" ;
  	double wind_h(veg_class, lat, lon) ;
  		wind_h:_FillValue = 9.96920996838687e+36 ;
  		wind_h:long_name = "wind_h" ;
  		wind_h:units = "m" ;
  		wind_h:description = "Height at which wind speed is measured." ;
  	double RGL(veg_class, lat, lon) ;
  		RGL:_FillValue = 9.96920996838687e+36 ;
  		RGL:long_name = "RGL" ;
  		RGL:units = "W/m^2." ;
  		RGL:description = "Minimum incoming shortwave radiation at which there will be transpiration. For trees this is about 30 W/m^2, for crops about 100 W/m^2." ;
    double rad_atten(veg_class, lat, lon) ;
  		rad_atten:_FillValue = 9.96920996838687e+36 ;
  		rad_atten:long_name = "rad_atten" ;
  		rad_atten:units = "fraction" ;
  		rad_atten:description = "Radiation attenuation factor. Normally set to 0.5, though may need to be adjusted for high latitudes." ;
  	double wind_atten(veg_class, lat, lon) ;
  		wind_atten:_FillValue = 9.96920996838687e+36 ;
  		wind_atten:long_name = "wind_atten" ;
  		wind_atten:units = "fraction" ;
  		wind_atten:description = "Wind speed attenuation through the overstory. The default value has been 0.5." ;
  	double trunk_ratio(veg_class, lat, lon) ;
  		trunk_ratio:_FillValue = 9.96920996838687e+36 ;
  		trunk_ratio:long_name = "trunk_ratio" ;
  		trunk_ratio:units = "fraction" ;
  		trunk_ratio:description = "Ratio of total tree height that is trunk (no branches). The default value has been 0.2." ;
  	double albedo(veg_class, month, lat, lon) ;
  		albedo:_FillValue = 9.96920996838687e+36 ;
  		albedo:long_name = "albedo" ;
  		albedo:units = "fraction" ;
  		albedo:description = "Shortwave albedo for vegetation type" ;
  	double veg_rough(veg_class, month, lat, lon) ;
  		veg_rough:_FillValue = 9.96920996838687e+36 ;
  		veg_rough:long_name = "veg_rough" ;
  		veg_rough:units = "m" ;
  		veg_rough:description = "Vegetation roughness length (typically 0.123 * vegetation height)" ;
  	double displacement(veg_class, month, lat, lon) ;
  		displacement:_FillValue = 9.96920996838687e+36 ;
  		displacement:long_name = "displacement" ;
  		displacement:units = "m" ;
  		displacement:description = "Vegetation displacement height (typically 0.67 * vegetation height)" ;
      // global attributes:
      		:description = "VIC parameter file" ;
      		:history = "Created: Mon Mar 28 09:32:46 2016\n/Users/jhamman/anaconda/lib/python3.4/site-packages/ipykernel/__main__.py -f /Users/jhamman/Library/Jupyter/runtime/kernel-b8627880-fe9e-4d1a-b2d7-fa3fb6a23e22.json\n" ;
      		:source = "/Users/jhamman/anaconda/lib/python3.4/site-packages/ipykernel/__main__.py" ;
      		:username = "jhamman" ;
      		:host = "Joes-iMac.local" ;
      }
```
