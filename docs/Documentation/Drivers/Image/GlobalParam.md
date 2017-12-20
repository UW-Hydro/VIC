# VIC Run-Time Options - Global Parameter File

The global parameter file serves two main purposes:

1.  Tells VIC the names, locations, and formats of input and output files
2.  Defines global parameters of the simulation (known as _run-time_ options)

The order of the options in the global parameter file is not important, but the complete option name must be followed by the required option type information. To help in understanding this file, an [example file](#example-global-parameter-file) has been attached at the bottom of this page.

# Define Simulation Parameters

The following options determine the type of simulation that will be performed.

## Main Simulation Parameters

| Name                 | Type    | Units  | Description                                                                                                                                                                                      |
|----------------------|---------|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| NLAYER               | integer | N/A    | Number of   moisture layers used by the model                                                                                                                                                    |
| NODES                | integer | N/A    | Number of   thermal solution nodes in the soil column                                                                                                                                            |
| MODEL_STEPS_PER_DAY  | integer | steps  | Number of   simulation time steps per day. NOTE: MODEL_STEPS_PER_DAY should be > 4 for   FULL_ENERGY=TRUE or FROZEN_SOIL=TRUE.                                                                   |
| SNOW_STEPS_PER_DAY   | integer | steps  | Number of   time steps per day used to solve the snow model (if MODEL_STEPS_PER_DAY >   1, SNOW_STEPS_PER_DAY should = MODEL_STEPS_PER_DAY)                                                      |
| RUNOFF_STEPS_PER_DAY | integer | steps  | Number of   time steps per day used to solve the runoff model (should be >=   MODEL_STEPS_PER_DAY)                                                                                               |
| STARTYEAR            | integer | year   | Year   model simulation starts. **NOTE**: STARTYEAR, STARTMONTH, STARTDAY and STARTSEC togetehr specify the begenning time point of the first simulation time step. |
| STARTMONTH           | integer | month  | Month   model simulation starts                                                                                                                                                                  |
| STARTDAY             | integer | day    | Day model   simulation starts                                                                                                                                                                    |
| STARTSEC             | integer | second | Second   model simulation starts (e.g., STARTSEC=0 for starting from the beginning of   a day; STARTSEC=3600 for starting from the beginning of the second hour of a   day) <br><br>Default = 0. |
| **EITHER:**          |         |        | *Note:*   **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but   **not both**                                                                                                |
| NRECS                | integer | N/A    | Number of   time steps over which to run model. ***NOTE: The number of records must be   defined such that the model completes the last day.***                                                  |
| **OR:**              |         |        | *Note:*   **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but   **not both**                                                                                                |
| ENDYEAR              | integer | year   | Year   model simulation ends. **NOTE**: ENDYEAR, ENDMONTH and ENDDAY together specify the last day on which the simulation runs. For sub-daily model run timestep, VIC always runs until the last timestep of this day. |
| ENDMONTH             | integer | month  | Month   model simulation ends                                                                                                                                                                    |
| ENDDAY               | integer | day    | Day model   simulation ends                                                                                                                                                                      |
| CALENDAR             | string  | N/A    | Calendar   to use. Valid calendars: 'standard', 'gregorian', 'proleptic_gregorian'   'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'. Must match the time variable type in the forcing netCDF file. |
| OUT_TIME_UNITS       | string  | N/A    | Units for output time variables. Valid options: 'SECONDS', 'MINUTES', 'HOURS', 'DAYS'. <br><br>Default = DAYS.                                                                                       |


## Define Energy Balance Parameters

The following options determine the method of resolving the surface energy balance.

| Name          | Type      | Units             | Description                                                                                                                                                                                                                                                                                   |
|-------------- |--------   |---------------    |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    |
| FULL_ENERGY   | string    | TRUE or FALSE     | Option for computing land surface temperature (soil or snowpack surface). <li>**TRUE** = compute (via iteration) the temperature that balances the surface energy budget.  <li>**FALSE** = set surface temperature equal to air temperature.  <br><br>Default = False.                                                |
| CLOSE_ENERGY  | string    | TRUE or FALSE     | Option for controlling links between the energy balances of the surface and the canopy. <li>**TRUE** = iterate between the canopy and surface energy balances until they are consistent. <li>**FALSE** = compute the surface and canopy energy balances separately, once per time step.  <br><br>Default = FALSE.     |

## Define Soil Temperature Parameters

The following options determine the type of simulation that will be performed.

| Name              | Type              | Units                              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|-------------------|-------------------|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FROZEN_SOIL       | string            | TRUE or FALSE                      | Option for handling the water/ice phase change in frozen soils. TRUE = account for water/ice phase change (including latent heat). FALSE = soil moisture always remains liquid, even when below 0 C; no latent heat effects and ice content is always 0. Default = FALSE. Note: to activate this option, the user must also set theFS_ACTIVE flag to 1 in the soil parameter file for each grid cell where this option is desired. In other words, the user can choose for some grid cells (e.g. cold ones) to compute ice contents and for others (e.g. warm ones) to skip the extra computation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| QUICK_FLUX        | string            | TRUE or FALSE                      | Option for computing the soil vertical temperature profile. TRUE = use the approximate method described by Liang et al. (1999) to compute soil temperatures and ground heat flux; this method ignores water/ice phase changes. FALSE = use the finite element method described in Cherkauer and Lettenmaier (1999) to compute soil temperatures and ground heat flux; this method is appropriate for accounting for water/ice phase changes. Default = FALSE (i.e. use Cherkauer and Lettenmaier (1999)) when running FROZEN_SOIL; and TRUE (i.e. use Liang et al. (1999)) in all other cases.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| IMPLICIT          | string            | TRUE or FALSE                      | If TRUE the model will use an implicit solution for the soil heat flux equation of Cherkauer and Lettenmaier (1999)(QUICK_FLUX is FALSE), otherwise uses original explicit solution. When QUICK_FLUX is TRUE the implicit solution has no effect. The user can override this option by setting IMPLICIT to FALSE in the global parameter file. The implicit solution is guaranteed to be stable for all combinations of time step and thermal node spacing; the explicit solution is only stable for some combinations. If the user sets IMPLICIT to FALSE, VIC will check the time step, node spacing, and soil thermal properties to confirm stability. If the explicit solution will not be stable, VIC will exit with an error message. Default = TRUE.                                                                                                                                                                                                                                                                                                                                         |
| QUICK_SOLVE       | string            | TRUE or FALSE                      | This option is a hybrid of QUICK_FLUX TRUE and FALSE. If TRUE model will use the method described by Liang et al. (1999)to compute ground heat flux during the surface energy balance iterations, and then will use the method described in Cherkauer and Lettenmaier (1999) for the final solution step. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| NOFLUX            | string            | TRUE or FALSE                      | If TRUE model will use a no flux bottom boundary with the finite difference soil thermal solution (i.e. QUICK_FLUX = FALSE or FULL_ENERGY = TRUE or FROZEN_SOIL = TRUE). Default = FALSE (i.e., use a constant temperature bottom boundary condition).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| EXP_TRANS         | string            | TRUE or FALSE                      | If TRUE the model will exponentially distributes the thermal nodes in the Cherkauer and Lettenmaier (1999) finite difference algorithm, otherwise uses linear distribution. (This is only used if FROZEN_SOIL = TRUE). Default = TRUE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| GRND_FLUX_TYPE    | string            | N/A                                | Options for handling ground flux:GF_406 = use (flawed) formulas for ground flux, deltaH, and fusion as in VIC 4.0.6 and earlier.GF_410 = use formulas from VIC 4.1.0. NOTE: this option exists for backwards compatibility with earlier releases and likely will be removed in later releases. Default = GF_410.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| TFALLBACK         | string            | TRUE or FALSE                      | Options for handling failures of T iterations to converge. FALSE = if T iteration fails to converge, report an error. TRUE = if T iteration fails to converge, use the previous time step's T value. This option affects the temperatures of canopy air, canopy snow, ground snow pack, ground surface, and soil T nodes. If TFALLBACK is TRUE, VIC will report the total number of instances in which the previous step's T was used, at the end of each grid cell's simulation. In addition, a time series of when these instances occurred (averaged across all veg tile/snow band combinations) can be written to the output files, using the following output variables:OUT_TFOL_FBFLAG = time series of T fallbacks in canopy snow T solution.OUT_TCAN_FBFLAG = time series of T fallbacks in canopy air T solution. OUT_SNOWT_FBFLAG = time series of T fallbacks in snow pack surface T solution.OUT_SURFT_FBFLAG = time series of T fallbacks in ground surface T solution.OUT_SOILT_FBFLAG = time series of T fallbacks in soil node T solution (one time series per node). Default = TRUE. |
| SHARE_LAYER_MOIST | string            | TRUE or FALSE                      | If TRUE, then *if* the soil moisture in the layer that contains more than half of the roots is above the critical point, then the plant's roots in the drier layers can access the moisture of the wetter layer so that the plant does not experience moisture limitation. <br> If FALSE or all of the soil layer moistures are below the critical point, transpiration in each layer is limited by the layer's soil moisture. <br><br> Default: TRUE.  |
| SPATIAL_FROST     | string (+integer) | string: TRUE or FALSE integer: N/A | Option to allow spatial heterogeneity in soil temperature:FALSE = Assume soil temperature is horizontally constant (only varies with depth). TRUE = Assume soil temperatures at each given depth are distributed horizontally with a uniform (linear) distribution, so that even when the mean temperature is below freezing, some portion of the soil within the grid cell at that depth could potentially be above freezing. This requires specifying a frost slope value as an extra field in the soil parameter file, so that the minimum/maximum temperatures can be computed from the mean value. The maximum and minimum temperatures will be set to mean temperature +/- frost_slope.If TRUE is specified, you must follow this with an integer value for Nfrost, the number of frost sub-areas (each having a distinct temperature). Default = FALSE.                                                                                                                                                                                                                                       |

## Precipitation (Rain and Snow) Parameters

Generally these default values do not need to be overridden.

| Name                  | Type              | Units                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|-----------------------|-------------------|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| SNOW_DENSITY          | string            | N/A                   | Options for computing snow density:DENS_BRAS = Use traditional VIC algorithm taken from Bras, 1990.DENS_SNTHRM = Use algorithm taken from SNTHRM model.Default = DENS_BRAS.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| BLOWING               | string            | TRUE or FALSE         | If TRUE, compute evaporative fluxes due to blowing snow. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| BLOWING_VAR_THRESHOLD | string            | TRUE or FALSE   | If TRUE, a variable shear stress threshold is used to determine the blowing snow flux. If FALSE, then a fixed threshold is used. See Li and Pomeroy (1997) for details. <br><br>Default: TRUE. |
| BLOWING_CALC_PROB     | string            | TRUE or FALSE  | If TRUE, the probability of occurrence of blowing snow is calculated as a function of environmental conditions. If FALSE, then the probability is set to 1. See Lu and Pomeroy (1997) for details. <br><br>Default: TRUE.  |
| BLOWING_SIMPLE        | string            | TRUE or FALSE   | If TRUE, the sublimation flux of blowing snow is calculated as a function vapor pressure and wind speed. If FALSE, then additional calculations are made to account for a saltation and suspension layer. See Lu and Pomeroy (1997) for details. <br><br>Default: FALSE. |
| BLOWING_FETCH         | string            | TRUE or FALSE   | This option is only used when BLOWING_SIMPLE is set to FALSE. When this option is set to TRUE, the fetch is accounted for in the calculation of the sublimation flux from blowing snow. If FALSE then the fetch is not used. See Lu and Pomeroy (1997) for details. <br><br> Default: TRUE. |
| BLOWING_SPATIAL_WIND  | string            | TRUE or FALSE   | If TRUE, multiple wind speed ranges, calculated according to a probability distribution, are used to determine the sublimation flux from blowing snow. If FALSE, then a single wind speed is used. See Lu and Pomeroy (1997) for details. <br><br>Default: TRUE. |
| COMPUTE_TREELINE      | string or integer | FALSE or veg class id | Options for handling above-treeline vegetation:FALSE = Do not compute treeline or replace vegetation above the treeline.CLASS_ID = Compute the treeline elevation based on average July temperatures; for those elevation bands with elevations above the treeline (or the entire grid cell if SNOW_BAND == 1 and the grid cell elevation is above the tree line), if they contain vegetation tiles having overstory, replace that vegetation with the vegetation having id CLASS_ID in the vegetation library. NOTE 1: You MUST supply VIC with a July average air temperature, in the optional July_Tavg field, AND set theJULY_TAVG_SUPPLIED option to TRUE so that VIC can read the soil parameter file correctly. NOTE 2: If LAKES=TRUE, COMPUTE_TREELINE MUST be FALSE.Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| CORRPREC              | string            | TRUE or FALSE         | If TRUE correct precipitation for gauge undercatch. NOTE: This option is not supported when using snow/elevation bands. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| MAX_SNOW_TEMP         | float             | deg C                 | Maximum temperature at which snow can fall. Default = 0.5 C.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| MIN_RAIN_TEMP         | float             | deg C                 | Minimum temperature at which rain can fall. Default = -0.5 C.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| SPATIAL_SNOW          | string            | TRUE or FALSE         | Option to allow spatial heterogeneity in snow water equivalent (yielding partial snow coverage) when the snow pack is melting:FALSE = Assume snow water equivalent is constant across grid cell. TRUE = Assume snow water equivalent is distributed horizontally with a uniform (linear) distribution, so that some portion of the grid cell has 0 snow pack. This requires specifying the max_snow_distrib_slope value as an extra field in the soil parameter file. NOTE: max_snow_distrib_slope should be set to twice the desired minimum spatial average snow pack depth [m]. I.e., if we define depth_thresh to be the minimum spatial average snow depth below which coverage < 1.0, then max_snow_distrib_slope = 2*depth_thresh. NOTE: Partial snow coverage is only computed when the snow pack has started melting and the spatial average snow pack depth <= max_snow_distrib_slope/2. During the accumulation season, coverage is 1.0. Even after the pack has started melting and depth <= max_snow_distrib_slope/2, new snowfall resets coverage to 1.0, and the previous partial coverage is stored. Coverage remains at 1.0 until the new snow has melted away, at which point the previous partial coverage is recovered. Default = FALSE. |

## Turbulent Flux Parameters

Generally these default values do not need to be overridden.

| Name                  | Type      | Units     | Description                                                                                                                                                                                   |
|---------------------  |--------   |-------    |--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   |
| MIN_WIND_SPEED        | float     | m/s       | Minimum allowable wind speed. <br><br>Default = 0.1 m/s. |
| AERO_RESIST_CANSNOW   | string    | N/A       | Options for aerodynamic resistance in snow-filled canopy:  <li>**AR_406** = Multiply by 10 for latent heat, but do NOT multiply by 10 for sensible heat. When no snow in canopy, use surface aero_resist instead of overstory aero_resist. (As in VIC 4.0.6).  <li>**AR_406_LS** = Multiply by 10 for both latent and sensible heat. When no snow in canopy, use surface aero_resist instead of overstory aero_resist.  <li>**AR_406_FULL** = Multiply by 10 for both latent and sensible heat. Always use overstory aero_resist (snow or bare).  <li>**AR_410** = Apply stability correction, instead of multiplying by 10, for both latent and sensible heat. Always use overstory aero_resist (snow or bare). <br><br>*NOTE*: this option exists for backwards compatibility with earlier releases and likely will be removed in later releases. <br><br>Default = AR_406_FULL. |

## Meteorological Forcing Disaggregation Parameters

Generally these default values do not need to be overridden.

| Name            | Type   | Units         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|-----------------|--------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| PLAPSE          | string | TRUE or FALSE | Options for computing grid cell average surface atmospheric pressure (and density) when it is not explicitly supplied as a meteorological forcing:  <li>**FALSE** = Set surface atmospheric pressure to constant 95.5 kPa (as in earlier releases).  <li>**TRUE** = Lapse surface atmospheric pressure (and density) from sea level to the grid cell average elevation.  <br><br>*NOTE 1*: air pressure is already lapsed to grid cell or band elevation when computing latent heat; this option only affects computation of sensible heat.  <br><br>*NOTE 2*: this option exists for backwards compatibility with earlier releases and likely will be removed in later releases (the TRUE option will become the standard behavior).  <br><br>Default = TRUE. |
| SW_PREC_THRESH  | float  | mm            | Minimum daily precipitation, above which incoming shortwave is dimmed by 25%, when shortwave is not supplied as a forcing but instead is estimated from daily temperature range. <br><br>*Note*: This option's purpose is to avoid erroneous dimming of estimated shortwave when using forcings that have been aggregated or re-sampled from a different resolution. Re-sampling can sometimes smear small amounts of precipitation from neighboring cells into cells that originally had no precipitation. The appropriate value of SW_PREC_THRESH must be found through examination of the forcings.  <br><br>Default = 0 mm (any precipitation causes dimming)  |
| MTCLIM_SWE_CORR | string | TRUE or FALSE | This controls VIC's estimates of incoming shortwave (when shortwave is not supplied as a forcing) in the presence of snow. When shortwave is supplied as a forcing, this option is ignored.  <li>**TRUE** = Adjust incoming shortwave for snow albedo effect.  <li>**FALSE** = Do not adjust shortwave (as in earlier releases).  <br><br>Default = TRUE.  |
| VP_ITER         | string | N/A           | This controls VIC's iteration between estimates of shortwave and vapor pressure:  <li>**VP_ITER_NEVER** = Never iterate; make estimates separately.  <li>**VP_ITER_ALWAYS** = Always iterate once (as in previous releases).  <li>**VP_ITER_ANNUAL** = Iterate once for arid climates (based on annual Precip/PET ratio) and never for humid climates.  <li>**VP_ITER_CONVERGE** = Always iterate until shortwave and vp stabilize.  <br><br>Default = VP_ITER_ALWAYS.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| VP_INTERP       | string | TRUE or FALSE | This controls sub-daily humidity estimates:  <li>**TRUE** = Interpolate daily VP estimates linearly between sunrise of one day to the next.  <li>**FALSE** = Hold VP constant for entire day (as in previous releases).  <br><br>Default = TRUE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| LW_TYPE         | string | N/A           | This controls the algorithm used to estimate clear-sky longwave radiation:  <li>**LW_TVA** = Tennessee Valley Authority algorithm (1972) (this option is what previous releases used)  <li>**LW_ANDERSON** = Algorithm of Anderson (1964)  <li>**LW_BRUTSAERT** = Algorithm of Brutsaert (1975)  <li>**LW_SATTERLUND** = Algorithm of Satterlund (1979)  <li>**LW_IDSO** = Algorithm of Idso (1981)  <li>**LW_PRATA** = Algorithm of Prata (1996)  <br><br>Default = LW_TVA.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| LW_CLOUD        | string | N/A           | This controls the algorithm used to estimate the influence of clouds on total longwave:  <li>**LW_CLOUD_BRAS** = Method from Bras textbook (this option is what previous releases used)  <li>**LW_CLOUD_DEARDORFF** = Algorithm of Deardorff (1978)  <br><br>Default = LW_CLOUD_DEARDORFF.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |

## Carbon Parameters

The following options only apply to carbon cycling.

| Name          | Type              | Units                     | Description                                                                                                                                                                           |
|-------------- |---------------    |-----------------------    |-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    |
| CARBON        | string            | TRUE or FALSE             | Options for handling carbon cycle: <li>**FALSE** = do not simulate carbon cycle <li>**TRUE** = simulate carbon cycle <br><br>Default = FALSE.                                                                                                                                                                      |
| RC_MODE       | string            | RC_JARVIS or RC_PHOTO     | Determines how canopy resistance is computed. Options for RC_MODE: <li>**RC_JARVIS** = VIC computes canopy resistance by applying resistance factors to the veg class's minimum canopy resistance listed in the veg library file. <li>**RC_PHOTO** = VIC computes canopy resistance by applying resistance factors to the canopy resistance corresponding to photosynthetic demand (in the absence of moisture limitation). <br><br>Default = RC_JARVIS.                                                                                                                                                                  |
| VEGLIB_PHOTO  | TRUE or FALSE     | string                    | Tells VIC about the contents of the veg library file. Options for VEGLIB_PHOTO: <li>**FALSE** = veg library file does not contain photosynthesis parameters. <li>**TRUE** = veg library file contains photosynthesis parameters. <br><br>Default = FALSE                                                                                                                                                                       |

## Miscellaneous Parameters

Generally these default values do not need to be overridden.

| Name              | Type      | Units             | Description                                                                                                                                                                                                                                                                                                                                                               |
|-----------------  |--------   |---------------    |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  |
| CONTINUEONERROR   | string    | TRUE or FALSE     | Options for handling fatal errors:. <li>**FALSE** = if simulation of a grid cell encounters an error, exit VIC. <li>**TRUE** = if simulation of a grid cell encounters an error, move to next grid cell. <br><br>*NOTE*: in either case, if a grid cell encounters a fatal error, the output files for that grid cell will likely be incomplete. But since most fatal errors are the result of failure of the temperature iteration to converge, seting the TFALLBACK option to TRUE should eliminate most fatal errors. See the section on Soil Temperature Options for more information.. <br><br>Default = TRUE.                                                                                                                                                                                                                                                                                                                                                           |

# Define State Files

The following options control input and output of state files.

| Name         | Type    | Units         | Description                                                                                                                                                                                                                                                                                 |
|--------------|---------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| INIT_STATE   | string  | path/filename | Full path and filename of initial state netCDF file. NOTE: if INIT_STATE is not specified, VIC will take initial soil moistures from the soil parameter file and set all other state variables to a default state.                                                                                 |
| STATENAME    | string  | path/filename | Path and file prefix of the state file to be created on the specified date. The date and second within the simulation at which the state is saved (as well as ".nc") will be appended to the file prefix to form a complete file name. *NOTE*: if STATENAME is not specified, VIC will not save its state in a statefile. |
| STATEYEAR    | integer | year          | Year at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEYEAR will be ignored. Also, STATEYEAR, STATEMONTH, STATEDAY and STATESEC together specify the exact time point when the model state is to be saved. A typical use of this time point is the end of the whole simulation period. For example, if the last day of simulation is 2000-01-10 (VIC always runs until the last time step of this day in the case of sub-daily time step), then user should specify STATEYEAR=2000, STATEMONTH=1, STATEDAY=11, STATESEC=0 to save the state of the end of simulation. [See here for more detail](StateFile.md) |
| STATEMONTH   | integer | month         | Month at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEMONTH will be ignored.                                                                                                                                                                   |
| STATEDAY     | integer | day           | Day at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEDAY will be ignored.                                                                                                                                                                       |
| STATESEC     | integer | second        | Second at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATESEC will be ignored.                                                                                                                                                                    |
| STATE_FORMAT | string  | N/A           | Output state netCDF file format. Valid options: NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET, NETCDF4_CLASSIC, NETCDF4. *NOTE*: if STATENAME is not specified, STATE_FORMAT will be ignored.                                                                                                       |

# Define Meteorological and Vegetation Forcing Files

This section describes how to define the forcing files needed by the VIC model. VIC handles vegetation historical timeseries (LAI, albedo, vegetation canopy cover fraction) similarly to meteorological forcings (with some exceptions; see below).

In image driver, the meteorological forcings for the whole domain are stored in netCDF files, with the data in each year stored in separate netCDF files, with the same filepath and prefix and ended in "YYYY.nc".

VIC will allow forcing data to be stored in two different sets of files with different file prefix (e.g., precip and wind speed in one file, tmin and tmax in another file; or meteorological variables in one file, and vegetation timeseries in another file). Note that if you are using two sets of forcing files, the parameters for the first file must be defined before those for the second.

| Name          | Type          | Units                    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|---------------|---------------|--------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FORCING1      | string        | pathname and file prefix | First forcing netCDF file name, always required. This must precede all other forcing parameters used to define the first forcing file. Forcing netCDF files must be in the format of one file for one year of forcing data, with file name ended in "<YYYY>.nc" where <YYYY> is the year. Only need to include the directory and file name prefix in the variable FORCING1; "<YYYY>.nc" will be automatically appended during VIC run.                                                                      |
| FORCING2      | string        | pathname and file prefix | Second forcing file name, or FALSE if only one file used.This must precede all other forcing parameters used to define the second forcing file, and follow those used to define the first forcing file.                                                                                                                                                                                                                                                                                                     |
| FORCE_TYPE    | string string | N/A                      | Defines what forcing types are read from the file, followed by corresponding netCDF variable name (separated by space or tab). The required forcing types are: AIR_TEMP, PREC, PRESSURE, SWDOWN, LWDOWN, VP, WIND.                                                                                                                                                                                                                                                                                          |
| WIND_H        | float         | m                        | Height of wind speed measurement over bare soil and snow cover. Wind measurement height over vegetation is now read from the vegetation library file for all types, the value in the global file only controls the wind height over bare soil and over the snow pack when a vegetation canopy is not defined. *Note*: in image driver, this global parameter is only used in precipitation correction (if enabled); wind measurement height over bare soil is actually read from the parameter netCDF file. |
| CANOPY_LAYERS | int           | N/A                      | Number of canopy layers in the model. Default: 3.                                                                                                                                                                                                                                                                                                                                                                                                                                                           |

- If using one forcing file, use only FORCING1, if using two forcing files, define all parameters for FORCING1, and then define all forcing parameters for FORCING2\. All parameters need to be defined for both forcing files when a second file is used.

See the [example file](#example-global-parameter-file) at the end of this page for an example forcing parameter setup.


# Define Domain file

The folloiwng options describe the input domain file information. See [Domain File](Domain.md) for details about the domain file.

| Name        | Type          | Units    | Description                                                                                                                                  |
|-------------|---------------|----------|----------------------------------------------------------------------------------------------------------------------------------------------|
| DOMAIN      | string        | pathname | Domain netCDF file path.                                                                                                                     |
| DOMAIN_TYPE | string string | N/A      | Domain variable type, followed by corresponding netCDF variable name. Domain variable types include: LAT, LON, MASK, AREA, FRAC, YDIM, XDIM. |


# Define Parameter Files

The following options describe the input parameter files.

| Name               | Type   | Units         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------------|--------|---------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| PAREMETERS         | string | path/filename | Parameter netCDF file path, including soil parameters. vegetation library, vegetation parameters and snow band information (if SNOW_BAND=TRUE).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| BASEFLOW           | string | N/A           | This option describes the form of the baseflow parameters in the soil parameter file. Valid options: ARNO, NIJSSEN2001. See classic driver global parameter file for detail (../Classic/GlobalParam.md).                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| JULY_TAVG_SUPPLIED | string | TRUE or FALSE | If TRUE then VIC will expect an additional variable in the parameter file (July_Tavg) to contain the grid cell's average July temperature. *NOTE*: Supplying July average temperature is only required if the COMPUTE_TREELINE option is set to TRUE. <br><br>Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                            |
| ORGANIC_FRACT      | string | TRUE or FALSE | TRUE = the parameter file contains extra variables: the organic fraction, and the bulk density and soil particle density of the organic matter in each soil layer. FALSE = the parameter file does not contain any information about organic soil, and organic fraction should be assumed to be 0. <br><br>Default = FALSE.                                                                                                                                                                                                                                                                                                                                               |
| ALB_SRC         | string          | N/A                 | This option tells VIC where to look for ALBEDO values: if FROM_VEGLIB or FROM_VEGPARAM = Use the ALBEDO values from the parameter file. FROM_VEGPARAM = Use the ALBEDO values listed in the vegetation parameter file. If FROM_VEGHIST = Use the ALBEDO values from the veg_hist forcing files. Note: for this to work, ALBEDO must be supplied in the veg_hist files and listed in the global parameter file as one of the variables in the files. <br><br>Default = FROM_VEGLIB.  |
| LAI_SRC            | string | N/A           | This option tells VIC where to look for LAI values: FROM_VEGLIB = Use the LAI values listed in the vegetation library file. FROM_VEGPARAM = Use the LAI values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_LAI must be TRUE. FROM_VEGHIST = Use the LAI values listed in the veg_hist forcing files. Note: for this to work, LAI must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. Default = FROM_VEGLIB.                                                                                                                                                                |
| FCAN_SRC           | string | N/A           | This option tells VIC where to look for FCANOPY values: FROM_DEFAULT = Set FCANOPY to 1.0 for all veg classes, all times, and all locations. FROM_VEGLIB = Use the FCANOPY values listed in the vegetation library file. Note: for this to work, VEGLIB_FCANOPY must be TRUE.. FROM_VEGPARAM = Use the FCANOPY values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_FCANOPY must be TRUE. FROM_VEGHIST = Use the FCANOPY values listed in the veg_hist forcing files. Note: for this to work, FCANOPY must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. Default = FROM_DEFAULT. |
| SNOW_BAND          | string | TRUE or FALSE | Whether multiple snowband option is enables. If TRUE, snowband information is specified in the parameter netCDF file. <br><br>Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| CONSTANTS          | string | path/filename | Constants / Parameters file name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

# Lake Parameters

The following options only take effect when the lake model is running.

| Name          | Type      | Units                                         | Description                                                                                                                                                                                                                                               |
|-------------- |--------   |---------------------------------------------  |--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   |
| LAKES         | string    | FALSE or path/filename                        | Options for handling lakes: <li>**FALSE** = do not simulate lakes lake parameter path/filename = simulate lakes and read the given file for lake model parameters <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| LAKE_PROFILE  | string    | TRUE or FALSE                                 | Options for describing lake profile: <li>**FALSE** = VIC computes a parabolic depth-area profile for the lake basin <li>**TRUE** = VIC reads user-specified depth-area profile from the lake parameter file <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| EQUAL_AREA    | string    | TRUE or FALSE                                 | Options for computing grid cell areas: <li>**FALSE** = Grid cell boundaries and centers fall on a regular lat-lon grid, i.e. grid cells appear as squares when plotted in geographic projection; this means that the grid cells do not have equal areas. <li>**TRUE** = Grid cells have equal area, i.e. they appear as squares when plotted in an equal-area projection; this means that their boundaries do not fall on a regular lat-lon grid and the cell centers are not equally-spaced in latitude and longitude. <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| RESOLUTION    | float     | decimal degrees of latitude or area in km<sup>2</sup> | Options for grid cell resolution: <li>If **EQUAL_AREA = FALSE**: width of grid cells, in decimal degrees latitude or longitude. <li>If **EQUAL_AREA = TRUE**: area of grid cells, in km<sup>2</sup>. <br><br>Default = none; this MUST be set by the user to match the grid cell size if the lake model is running.                                                                                                                                                    |

# Define Output Files

The following options describe the location of the log and model history files. Click [here](OutputFormatting.md) for more information.

| Name                  | Type      | Units             | Description                                                                        |
|---------------------- |---------  |---------------    |----------------------------------------------------------------------------------- |
| LOG_DIR               | string    | path name         | Name of directory where log files should be written (optional, default is stdout)  |
| RESULT_DIR            | string    | path name         | Name of directory where model results are written                                  |

The following options describe the settings for each output stream:

| Name       | Type                                 | Units                                | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
|------------|--------------------------------------|--------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| OUTFILE*   | string                               | prefix                               | Information about this output file: Prefix of the output file (to which the first time step of the results contained in a file as well as suffix ".nc" will be appended). This should be specified once for each stream of output files.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| AGGFREQ    | string [integer/string]              | frequency count                      | Describes aggregation frequency for output stream. Valid options for frequency are: NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. Count may be an positive integer or a string with date format YYYY-MM-DD[-SSSSS] in the case of DATE. Default frequency is NDAYS. <bar><br>Default count is 1.                                                                                                                                                                                                                                                                                                                                                                                                    |
| HISTFREQ   | string [integer/string]              | frequency count                      | Describes the frequency/length of output results to be put in an individual file. Valid options are: NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. <br><br>Default is to output all results to one single file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| COMPRESS   | string/integer                       | TRUE, FALSE, or lvl                  | if TRUE or > 0 compress input and output files when done (uses gzip), if an integer [1-9] is supplied, it is used to set thegzip compression level                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| OUT_FORMAT | string                               | N/A                                  | Output netCDF format. Valid options:NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET, NETCDF4_CLASSIC, NETCDF4.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| OUTVAR*    | string string string integer string  | name format type multiplier aggtype  | Information about this output variable: <br>Name (must match a name listed in vic_driver_shared_all.h) <br>Output format (not used in image driver, replaced by "*") <br>Data type (one of: OUT_TYPE_DEFAULT, OUT_TYPE_CHAR, OUT_TYPE_SINT, OUT_TYPE_USINT, OUT_TYPE_INT, OUT_TYPE_FLOAT,OUT_TYPE_DOUBLE) <br>Multiplier - number to multiply the data with in order to recover the original values (only valid with OUT_FORMAT=BINARY) <br>Aggregation method - temporal aggregation method to use (one of: AGG_TYPE_DEFAULT, AGG_TYPE_AVG, AGG_TYPE_BEG, AGG_TYPE_END, AGG_TYPE_MAX, AGG_TYPE_MIN, AGG_TYPE_SUM) This should be specified once for each output variable. [Click here for more information](OutputFormatting.md). |

 - *Note: `OUTFILE`, and `OUTVAR` are optional; if omitted, traditional output files are produced. [Click here for details on using these instructions](OutputFormatting.md).*

## Example Global Parameter File:
```
#######################################################################
# VIC Model Parameters - 5.0
#######################################################################
# $Id$
#######################################################################
# Simulation Parameters
#######################################################################
MODEL_STEPS_PER_DAY   24  # number of model time steps in 24 hour period
SNOW_STEPS_PER_DAY    24  # number of snow model time steps in 24 hour period
RUNOFF_STEPS_PER_DAY  24  # number of runoff time steps in 24 hour period

STARTYEAR   1949 # year model simulation starts
STARTMONTH  1   # month model simulation starts
STARTDAY    1   # day model simulation starts
ENDYEAR     1949
ENDMONTH    1
ENDDAY      10
CALENDAR    PROLEPTIC_GREGORIAN

FULL_ENERGY FALSE   # calculate full energy balance
FROZEN_SOIL FALSE   # calculate frozen soils

#######################################################################
# DOMAIN INFO
#######################################################################
DOMAIN         domain.Stehekin.nc
DOMAIN_TYPE    LAT     lat
DOMAIN_TYPE    LON     lon
DOMAIN_TYPE    MASK    mask
DOMAIN_TYPE    AREA    area
DOMAIN_TYPE    FRAC    frac
DOMAIN_TYPE    YDIM    lat
DOMAIN_TYPE    XDIM    lon

#######################################################################
# Energy Balance Parameters
#######################################################################
FULL_ENERGY     FALSE   # TRUE = calculate full energy balance; FALSE = compute water balance only.  Default = FALSE.
#CLOSE_ENERGY   FALSE   # TRUE = all energy balance calculations (canopy air, canopy snow, ground snow,
                        # and ground surface) are iterated to minimize the total column error.  Default = FALSE.

#######################################################################
# Soil Temperature Parameters
# VIC will choose appropriate value for QUICK_FLUX depending on values of FULL_ENERGY and FROZEN_SOIL; the user should only need to override VIC's choices in special cases.
# The other options in this section are only applicable when FROZEN_SOIL is TRUE and their values depend on the application.
#######################################################################
FROZEN_SOIL FALSE   # TRUE = calculate frozen soils.  Default = FALSE.
#QUICK_FLUX FALSE   # TRUE = use simplified ground heat flux method of Liang et al (1999); FALSE = use finite element method of Cherkauer et al (1999)
#IMPLICIT   TRUE    # TRUE = use implicit solution for soil heat flux equation of Cherkauer et al (1999), otherwise uses original explicit solution.  Default = TRUE.
#QUICK_SOLVE    FALSE   # TRUE = Use Liang et al., 1999 formulation for iteration, but explicit finite difference method for final step.
#NO_FLUX        FALSE   # TRUE = use no flux lower boundary for ground heat flux computation; FALSE = use constant flux lower boundary condition.  If NO_FLUX = TRUE, QUICK_FLUX MUST = FALSE.  Default = FALSE.
#EXP_TRANS  TRUE    # TRUE = exponentially distributes the thermal nodes in the Cherkauer et al. (1999) finite difference algorithm, otherwise uses linear distribution.  Default = TRUE.
#GRND_FLUX_TYPE GF_410  # Options for ground flux:
#           # GF_406 = use (flawed) formulas for ground flux, deltaH, and fusion from VIC 4.0.6 and earlier;
#           # GF_410 = use formulas from VIC 4.1.0 (ground flux, deltaH, and fusion are correct; deltaH and fusion ignore surf_atten);
#           # Default = GF_410
#TFALLBACK  TRUE    # TRUE = when temperature iteration fails to converge, use previous time step's T value
#SPATIAL_FROST  FALSE   (Nfrost)    # TRUE = use a uniform distribution to simulate the spatial distribution of soil frost; FALSE = assume that the entire grid cell is frozen uniformly.  If TRUE, then replace (Nfrost) with the number of frost subareas, i.e., number of points on the spatial distribution curve to simulate.  Default = FALSE.

#######################################################################
# Precip (Rain and Snow) Parameters
# Generally these default values do not need to be overridden
#######################################################################
#SNOW_DENSITY   DENS_BRAS   # DENS_BRAS = use traditional VIC algorithm taken from Bras, 1990; DENS_SNTHRM = use algorithm taken from SNTHRM model.
#BLOWING        FALSE   # TRUE = compute evaporative fluxes due to blowing snow
#COMPUTE_TREELINE   FALSE   # Can be either FALSE or the id number of an understory veg class; FALSE = turn treeline computation off; VEG_CLASS_ID = replace any overstory veg types with the this understory veg type in all snow bands for which the average July Temperature <= 10 C (e.g. "COMPUTE_TREELINE 10" replaces any overstory veg cover with class 10)
#CORRPREC   FALSE   # TRUE = correct precipitation for gauge undercatch
#SPATIAL_SNOW   FALSE   # TRUE = use a uniform distribution to simulate the partial coverage of the
                        # surface by a thin snowpack.  Coverage is assumed to be uniform after snowfall
                        # until the pack begins to melt.  If TRUE, VIC will expect an additional column
                        # in the soil paramter file containing the snow distibution slope parameter
                        # (= 2 * snow depth below which coverage < 1).

#######################################################################
# Turbulent Flux Parameters
# Generally these default values do not need to be overridden
#######################################################################
#MIN_WIND_SPEED 0.1 # minimum allowable wind speed (m/s)
#AERO_RESIST_CANSNOW    AR_406_FULL # Options for aerodynamic resistance in snow-filled canopy:
#           # AR_406    = multiply by 10 for latent heat but do NOT multiply by 10 for sensible heat and do NOT apply stability correction (as in VIC 4.0.6); when no snow in canopy, use surface aero_resist for ET.
#           # AR_406_LS     = multiply by 10 for latent heat AND sensible heat and do NOT apply stability correction; when no snow in canopy, use surface aero_resist for ET.
#           # AR_406_FULL   = multiply by 10 for latent heat AND sensible heat and do NOT apply stability correction; additionally, always use overstory aero_resist for ET (as in 4.1.0).
#           # AR_410    = apply stability correction but do NOT multiply by 10 (as in VIC 4.1.0); additionally, always use overstory aero_resist for ET (as in 4.1.0).
#           # Default   = AR_406_FULL

#######################################################################
# Meteorological Forcing Disaggregation Parameters
# Generally these default values do not need to be overridden
#######################################################################
#PLAPSE     TRUE    # This controls how VIC computes air pressure when air pressure is not supplied as an input forcing: TRUE = set air pressure to sea level pressure, lapsed to grid cell average elevation; FALSE = set air pressure to constant 95.5 kPa (as in all versions of VIC pre-4.1.1)
#SW_PREC_THRESH     0   # Minimum daily precip [mm] that can cause dimming of incoming shortwave; default = 0.
#MTCLIM_SWE_CORR    TRUE    # This controls VIC's estimates of incoming shortwave in the presence of snow; TRUE = adjust incoming shortwave for snow albedo effect; FALSE = do not adjust shortwave; default = TRUE
#VP_ITER        VP_ITER_ANNUAL  # This controls VIC's iteration between estimates of shortwave and vapor pressure:
#           # VP_ITER_NEVER = never iterate; make estimates separately
#           # VP_ITER_ALWAYS = always iterate once
#           # VP_ITER_ANNUAL = iterate once for arid climates based on annual Precip/PET ratio
#           # VP_ITER_CONVERGE = iterate until shortwave and vp stabilize
#           # default = VP_ITER_ALWAYS
#VP_INTERP  TRUE    # This controls sub-daily humidity estimates; TRUE = interpolate daily VP estimates linearly between sunrise of one day to the next; FALSE = hold VP constant for entire day
#LW_TYPE        LW_PRATA    # This controls the algorithm used to estimate clear-sky longwave radiation:
#           # LW_TVA = Tennessee Valley Authority algorithm (1972) (this was traditional VIC algorithm)
#           # other options listed in vic_driver_shared_all.h
#           # default = LW_PRATA
#LW_CLOUD   LW_CLOUD_DEARDORFF  # This controls the algorithm used to estimate the influence of clouds on total longwave:
#           # LW_CLOUD_BRAS = method from Bras textbook (this was the traditional VIC algorithm)
#           # LW_CLOUD_DEARDORFF = method of Deardorff (1978)
#           # default = LW_CLOUD_DEARDORFF

#######################################################################
# Carbon Cycle Parameters
#######################################################################
#CARBON         FALSE       # TRUE = simulate carbon cycle; FALSE = do not simulate carbon cycle.  Default = FALSE.
#VEGLIB_PHOTO   FALSE       # TRUE = photosynthesis parameters are included in the veg library file.  Default = FALSE.
#RC_MODE    RC_JARVIS   # RC_JARVIS = canopy resistance computed by applying resistance factors to the veg class's minimum resistance, listed in the veg library
                            # RC_PHOTO = canopy resistance computed by applying resistance factors to the minimum resistance required by current photosynthetic demand.  Default = RC_JARVIS.

#######################################################################
# Miscellaneous Simulation Parameters
# Generally these default values do not need to be overridden
#######################################################################
#CONTINUEONERROR    TRUE    # TRUE = if simulation aborts on one grid cell, continue to next grid cell

#######################################################################
# State Files and Parameters
#######################################################################
#INIT_STATE  # Initial state path/file
#STATENAME output/image/Stehekin/states  # Output state file path/prefix. The time (STATEYEAR,STATEMONTH,STATEDAY,STATESEC) will be appended to the prefix automatically in the format yyyymmdd.
#STATEYEAR   1949    # year to save model state
#STATEMONTH  1  # month to save model state
#STATEDAY    10  # day to save model state
#STATESEC    82800  # second to save model state
#STATE_FORMAT           NETCDF4_CLASSIC  # State file format, valid options:
#NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET, NETCDF4_CLASSIC, NETCDF4

#######################################################################
# Forcing Files and Parameters
# netcdf forcing files will be of the form: <FORCING1>YYYY.nc
#######################################################################
FORCING1      forcings/Stehekin_image_test.forcings_10days.
FORCE_TYPE    AIR_TEMP -   tas    # Average air temperature, K
FORCE_TYPE    PREC -       prcp   # Total precipitation (rain and snow), kg/m2/s
FORCE_TYPE    PRESSURE -   pres   # Atmospheric pressure, Pa
FORCE_TYPE    SWDOWN       dswrf  # Incoming shortwave, W/m2
FORCE_TYPE    LWDOWN ---     dlwrf  # Incoming longwave radiation, W/m2
FORCE_TYPE    VP           vp   # Vapor pressure, kPa
FORCE_TYPE    WIND         wind   # Wind speed, m/s
# WIND_H        10.0                # height of wind speed measurement. NOTE: in image driver, this global parameter is only used for precipitation correction (if enabled); wind measurement height over bare soil is read from the parameter netCDF file.

#######################################################################
# Land Surface Files and Parameters
#######################################################################
PARAMETERS      params/Stehekin.params.nc
SNOW_BAND       TRUE
BASEFLOW        ARNO
JULY_TAVG_SUPPLIED  FALSE
ORGANIC_FRACT       FALSE
LAI_SRC             FROM_VEGPARAM
NODES           3  # number of soil thermal nodes

#######################################################################
# Lake Simulation Parameters
# These need to be un-commented and set to correct values only when running lake model (LAKES is not FALSE)
#######################################################################
#LAKES      (put lake parameter path/file here) # Lake parameter path/file
#LAKE_PROFILE   FALSE   # TRUE = User-specified depth-area parameters in lake parameter file; FALSE = VIC computes a parabolic depth-area profile
#EQUAL_AREA FALSE   # TRUE = grid cells are from an equal-area projection; FALSE = grid cells are on a regular lat-lon grid
#RESOLUTION 0.125   # Grid cell resolution (degrees if EQUAL_AREA is FALSE, km^2 if EQUAL_AREA is TRUE); ignored if LAKES is FALSE

#######################################################################
# Output Files and Parameters
#######################################################################
LOG_DIR         (put the log directory path here)       # Log directory path
RESULT_DIR      (put the result directory path here)    # Results directory path

#######################################################################
#
# Output File Contents
#
# You can specify your output file names and contents in the global param file
# (see the VIC documentation for more information).
#
# If you do not specify file names and contents in the global param
# file, VIC will produce the same set of output files that it has
# produced in earlier versions, namely "fluxes" and "snow" files, plus
# "fdepth" files if FROZEN_SOIL is TRUE and "snowband" files if
# snowbands are specified.  These files will have the same contents and
# format as in earlier versions.
#
# Format:
#
# Output File Contents
# OUTFILE _prefix_
# OUTFREQ         _freq_          _VALUE_
# HISTFREQ        _freq_          _VALUE_
# COMPRESS        _compress_
# OUT_FORMAT      _nc_format_
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
#
# OUTFILE _prefix_
# OUTFREQ         _freq_          _VALUE_
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
# OUTVAR  _varname_   [_format_  [_type_ [_multiplier_ [_aggtype_]]]]
#
#
# _prefix_     = name of the output file, NOT including the date stamp or the suffix
# _freq_       = Describes aggregation frequency for output stream. Valid
#                options for frequency are:
#                  NEVER     = never write to history file
#                  NSTEPS    = write to history every _value_ steps
#                  NSECONDS  = write to history every _value_ seconds
#                  NMINUTES  = write to history every _value_ minutes
#                  NHOURS    = write to history every _value_ hours
#                  NDAYS     = write to history every _value_ days
#                  NMONTHS   = write to history every _value_ months
#                  NYEARS    = write to history every _value_ years
#                  DATE      = write to history on the date: _value_
#                  END       = write to history at the end of the simulation
# _value_      = integer describing the number of _freq_ intervals to pass
#                before writing to the history file.
# _compress_   = netCDF gzip compression option.  TRUE, FALSE, or integer between 1-9.
# _nc_format_  = netCDF format. NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET,
#                NETCDF4_CLASSIC, or NETCDF4
# _varname_    = name of the variable (this must be one of the
#                output variable names listed in vic_driver_shared_all.h.)
#
# _format_     = not used in image driver, replace with *
#
# _type_, and _multiplier_, and _aggtype_ are optional.
# If these are omitted, the default values will be used.
#
# _type_       = data type code. Must be one of:
#                  OUT_TYPE_DOUBLE = double-precision floating point
#                  OUT_TYPE_FLOAT  = single-precision floating point
#                  OUT_TYPE_INT    = integer
#                  OUT_TYPE_USINT  = unsigned short integer
#                  OUT_TYPE_SINT   = short integer
#                  OUT_TYPE_CHAR   = char
#                  *               = use the default type
# _multiplier_ = (for binary output files) factor to multiply
#                the data by before writing, to increase precision.
#                  *    = use the default multiplier for this variable
# _aggtype_    = Aggregation method to use for temporal aggregation. Valid
#                options for aggtype are:
#                  AGG_TYPE_DEFAULT = default aggregation type for variable
#                  AGG_TYPE_AVG     = average over aggregation window
#                  AGG_TYPE_BEG     = beginning of aggregation window
#                  AGG_TYPE_END     = end of aggregation window
#                  AGG_TYPE_MAX     = maximum in aggregation window
#                  AGG_TYPE_MIN     = minimum in aggregation window
#                  AGG_TYPE_SUM     = sum over aggregation window
#
#
#######################################################################
```
