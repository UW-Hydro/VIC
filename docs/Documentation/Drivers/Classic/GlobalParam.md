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
| STARTYEAR            | integer | year   | Year   model simulation starts. **NOTE**: STARTYEAR, STARTMONTH, STARTDAY and STARTSEC together specify the begenning time point of the first simulation time step.  |
| STARTMONTH           | integer | month  | Month   model simulation starts                                                                                                                                                                  |
| STARTDAY             | integer | day    | Day model   simulation starts                                                                                                                                                                    |
| STARTSEC             | integer | second | Second   model simulation starts (e.g., STARTSEC=0 for starting from the beginning of   a day; STARTSEC=3600 for starting from the beginning of the second hour of a   day) <br><br>Default = 0. |
| **EITHER:**          |         |        | *Note:*   **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but   **not both**                                                                                                |
| NRECS                | integer | N/A    | Number of   time steps over which to run model. ***NOTE: The number of records must be   defined such that the model completes the last day.***                                                  |
| **OR:**              |         |        | *Note:*   **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but   **not both**                                                                                                |
| ENDYEAR              | integer | year   | Year   model simulation ends. **NOTE**: ENDYEAR, ENDMONTH and ENDDAY together specify the last day on which the simulation runs. For sub-daily model run timestep, VIC always runs until the last timestep of this day. |
| ENDMONTH             | integer | month  | Month   model simulation ends                                                                                                                                                                    |
| ENDDAY               | integer | day    | Day model   simulation ends                                                                                                                                                                      |
| CALENDAR             | string  | N/A    | Calendar   to use. Valid calendars: 'standard', 'gregorian', 'proleptic_gregorian'   'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day'.                                            |
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
| SHARE_LAYER_MOIST | string            | TRUE or FALSE                      | If TRUE, then *if* the soil moisture in the layer that contains more than half of the roots is above the critical point, then the plant's roots in the drier layers can access the moisture of the wetter layer so that the plant does not experience moisture limitation. <br> If FALSE or all of the soil layer moistures are below the critical point, transpiration in each layer is limited by the layer's soil moisture. <br><br> Default: TRUE.              |
| SPATIAL_FROST     | string (+integer) | string: TRUE or FALSE integer: N/A | Option to allow spatial heterogeneity in soil temperature:FALSE = Assume soil temperature is horizontally constant (only varies with depth). TRUE = Assume soil temperatures at each given depth are distributed horizontally with a uniform (linear) distribution, so that even when the mean temperature is below freezing, some portion of the soil within the grid cell at that depth could potentially be above freezing. This requires specifying a frost slope value as an extra field in the soil parameter file, so that the minimum/maximum temperatures can be computed from the mean value. The maximum and minimum temperatures will be set to mean temperature +/- frost_slope.If TRUE is specified, you must follow this with an integer value for Nfrost, the number of frost sub-areas (each having a distinct temperature). Default = FALSE.                                                                                                                                                                                                                                       |

## Precipitation (Rain and Snow) Parameters

Generally these default values do not need to be overridden.

| Name                  | Type              | Units                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|-----------------------|-------------------|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| SNOW_DENSITY          | string            | N/A                   | Options for computing snow density:DENS_BRAS = Use traditional VIC algorithm taken from Bras, 1990.DENS_SNTHRM = Use algorithm taken from SNTHRM model.Default = DENS_BRAS.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| BLOWING               | string            | TRUE or FALSE         | If TRUE, compute evaporative fluxes due to blowing snow. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| BLOWING_VAR_THRESHOLD | string            | TRUE or FALSE    |  If TRUE, a variable shear stress threshold is used to determine the blowing snow flux. If FALSE, then a fixed threshold is used. See Li and Pomeroy (1997) for details. <br><br>Default: TRUE. |
| BLOWING_CALC_PROB     | string            | TRUE or FALSE   | If TRUE, the probability of occurrence of blowing snow is calculated as a function of environmental conditions. If FALSE, then the probability is set to 1. See Lu and Pomeroy (1997) for details. <br><br>Default: TRUE.  |
| BLOWING_SIMPLE        | string            | TRUE or FALSE  | If TRUE, the sublimation flux of blowing snow is calculated as a function vapor pressure and wind speed. If FALSE, then additional calculations are made to account for a saltation and suspension layer. See Lu and Pomeroy (1997) for details. <br><br>Default: FALSE. |
| BLOWING_FETCH         | string            | TRUE or FALSE   | This option is only used when BLOWING_SIMPLE is set to FALSE. When this option is set to TRUE, the fetch is accounted for in the calculation of the sublimation flux from blowing snow. If FALSE then the fetch is not used. See Lu and Pomeroy (1997) for details. <br><br> Default: TRUE. |
| BLOWING_SPATIAL_WIND  | string            | TRUE or FALSE  | If TRUE, multiple wind speed ranges, calculated according to a probability distribution, are used to determine the sublimation flux from blowing snow. If FALSE, then a single wind speed is used. See Lu and Pomeroy (1997) for details. <br><br>Default: TRUE. |
| COMPUTE_TREELINE      | string or integer | FALSE or veg class id | Options for handling above-treeline vegetation:FALSE = Do not compute treeline or replace vegetation above the treeline.CLASS_ID = Compute the treeline elevation based on average July temperatures; for those elevation bands with elevations above the treeline (or the entire grid cell if SNOW_BAND == 1 and the grid cell elevation is above the tree line), if they contain vegetation tiles having overstory, replace that vegetation with the vegetation having id CLASS_ID in the vegetation library. NOTE 1: You MUST supply VIC with a July average air temperature, in the optional July_Tavg field, AND set theJULY_TAVG_SUPPLIED option to TRUE so that VIC can read the soil parameter file correctly. NOTE 2: If LAKES=TRUE, COMPUTE_TREELINE MUST be FALSE.Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| CORRPREC              | string            | TRUE or FALSE         | If TRUE correct precipitation for gauge undercatch. NOTE: This option is not supported when using snow/elevation bands. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
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

| Name         | Type    | Units           | Description                                                                                                                                                                                                                                                                                 |
|--------------|---------|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| INIT_STATE   | string  | path/filename   | Full path and filename of initial state file. *NOTE*: if INIT_STATE is not specified, VIC will take initial soil moistures from the soil parameter file and set all other state variables to a default state.                                                                               |
| STATENAME    | string  | path/filename   | Path and file prefix of the state file to be created on the specified date. The date and second within the simulation at which the state is saved will be appended to the file prefix to form a complete file name. *NOTE*: if STATENAME is not specified, VIC will not save its state in a statefile. |
| STATEYEAR    | integer | year            | Year at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEYEAR will be ignored. Also, STATEYEAR, STATEMONTH, STATEDAY and STATESEC together specify the exact time point when the model state is to be saved. A typical use of this time point is the end of the whole simulation period. For example, if the last day of simulation is 2000-01-10 (VIC always runs until the last time step of this day in the case of sub-daily time step), then user should specify STATEYEAR=2000, STATEMONTH=1, STATEDAY=11, STATESEC=0 to save the state of the end of simulation. [See here for more detail](StateFile.md) |
| STATEMONTH   | integer | month           | Month at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEMONTH will be ignored.                                                                                                                                                                   |
| STATEDAY     | integer | day             | Day at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATEDAY will be ignored.                                                                                                                                                                       |
| STATESEC     | integer | second          | Second at which model simulation state should be saved. *NOTE*: if STATENAME is not specified, STATESEC will be ignored.                                                                                                                                                                    |
| STATE_FORMAT | string  | BINARY OR ASCII | If ASCII, VIC reads/writes the intial/output state files in ASCII format. If BINARY, VIC reads/writes intial/output state files in binary format. NOTE: if INIT_STATE or STATENAME are not specified, STATE_FORMAT will be ignored.                                                         |

# Define Meteorological and Vegetation Forcing Files

This section describes how to define the forcing files needed by the VIC model.  VIC handles vegetation historical timeseries (LAI, albedo, vegetation canopy cover fraction) similarly to meteorological forcings (with some exceptions; see below).

Unlike model parameters, for which 1 file contains data for all grid cells, the meteorological forcings are stored as a separate time series for each grid cell. The time step length of the input forcings must match the time step length at which VIC is running. Input files can be ASCII or Binary (signed or unsigned short ints) column formatted. Columns in the file must be in the same order as they are defined in the global control file.

VIC will allow forcing data to be stored in two different files per grid cell (e.g., precip and wind speed in one file, tmin and tmax in another file; or meteorological variables in one file, and vegetation timeseries in another file). Note that if you are using two forcing files per grid cell, the parameters for the first file must be defined before those for the second. **Bold** numbers indicate the order in which these values should be defined, after each forcing file (`FORCING1` or `FORCING2`). Options that do not have a bold number apply to both forcing file types and should appear after the numbered options.

All FORCING filenames are actually the pathname, and prefix for gridded data types: ex. `DATA/forcing_YY.YYY_XX.XXX`. Latitude and longitude index suffix is added by VIC based on the `GRID_DECIMAL` parameter defined above, and the latitude and longitude values defined in the [soil parameter file](SoilParam.md).

| Name                | Type              | Units                       | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
|---------------------|-------------------|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| FORCING1            | string            | pathname and file prefix    | First forcing file name, always required. This must precede all other forcing parameters used to define the first forcing file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| FORCING2            | string            | pathname and file prefix    | Second forcing file name, or FALSE if only one file used.This must precede all other forcing parameters used to define the second forcing file, and follow those used to define the first forcing file.                                                                                                                                                                                                                                                                                                                                                                                                                        |
| FORCE_FORMAT        | string            | BINARY or ASCII             | Defines the format type for the forcing files.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| FORCE_ENDIAN        | string            | BIG or LITTLE               | Identifies the architecture of the machine on which the binary forcing files were created:BIG = big-endian (e.g. SUN).LITTLE = little-endian (e.g. PC/linux). Model will identify the endian of the current machine, and swap bytes if necessary. Required for binary forcing file, not used for ASCII forcing file.                                                                                                                                                                                                                                                                                                           |
| N_TYPES             | int               | N/A                         | Number of columns in the current data file, with the following exception: for the vegetation history variables ALBEDO, LAI, and FCANOPY, there must be multiple columns for these variables, one per vegetation tile. In this case, ALBEDO, LAI, and FCANOPY each count as only 1 variable despite covering multiple columns.                                                                                                                                                                                                                                                                                            |
| FORCE_TYPE          | stringstringfloat | VarName(un)signedmultiplier | Defines what forcing types are read from the file, and in what order. For ASCII file only the forcing type needs to be defined, but for Binary file each line must also define whether the column is SIGNED or UNSIGNED short int and by what factor values are multiplied before being written to output. Note: Unlike other variables, ALBEDO, LAI, and FCANOPY, each span multiple columns, one column per veg tile. This will generally vary from one grid cell to the next as the number of veg tiles varies. However, ALBEDO, LAI, and FCANOPY should each have only one FORCE_TYPE entry. [Click here for details](ForcingData.md). |
| FORCE_STEPS_PER_DAY | integer           | steps                       | Number of timesteps per day in forcing file (must be = SNOW_STPES_PER_DAY)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| FORCEYEAR           | integer           | year                        | Year meteorological forcing files start                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| FORCEMONTH          | integer           | month                       | Month meteorological forcing files start                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| FORCEDAY            | integer           | day                         | Day meteorological forcing files start                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| FORCESEC            | integer           | second                      | Second meteorological forcing files start. <br><br> Default: 0.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| GRID_DECIMAL        | integer           | N/A                         | Number of decimals to use in gridded file name extensions                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| WIND_H              | float             | m                           | Height of wind speed measurement over bare soil and snow cover. Wind measurement height over vegetation is now read from the vegetation library file for all types, the value in the global file only controls the wind height over bare soil and over the snow pack when a vegetation canopy is not defined.                                                                                                                                                                                                                                                                                                                  |
| CANOPY_LAYERS       | int               | N/A                         | Number of canopy layers in the model. Default: 3.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |

- If using one forcing file, use only FORCING1, if using two forcing files, define all parameters for FORCING1, and then define all forcing parameters for FORCING2\. All parameters need to be defined for both forcing files when a second file is used.

_Examples._ a standard four column daily forcing data file will be defined as:

## ASCII File

    FORCING1             forcings/full_data_
    FORCE_FORMAT         ASCII
    FORCE_TYPE           PREC
    FORCE_TYPE           AIR_TEMP
    FORCE_TYPE           SWDOWN
    FORCE_TYPE           LWDOWN
    FORCE_TYPE           SKIP  # This column is air density, which is not needed by VIC
    FORCE_TYPE           PRESSURE
    FORCE_TYPE           VP
    FORCE_TYPE           WIND
    FORCE_STEPS_PER_DAY  24    # Forcing time step length (hours)
    FORCEYEAR            1949  # Year of first forcing record
    FORCEMONTH           01    # Month of first forcing record
    FORCEDAY             01    # Day of first forcing record
    GRID_DECIMAL         4     # Number of digits after decimal point in forcing file names
    WIND_H               10.0  # height of wind speed measurement (m)

## Binary File

    FORCING1             forcings/full_data_
    FORCE_FORMAT         BINARY
    FORCE_ENDIAN        LITTLE
    FORCE_TYPE           PREC       UNSIGNED    40
    FORCE_TYPE           AIR_TEMP   SIGNED      100
    FORCE_TYPE           SWDOWN     UNSIGNED    100
    FORCE_TYPE           LWDOWN     UNSIGNED    100
    FORCE_TYPE           PRESSURE   UNSIGNED    100
    FORCE_TYPE           VP         UNSIGNED    100
    FORCE_TYPE           WIND       UNSIGNED    100
    FORCE_STEPS_PER_DAY  24    # Forcing time step length (hours)
    FORCEYEAR            1949  # Year of first forcing record
    FORCEMONTH           01    # Month of first forcing record
    FORCEDAY             01    # Day of first forcing record
    GRID_DECIMAL         4     # Number of digits after decimal point in forcing file names
    WIND_H               10.0  # height of wind speed measurement (m)


# Define Parameter Files

The following options describe the input parameter files.

| Name               | Type            | Units               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------------|-----------------|---------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| SOIL               | string          | path/filename       | the Soil parameter file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| BASEFLOW           | string          | N/A                 | This option describes the form of the baseflow parameters in the soil parameter file:ARNO = fields 5-8 of the soil parameter file are the standard VIC baseflow parametersNIJSSEN2001 = fields 5-8 of the soil parameter file are the baseflow parameters from Nijssen et al (2001) Default = ARNO.                                                                                                                                                                                                                                                                                                                                                                       |
| JULY_TAVG_SUPPLIED | string          | TRUE or FALSE       | If TRUE then VIC will expect an additional column (July_Tavg) in the soil parameter file to contain the grid cell's average July temperature. If your soil parameter file contains this optional column, you MUST set JULY_TAVG_SUPPLIED to TRUE so that VIC can read the soil parameter file correctly. NOTE: Supplying July average temperature is only required if the COMPUTE_TREELINE option is set to TRUE. Default = FALSE.                                                                                                                                                                                                                                        |
| ORGANIC_FRACT      | string          | TRUE or FALSE       | TRUE = the soil parameter file contains 3*Nlayer extra columns, listing, for each layer: the organic fraction, and the bulk density and soil particle density of the organic matter in the soil layer. FALSE = the soil parameter file does not contain any information about organic soil, and organic fraction should be assumed to be 0. Default = FALSE.                                                                                                                                                                                                                                                                                                               |
| VEGLIB             | string          | path/filename       | Vegetation library file name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| VEGPARAM           | string          | path/filename       | Vegetation parameter file name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| ROOT_ZONES         | integer         | N/A                 | Number of defined root zones defined for root distribution.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| VEGPARAM_ALB       | string          | TRUE or FALSE       | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly ALBEDO values for each vegetation type for each grid cell. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ALB_SRC         | string          | N/A                 | This option tells VIC where to look for ALBEDO values: FROM_VEGLIB = Use the ALBEDO values listed in the vegetation library file. FROM_VEGPARAM = Use the ALBEDO values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_ALB must be TRUE. FROM_VEGHIST = Use the ALBEDO values listed in the veg_hist forcing files. Note: for this to work, ALBEDO must be supplied in the veg_hist files and listed in the global parameter file as one of the variables in the files. Default = FROM_VEGLIB.  |
| VEGPARAM_LAI       | string          | TRUE or FALSE       | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly LAI values for each vegetation type for each grid cell. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| LAI_SRC            | string          | N/A                 | This option tells VIC where to look for LAI values: FROM_VEGLIB = Use the LAI values listed in the vegetation library file. FROM_VEGPARAM = Use the LAI values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_LAI must be TRUE. FROM_VEGHIST = Use the LAI values listed in the veg_hist forcing files. Note: for this to work, LAI must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. Default = FROM_VEGLIB.                                                                                                                                                                 |
| VEGLIB_FCAN        | string          | TRUE or FALSE       | If TRUE the vegetation library file contains monthly FCANOPY values for each vegetation type for each grid cell (between the LAI and ALBEDO values). Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| VEGPARAM_FCAN      | string          | TRUE or FALSE       | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly FCANOPY values for each vegetation type for each grid cell. Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| FCAN_SRC           | string          | N/A                 | This option tells VIC where to look for FCANOPY values: FROM_DEFAULT = Set FCANOPY to 1.0 for all veg classes, all times, and all locations. FROM_VEGLIB = Use the FCANOPY values listed in the vegetation library file. Note: for this to work, VEGLIB_FCANOPY must be TRUE.. FROM_VEGPARAM = Use the FCANOPY values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_FCANOPY must be TRUE. FROM_VEGHIST = Use the FCANOPY values listed in the veg_hist forcing files. Note: for this to work, FCANOPY must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. Default = FROM_DEFAULT. |
| SNOW_BAND          | integer[string] | N/A [path/filename] | Maximum number of snow elevation bands to use, and the name (with path) of the snow elevation band file. For example: SNOW_BAND 5 path/filename. To turn off this feature, set the number of snow bands to 1 and do not follow this with a snow elevation band file name. Default = 1.                                                                                                                                                                                                                                                                                                                                                                                    |
| CONSTANTS          | string          | path/filename       | Constants / Parameters file name                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

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

| Name        | Type      | Units             | Description                                                                        |
|------------ |---------  |---------------    |----------------------------------------------------------------------------------- |
| OUTFILE\*   | string    | prefix            | Information about this output file: <br>Prefix of the output file (to which the lat and lon will be appended3) <br> This should be specified once for each output file. [Click here for more information.](OutputFormatting.md) |
| AGGFREQ     | string <br> [integer/string]   | frequency <br> count | Describes aggregation frequency for output stream.  Valid options for frequency are: NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. Count may be an positive integer or a string with date format YYYY-MM-DD[-SSSSS] in the case of DATE. <br> Default `frequency` is `NDAYS`. Default `count` is 1. |
| COMPRESS    | string/integer | TRUE, FALSE, or lvl | if TRUE or > 0 compress input and output files when done (uses `gzip`), if an integer [1-9] is supplied, it is used to set the`gzip` compression level |
| OUT_FORMAT  | string    | BINARY OR ASCII   | If BINARY write output files in binary (default is ASCII).                                                                                                                                  |
| OUTVAR\*    | <br> string <br> string <br> string <br> integer <br> string <br> | <br> name <br> format <br> type <br> multiplier <br> aggtype <br> | Information about this output variable:<br>Name (must match a name listed in vic_driver_shared_all.h) <br> Output format (C fprintf-style format code) (only valid with OUT_FORMAT=ASCII) <br>Data type (one of: OUT_TYPE_DEFAULT, OUT_TYPE_CHAR, OUT_TYPE_SINT, OUT_TYPE_USINT, OUT_TYPE_INT, OUT_TYPE_FLOAT,OUT_TYPE_DOUBLE) <br> Multiplier - number to multiply the data with in order to recover the original values (only valid with OUT_FORMAT=BINARY) <br> Aggregation method - temporal aggregation method to use (one of: AGG_TYPE_DEFAULT, AGG_TYPE_AVG, AGG_TYPE_BEG, AGG_TYPE_END, AGG_TYPE_MAX, AGG_TYPE_MIN, AGG_TYPE_SUM) <br> <br> This should be specified once for each output variable. [Click here for more information.](OutputFormatting.md)|

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
NLAYER      3   # number of soil layers
NODES       10  # number of soil thermal nodes
MODEL_STEPS_PER_DAY  24   # number of model time steps per day (set to 1 if FULL_ENERGY = FALSE, set to > 4 if FULL_ENERGY = TRUE)
SNOW_STEPS_PER_DAY  24   # number of time steps per day for which to solve the snow model (should = MODEL_STEPS_PER_DAY if MODEL_STEPS_PER_DAY > 1)
RUNOFF_STEPS_PER_DAY  24   # number of time steps per day for which to solve the runoff model (should be >= MODEL_STEPS_PER_DAY)
STARTYEAR   2000    # year model simulation starts
STARTMONTH  01  # month model simulation starts
STARTDAY    01  # day model simulation starts
ENDYEAR     2000    # year model simulation ends
ENDMONTH    12  # month model simulation ends
ENDDAY      31  # day model simulation ends

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
#MAX_SNOW_TEMP  0.5 # maximum temperature (C) at which snow can fall
#MIN_RAIN_TEMP  -0.5    # minimum temperature (C) at which rain can fall
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
#INIT_STATE (put the initial state path/filename here)  # Initial state path/file
#STATENAME  (put the path/prefix of output state file here) # Output state file path/prefix.  The date (STATEYEAR,STATEMONTH,STATEDAY) will be appended to the prefix automatically in the format yyyymmdd.
#STATEYEAR  2000    # year to save model state
#STATEMONTH 12  # month to save model state
#STATEDAY   31  # day to save model state
#STATE_FORMAT   ASCII  # BINARY OR ASCII

#######################################################################
# Forcing Files and Parameters
#
#       All FORCING filenames are actually the pathname, and prefix
#               for gridded data types: ex. DATA/forcing_
#               Latitude and longitude index suffix is added by VIC
#
#   There must be 1 FORCE_TYPE entry for each variable (column) in the forcing file
#
#   If FORCE_TYPE is BINARY, each FORCE_TYPE must be followed by:
#           SIGNED/UNSIGNED SCALE_FACTOR
#       For example (BINARY):
#           FORCE_TYPE  PREC    UNSIGNED    40
#       or (ASCII):
#           FORCE_TYPE  PREC
#######################################################################
FORCING1             forcings/full_data_
FORCE_FORMAT         ASCII
FORCE_TYPE           PREC
FORCE_TYPE           AIR_TEMP
FORCE_TYPE           SWDOWN
FORCE_TYPE           LWDOWN
FORCE_TYPE           SKIP  # This column is air density, which is not needed by VIC
FORCE_TYPE           PRESSURE
FORCE_TYPE           VP
FORCE_TYPE           WIND
FORCE_STEPS_PER_DAY  24    # Forcing time step length (hours)
FORCEYEAR            1949  # Year of first forcing record
FORCEMONTH           01    # Month of first forcing record
FORCEDAY             01    # Day of first forcing record
GRID_DECIMAL         4     # Number of digits after decimal point in forcing file names
WIND_H               10.0  # height of wind speed measurement (m)

#######################################################################
# Land Surface Files and Parameters
#######################################################################
SOIL            (put the soil parameter path/file here) # Soil parameter path/file
BASEFLOW    ARNO    # ARNO = columns 5-8 are the standard VIC baseflow parameters; NIJSSEN2001 = columns 5-8 of soil file are baseflow parameters from Nijssen et al (2001)
JULY_TAVG_SUPPLIED  FALSE   # TRUE = final column of the soil parameter file will contain average July air temperature, for computing treeline; this will be ignored if COMPUTE_TREELINE is FALSE; FALSE = compute the treeline based on the average July air temperature of the forcings over the simulation period
ORGANIC_FRACT   FALSE   # TRUE = simulate organic soils; soil param file contains 3*Nlayer extra columns, listing for each layer the organic fraction, and the bulk density and soil particle density of the organic matter in the soil layer; FALSE = soil param file does not contain any information about organic soil, and organic fraction should be assumed to be 0
VEGLIB          (put the veg library path/file here)    # Veg library path/file
VEGPARAM        (put the veg parameter path/file here)  # Veg parameter path/file
ROOT_ZONES      3   # Number of root zones (must match format of veg param file)
#VEGLIB_FCAN    FALSE   # TRUE = veg lib file contains 12 monthly values of partial vegcover fraction for each veg class, between the LAI and albedo values
#VEGPARAM_LAI   TRUE    # TRUE = veg param file contains LAI information; FALSE = veg param file does NOT contain LAI information
#VEGPARAM_ALB   FALSE    # TRUE = veg param file contains albedo information; FALSE = veg param file does NOT contain albedo information
#VEGPARAM_FCAN  FALSE    # TRUE = veg param file contains veg_cover information; FALSE = veg param file does NOT contain veg_cover information
#LAI_SRC    FROM_VEGLIB    # FROM_VEGPARAM = read LAI from veg param file; FROM_VEGLIB = read LAI from veg library file
#ALB_SRC    FROM_VEGLIB    # FROM_VEGPARAM = read albedo from veg param file; FROM_VEGLIB = read albedo from veg library file
#FCAN_SRC   FROM_VEGLIB    # FROM_VEGPARAM = read fcanopy from veg param file; FROM_VEGLIB = read fcanopy from veg library file
SNOW_BAND   1   # Number of snow bands; if number of snow bands > 1, you must insert the snow band path/file after the number of bands (e.g. SNOW_BAND 5 my_path/my_snow_band_file)

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
#   OUTFILE       <prefix>
#   AGGFREQ       <freq>            <value>
#   COMPRESS      <compress>
#   OUT_FORMAT    <file_format>
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#
#   OUTFILE       <prefix>
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#   OUTVAR        <varname>       [ <format>       [ <type>  [ <multiplier>   [ <aggtype>]]]]
#
#
# where
#   <prefix>     = name of the output file, NOT including latitude
#                  and longitude
#   <freq>       = Describes aggregation frequency for output stream. Valid
#                  options for frequency are:
#                    NEVER     = never write to history file
#                    NSTEPS    = write to history every <value> steps
#                    NSECONDS  = write to history every <value> seconds
#                    NMINUTES  = write to history every <value> minutes
#                    NHOURS    = write to history every <value> hours
#                    NDAYS     = write to history every <value> days
#                    NMONTHS   = write to history every <value> months
#                    NYEARS    = write to history every <value> years
#                    DATE      = write to history on the date: <value>
#                    END       = write to history at the end of the simulation
#   <value>      = integer or date string (YYYY-MM-DD) describing the number
#                  of <freq> intervals to pass before writing to the history file.
#   <compress>   = gzip compression option.  TRUE, FALSE, or integer between 1-9.
#   <varname>    = name of the variable (this must be one of the
#                  output variable names listed in vic_driver_shared_all.h.)
#   <format>     = (for ascii output files) fprintf format string,
#                  e.g.
#                    %.4f = floating point with 4 decimal places
#                    %.7e = scientific notation w/ 7 decimal places
#                    *    = use the default format for this variable
#
#   <format>, <type>, and <multiplier> are optional.  For a given
#   variable, you can specify either NONE of these, or ALL of
#   these.  If these are omitted, the default values will be used.
#
#   <type>       = (for binary output files) data type code.
#                  Must be one of:
#                    OUT_TYPE_DOUBLE = double-precision floating point
#                    OUT_TYPE_FLOAT  = single-precision floating point
#                    OUT_TYPE_INT    = integer
#                    OUT_TYPE_USINT  = unsigned short integer
#                    OUT_TYPE_SINT   = short integer
#                    OUT_TYPE_CHAR   = char
#                    *               = use the default type
#   <multiplier> = (for binary output files) factor to multiply
#                  the data by before writing, to increase precision.
#                    *    = use the default multiplier for this variable
#   <aggtype>    = Aggregation method to use for temporal aggregation. Valid
#                  options for aggtype are:
#                    AGG_TYPE_DEFAULT = default aggregation type for variable
#                    AGG_TYPE_AVG     = average over aggregation window
#                    AGG_TYPE_BEG     = beginning of aggregation window
#                    AGG_TYPE_END     = end of aggregation window
#                    AGG_TYPE_MAX     = maximum in aggregation window
#                    AGG_TYPE_MIN     = minimum in aggregation window
#                    AGG_TYPE_SUM     = sum over aggregation window
#
#######################################################################
```
