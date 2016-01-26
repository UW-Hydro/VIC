# VIC Run-Time Options - Global Parameter File

The global parameter file serves two main purposes:

1.  Tells VIC the names, locations, and formats of input and output files
2.  Defines global parameters of the simulation (known as _run-time_ options)

The order of the options in the global parameter file is not important, but the complete option name must be followed by the required option type information. To help in understanding this file, an [example file](#example-global-parameter-file) has been attached at the bottom of this page.

# Define Simulation Parameters

The following options determine the type of simulation that wil3l be performed.

## Main Simulation Parameters

| Name          | Type      | Units     | Description                                                                                                                               |
|------------   |---------  |-------    |---------------------------------------------------------------------------------------------------------------------------------------    |
| NLAYER        | integer   | N/A       | Number of moisture layers used by the model                                                                                               |
| NODES         | integer   | N/A       | Number of thermal solution nodes in the soil column                                                                                       |
| MODEL_STEPS_PER_DAY | integer   | steps     | Number of simulation time steps per day. NOTE: MODEL_STEPS_PER_DAY should be > 4 for FULL_ENERGY=TRUE or FROZEN_SOIL=TRUE.             |
| SNOW_STEPS_PER_DAY  | integer   | steps     | Number of time steps per day used to solve the snow model (if MODEL_STEPS_PER_DAY > 1, SNOW_STEPS_PER_DAY should = MODEL_STEPS_PER_DAY)                 |
| RUNOFF_STEPS_PER_DAY  | integer   | steps     | Number of time steps per day used to solve the runoff model (should be >= MODEL_STEPS_PER_DAY)                 |
| STARTYEAR     | integer   | year      | Year model simulation starts                                                                                                              |
| STARTMONTH    | integer   | month     | Month model simulation starts                                                                                                             |
| STARTDAY      | integer   | day       | Day model simulation starts                                                                                                               |
| EITHER:       |           |           | *Note:* **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but **not both**                                                       |
| NRECS         | integer   | N/A       | Number of time steps over which to run model. ***NOTE: The number of records must be defined such that the model completes the last day.***     |
| OR:           |           |           | *Note:* **either** NRECS or ENDYEAR, ENDMONTH, and ENDDAY must be specified, but **not both**                                                       |
| ENDYEAR       | integer   | year      | Year model simulation ends                                                                                                                |
| ENDMONTH      | integer   | month     | Month model simulation ends                                                                                                               |
| ENDDAY        | integer   | day       | Day model simulation ends                                                                                                                 |
| CALENDAR      | string    | N/A       | Calendar to use, Valid calendars: 'standard', 'gregorian', 'proleptic_gregorian' 'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day' |


## Define Energy Balance Parameters

The following options determine the method of resolving the surface energy balance.

| Name          | Type      | Units             | Description                                                                                                                                                                                                                                                                                   |
|-------------- |--------   |---------------    |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    |
| FULL_ENERGY   | string    | TRUE or FALSE     | Option for computing land surface temperature (soil or snowpack surface). <li>**TRUE** = compute (via iteration) the temperature that balances the surface energy budget.  <li>**FALSE** = set surface temperature equal to air temperature.  <br><br>Default = False.                                                |
| CLOSE_ENERGY  | string    | TRUE or FALSE     | Option for controlling links between the energy balances of the surface and the canopy. <li>**TRUE** = iterate between the canopy and surface energy balances until they are consistent. <li>**FALSE** = compute the surface and canopy energy balances separately, once per time step.  <br><br>Default = FALSE.     |

## Define Soil Temperature Parameters

The following options determine the type of simulation that will be performed.

| Name              | Type          | Units                     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|----------------   |------------   |-----------------------    |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FROZEN_SOIL       | string        | TRUE or FALSE             | Option for handling the water/ice phase change in frozen soils.  <li>**TRUE** = account for water/ice phase change (including latent heat).  <li>**FALSE** = soil moisture always remains liquid, even when below 0 C; no latent heat effects and ice content is always 0. <br><br>Default = FALSE. <br><br>*Note:* to activate this option, the user must **also** set the **FS_ACTIVE** flag to 1 in the soil parameter file for each grid cell where this option is desired. In other words, the user can choose for some grid cells (e.g. cold ones) to compute ice contents and for others (e.g. warm ones) to skip the extra computation. |
| QUICK_FLUX        | string        | TRUE or FALSE             | Option for computing the soil vertical temperature profile.  <li>**TRUE** = use the approximate method described by [Liang et al. (1999)](http://dx.doi.org/10.1029/94JD00483) to compute soil temperatures and ground heat flux; this method ignores water/ice phase changes. <li>**FALSE** = use the finite element method described in [Cherkauer and Lettenmaier (1999)](http://dx.doi.org/10.1029/1999JD900337) to compute soil temperatures and ground heat flux; this method is appropriate for accounting for water/ice phase changes.  <br><br>Default = FALSE (i.e. use [Cherkauer and Lettenmaier (1999)](http://dx.doi.org/10.1029/1999JD900337)) when running FROZEN_SOIL; and TRUE (i.e. use [Liang et al. (1999)](http://dx.doi.org/10.1029/94JD00483)) in all other cases. |
| IMPLICIT          | string        | TRUE or FALSE             | If TRUE the model will use an implicit solution for the soil heat flux equation of [Cherkauer and Lettenmaier (1999)](http://dx.doi.org/10.1029/1999JD900337) (QUICK_FLUX is FALSE), otherwise uses original explicit solution. When QUICK_FLUX is TRUE the implicit solution has no effect.  <br>The user can override this option by setting IMPLICIT to FALSE in the global parameter file. The implicit solution is guaranteed to be stable for all combinations of time step and thermal node spacing; the explicit solution is only stable for some combinations. If the user sets IMPLICIT to FALSE, VIC will check the time step, node spacing, and soil thermal properties to confirm stability. If the explicit solution will not be stable, VIC will exit with an error message. <br><br>Default = TRUE. |
| QUICK_SOLVE       | string        | TRUE or FALSE             | This option is a hybrid of QUICK_FLUX TRUE and FALSE. If TRUE model will use the method described by [Liang et al. (1999)](http://dx.doi.org/10.1029/94JD00483) to compute ground heat flux during the surface energy balance iterations, and then will use the method described in [Cherkauer and Lettenmaier (1999)](http://dx.doi.org/10.1029/1999JD900337) for the final solution step. <br><br>Default = FALSE.   |
| NOFLUX            | string        | TRUE or FALSE             | If TRUE model will use a no flux bottom boundary with the finite difference soil thermal solution (i.e. QUICK_FLUX = FALSE or FULL_ENERGY = TRUE or FROZEN_SOIL = TRUE). <br><br>Default = FALSE (i.e., use a constant temperature bottom boundary condition).    |
| EXP_TRANS         | string        | TRUE or FALSE             | If TRUE the model will exponentially distributes the thermal nodes in the [Cherkauer and Lettenmaier (1999)](http://dx.doi.org/10.1029/1999JD900337) finite difference algorithm, otherwise uses linear distribution. (This is only used if FROZEN_SOIL = TRUE).  <br><br>Default = TRUE.   |
| GRND_FLUX_TYPE    | string        | N/A                       | Options for handling ground flux: <li>**GF_406** = use (flawed) formulas for ground flux, deltaH, and fusion as in VIC 4.0.6 and earlier. <li>**GF_410** = use formulas from VIC 4.1.0. <br><br>*NOTE:* this option exists for backwards compatibility with earlier releases and likely will be removed in later releases. <br><br>Default = GF_410. |
| TFALLBACK         | string        | TRUE or FALSE             | Options for handling failures of T iterations to converge.  <li>**FALSE** = if T iteration fails to converge, report an error.  <li>**TRUE** = if T iteration fails to converge, use the previous time step's T value.  <br><br>This option affects the temperatures of canopy air, canopy snow, ground snow pack, ground surface, and soil T nodes. If TFALLBACK is TRUE, VIC will report the total number of instances in which the previous step's T was used, at the end of each grid cell's simulation. In addition, a time series of when these instances occurred (averaged across all veg tile/snow band combinations) can be written to the output files, using the following output variables: <li>OUT_TFOL_FBFLAG = time series of T fallbacks in canopy snow T solution. <li>OUT_TCAN_FBFLAG = time series of T fallbacks in canopy air T solution. OUT_SNOWT_FBFLAG = time series of T fallbacks in snow pack surface T solution. <li>OUT_SURFT_FBFLAG = time series of T fallbacks in ground surface T solution. <li>OUT_SOILT_FBFLAG = time series of T fallbacks in soil node T solution (one time series per node). <br><br>Default = TRUE. |
| SPATIAL_FROST     | string <br>(+integer) | string: TRUE or FALSE <br> integer: N/A  | Option to allow spatial heterogeneity in soil temperature: <li>**FALSE** = Assume soil temperature is horizontally constant (only varies with depth).  <li>**TRUE** = Assume soil temperatures at each given depth are distributed horizontally with a uniform (linear) distribution, so that even when the mean temperature is below freezing, some portion of the soil within the grid cell at that depth could potentially be above freezing. This requires specifying a frost slope value as an extra field in the soil parameter file, so that the minimum/maximum temperatures can be computed from the mean value. The maximum and minimum temperatures will be set to mean temperature +/- frost_slope. <li> If TRUE is specified, you must follow this with an integer value for Nfrost, the number of frost sub-areas (each having a distinct temperature). <br><br>Default = FALSE. |

## Precipitation (Rain and Snow) Parameters

Generally these default values do not need to be overridden.

| Name              | Type                  | Units                     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
|------------------ |-------------------    |-----------------------    |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    |
| SNOW_DENSITY      | string                | N/A                       | Options for computing snow density: <li>**DENS_BRAS** = Use traditional VIC algorithm taken from Bras, 1990. <li>**DENS_SNTHRM** = Use algorithm taken from SNTHRM model.<br><br>Default = DENS_BRAS. |
| BLOWING           | string                | TRUE or FALSE             | If TRUE, compute evaporative fluxes due to blowing snow. <br><br>Default = FALSE. |
| COMPUTE_TREELINE  | string or integer     | FALSE or veg class id     | Options for handling above-treeline vegetation: <li>**FALSE** = Do not compute treeline or replace vegetation above the treeline.  <li>**CLASS_ID** = Compute the treeline elevation based on average July temperatures; for those elevation bands with elevations above the treeline (or the entire grid cell if SNOW_BAND == 1 and the grid cell elevation is above the tree line), if they contain vegetation tiles having overstory, replace that vegetation with the vegetation having id CLASS_ID in the vegetation library. <br><br>*NOTE 1*: You MUST supply VIC with a July average air temperature, in the optional [July_Tavg](SoilParam.md#July_Tavg) field, AND set the [JULY_TAVG_SUPPLIED](#JULY_TAVG_SUPPLIED) option to TRUE so that VIC can read the soil parameter file correctly. <br><br>**NOTE 2**: If LAKES=TRUE, COMPUTE_TREELINE MUST be FALSE. <br>Default = FALSE.|
| CORRPREC          | string                | TRUE or FALSE             | If TRUE correct precipitation for gauge undercatch.  <br><br>***NOTE: This option is not supported when using snow/elevation bands.*** <br><br>Default = FALSE. |
| MAX_SNOW_TEMP     | float                 | deg C                     | Maximum temperature at which snow can fall. <br><br>Default = 0.5 C. |
| MIN_RAIN_TEMP     | float                 | deg C                     | Minimum temperature at which rain can fall. <br><br>Default = -0.5 C. |
| SPATIAL_SNOW      | string                | TRUE or FALSE             | Option to allow spatial heterogeneity in snow water equivalent (yielding partial snow coverage) when the snow pack is melting: <li>*FALSE* = Assume snow water equivalent is constant across grid cell. <li>*TRUE* = Assume snow water equivalent is distributed horizontally with a uniform (linear) distribution, so that some portion of the grid cell has 0 snow pack. This requires specifying the max_snow_distrib_slope value as an extra field in the [soil parameter file](SoilParam.md). <br><br>NOTE: max_snow_distrib_slope should be set to twice the desired minimum spatial average snow pack depth [m]. I.e., if we define depth_thresh to be the minimum spatial average snow depth below which coverage < 1.0, then max_snow_distrib_slope = 2*depth_thresh. <br><br>NOTE: Partial snow coverage is only computed when the snow pack has started melting and the spatial average snow pack depth <= max_snow_distrib_slope/2. During the accumulation season, coverage is 1.0. Even after the pack has started melting and depth <= max_snow_distrib_slope/2, new snowfall resets coverage to 1.0, and the previous partial coverage is stored. Coverage remains at 1.0 until the new snow has melted away, at which point the previous partial coverage is recovered. <br><br>Default = FALSE.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |

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

| Name                  | Type      | Units             | Description                                                                                                                                                                                               |
|-------------------    |---------  |---------------    |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| INIT_STATE            | string    | path/filename     | Full path and filename of initial state file. <br><br>*NOTE*: if INIT_STATE is not specified, VIC will take initial soil moistures from the soil parameter file and set all other state variables to a default state.                                             |
| STATENAME             | string    | path/filename     | Path and file prefix of the state file to be created on the specified date. The date within the simulation at which the state is saved will be appended to the file prefix to form a complete file name. <br><br>*NOTE*: if STATENAME is not specified, VIC will not save its state in a statefile.                                                                                                                          |
| STATEYEAR             | integer   | year              | Year at which model simulation state should be saved. <br><br>*NOTE*: if STATENAME is not specified, STATEYEAR will be ignored.                                                                                                                                           |
| STATEMONTH            | integer   | month             | Month at which model simulation state should be saved. <br><br>*NOTE*: if STATENAME is not specified, STATEMONTH will be ignored.                                                                                                                                          |
| STATEDAY              | integer   | day               | Day at which model simulation state should be saved. State will be saved at the end of the final timestep on this day. <br><br>*NOTE*: if STATENAME is not specified, STATEDAY will be ignored.                                                                                                                                            |
| BINARY_STATE_FILE     | string    | TRUE or FALSE     | If FALSE, VIC reads/writes the intial/output state files in ASCII format. If TRUE, VIC reads/writes intial/output state files in binary format. <br><br>*NOTE*: if INIT_STATE or STATENAME are not specified, BINARY_STATE_FILE will be ignored.                                                                                                                    |

# Define Meteorological and Vegetation Forcing Files

This section describes how to define the forcing files needed by the VIC model.  VIC handles vegetation historical timeseries (LAI, albedo, partial vegetation cover fraction) similarly to meteorological forcings (with some exceptions; see below).

Unlike model parameters, for which 1 file contains data for all grid cells, the meteorological forcings are stored as a separate time series for each grid cell. The time step length of the input forcings must match the time step length at which VIC is running. Input files can be ASCII or Binary (signed or unsigned short ints) column formatted. Columns in the file must be in the same order as they are defined in the global control file.

VIC will allow forcing data to be stored in two different files per grid cell (e.g., precip and wind speed in one file, tmin and tmax in another file; or meteorological variables in one file, and vegetation timeseries in another file). Note that if you are using two forcing files per grid cell, the parameters for the first file must be defined before those for the second. **Bold** numbers indicate the order in which these values should be defined, after each forcing file (FORCING1 or FORCING2). Options that do not have a bold number apply to both forcing file types and should appear after the numbered options.

All FORCING filenames are actually the pathname, and prefix for gridded data types: ex. DATA/forcing_YY.YYY_XX.XXX. Latitude and longitude index suffix is added by VIC based on GRID_DECIMAL parameter defined above, and the latitude and longitude values defined in the [soil parameter file](SoilParam.md).

| Name              | Type      | Units                     | Description                                                                                                                                                                                                                                                                                                                           |
|------------------ |---------  |-------------------------- |------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   |
| (1*) FORCING1     | string    | pathname and file prefix  | First forcing file name, always required. ***This must precede all other forcing parameters used to define the first forcing file.***                                                                                                                                                                                                 |
| (1*) FORCING2     | string    | pathname and file prefix  | Second forcing file name, or FALSE if only one file used. ***This must precede all other forcing parameters used to define the second forcing file, and follow those used to define the first forcing file.***                                                                                                                        |
| (2) FORCE_FORMAT  | string    | BINARY or ASCII           | Defines the format type for the forcing files.                                                                                                                                                                                                                                                                                        |
| (3)FORCE_ENDIAN   | string    | BIG or LITTLE             | Identifies the architecture of the machine on which the binary forcing files were created:  <li>**BIG** = big-endian (e.g. SUN).  <li>**LITTLE** = little-endian (e.g. PC/linux). Model will identify the endian of the current machine, and swap bytes if necessary. Required for binary forcing file, not used for ASCII forcing file. |
| (4) N_TYPES       | int       | N/A                       | Number of columns in the current data file, with the following exception: for the vegetation history variables ALBEDO, LAI_IN, and VEGCOVER, there must be multiple columns for these variables, one per vegetation tile. In this case, ALBEDO, LAI_IN, and VEGCOVER each count as only 1 variable despite covering multiple columns.    |
| (5) [FORCE_TYPE](ForcingData.md) | string<br>string<br>float | VarName<br>(un)signed<br>multiplier | Defines what forcing types are read from the file, and in what order. For ASCII file only the forcing type needs to be defined, but for Binary file each line must also define whether the column is SIGNED or UNSIGNED short int and by what factor values are multiplied before being written to output. Note: Unlike other variables, ALBEDO, LAI_IN, and VEGCOVER, each span multiple columns, one column per veg tile.  This will generally vary from one grid cell to the next as the number of veg tiles varies.  However, ALBEDO, LAI_IN, and VEGCOVER should each have only one FORCE_TYPE entry.  [Click here for details.](ForcingData.md) |
| (6) FORCE_STEPS_PER_DAY | integer   | steps                     | Number of timesteps per day in forcing file (must be >= 1)                                                                                                                                                                                                                                                                                  |
| (7) FORCEYEAR     | integer   | year                      | Year meteorological forcing files start                                                                                                                                                                                                                                                                                               |
| (8) FORCEMONTH    | integer   | month                     | Month meteorological forcing files start                                                                                                                                                                                                                                                                                              |
| (9) FORCEDAY      | integer   | day                       | Day meteorological forcing files start                                                                                                                                                                                                                                                                                                |
| GRID_DECIMAL      | integer   | N/A                       | Number of decimals to use in gridded file name extensions                                                                                                                                                                                                                                                                             |
| WIND_H            | float     | m                         | Height of wind speed measurement over bare soil and snow cover. ***Wind measurement height over vegetation is now read from the vegetation library file for all types, the value in the global file only controls the wind height over bare soil and over the snow pack when a vegetation canopy is not defined.***                   |
| ALMA_INPUT        | string    | TRUE or FALSE             | This option tells VIC the units to expect for the input variables:  <li>**FALSE** = Use standard VIC units: for moisture fluxes, use cumulative mm over the time step; for temperature, use degrees C;  <li>**TRUE** = Use the units of the ALMA convention: for moisture fluxes, use the average rate in mm/s (or kg/m<sup>2</sup>s) over the time step; for temperature, use degrees K;  <br><br>Default = FALSE. |

- If using one forcing file, use only FORCING1, if using two forcing files, define all parameters for FORCING1, and then define all forcing parameters for FORCING2\. All parameters need to be defined for both forcing files when a second file is used.

_Examples._ a standard four column daily forcing data file will be defined as:

## ASCII File

    FORCING1   FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     4
    FORCE_TYPE  PREC
    FORCE_TYPE  WIND
    FORCE_FORMAT    ASCII
    FORCE_STEPS_PER_DAY    24

## Binary File

    FORCING1   FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     4
    FORCE_TYPE  PREC    UNSIGNED    40
    FORCE_TYPE  WIND    SIGNED      100
    FORCE_FORMAT    BINARY
    FORCE_ENDIAN    LITTLE
    FORCE_STEPS_PER_DAY    24

## Vegetation Timeseries Variables

For each variable, there must be a separate column for each vegetation tile in the grid cell (which generally will vary from one grid cell to the next).  For example, if there are 3 vegetation tiles in a particular grid cell; and you wish to supply VIC with LAI, partial vegetation cover fraction, and albedo; the input file for the given cell should look like:

    LAI1 LAI2 LAI3 VEGCOVER1 VEGCOVER2 VEGCOVER3 ALBEDO1 ALBEDO2 ALBEDO3

where the 1, 2, and 3 correspond to the first, second, and third tiles listed in the vegetation parameter file, respectively; and the file should be described in the global parameter file as:

    FORCING2    FORCING_DATA/veg_hist/veg_hist__
    N_TYPES     3
    FORCE_TYPE  LAI_IN
    FORCE_TYPE  VEGCOVER
    FORCE_TYPE  ALBEDO
    FORCE_FORMAT    ASCII
    FORCE_STEPS_PER_DAY    1
    FORCEYEAR   1950
    FORCEMONTH  1
    FORCEDAY    1

NOTE that N_TYPES is 3 in the example above, not 9.  This is because N_TYPES only counts the number of different variable types, NOT the total number of columns.

# Define Parameter Files

The following options describe the input parameter files.

| Name                  | Type      | Units             | Description                                                                                                                                                                                                                                                                                               |
|--------------------   |---------  |---------------    |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| SOIL                  | string    | path/filename     | the Soil parameter file.                                                                                                                                                                                                                                                                                  |
| BASEFLOW              | string    | N/A               | This option describes the form of the baseflow parameters in the soil parameter file: <li>**ARNO** = fields 5-8 of the soil parameter file are the standard VIC baseflow parameters <li>**NIJSSEN2001** = fields 5-8 of the soil parameter file are the baseflow parameters from Nijssen et al (2001) <br><br>Default = ARNO. |
| JULY_TAVG_SUPPLIED    | string    | TRUE or FALSE     | If TRUE then VIC will expect an additional column (July_Tavg) in the soil parameter file to contain the grid cell's average July temperature. If your soil parameter file contains this optional column, you MUST set JULY_TAVG_SUPPLIED to TRUE so that VIC can read the soil parameter file correctly. <br><br>*NOTE*: Supplying July average temperature is only required if the COMPUTE_TREELINE option is set to TRUE. <br><br>Default = FALSE. |
| ORGANIC_FRACT         | string    | TRUE or FALSE     | (release 4.1.2 and later) <li>**TRUE** = the soil parameter file contains `3*Nlayer` extra columns, listing, for each layer: the organic fraction, and the bulk density and soil particle density of the organic matter in the soil layer. <li>**FALSE** = the soil parameter file does not contain any information about organic soil, and organic fraction should be assumed to be 0. <br><br>Default = FALSE. |
| VEGLIB                | string    | path/filename     | Vegetation library file name                                                                                                                                                                                                                                                                              |
| VEGPARAM              | string    | path/filename     | Vegetation parameter file name                                                                                                                                                                                                                                                                            |
| ROOT_ZONES            | integer   | N/A               | Number of defined root zones defined for root distribution.                                                                                                                                                                                                                                               |
| VEGPARAM_ALBEDO       | string    | TRUE or FALSE     | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly ALBEDO values for each vegetation type for each grid cell. <br><br>Default = FALSE. |
| ALBEDO_SRC            | string    | N/A               | This option tells VIC where to look for ALBEDO values: <li>**FROM_VEGLIB** = Use the ALBEDO values listed in the vegetation library file. <li>**FROM_VEGPARAM** = Use the ALBEDO values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_ALBEDO must be TRUE. <li>**FROM_VEGHIST** = Use the ALBEDO values listed in the veg_hist forcing files. Note: for this to work, ALBEDO must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. <br><br>Default = FROM_VEGLIB. |
| VEGPARAM_LAI          | string    | TRUE or FALSE     | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly LAI values for each vegetation type for each grid cell. <br><br>Default = FALSE. |
| LAI_SRC               | string    | N/A               | This option tells VIC where to look for LAI values: <li>**FROM_VEGLIB** = Use the LAI values listed in the vegetation library file. <li>**FROM_VEGPARAM** = Use the LAI values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_LAI must be TRUE. <li>**FROM_VEGHIST** = Use the LAI values listed in the veg_hist forcing files. Note: for this to work, LAI_IN must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. <br><br>Default = FROM_VEGLIB. |
| VEGLIB_VEGCOVER       | string    | TRUE or FALSE     | If TRUE the vegetation library file contains monthly VEGCOVER values for each vegetation type for each grid cell (between the LAI and ALBEDO values). <br><br>Default = FALSE. |
| VEGPARAM_VEGCOVER     | string    | TRUE or FALSE     | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly VEGCOVER values for each vegetation type for each grid cell. <br><br>Default = FALSE. |
| VEGCOVER_SRC          | string    | N/A               | This option tells VIC where to look for VEGCOVER values: <li>**FROM_DEFAULT** = Set VEGCOVER to 1.0 for all veg classes, all times, and all locations. <li>**FROM_VEGLIB** = Use the VEGCOVER values listed in the vegetation library file. Note: for this to work, VEGLIB_VEGCOVER must be TRUE.. <li>**FROM_VEGPARAM** = Use the VEGCOVER values listed in the vegetation parameter file. Note: for this to work, VEGPARAM_VEGCOVER must be TRUE. <li>**FROM_VEGHIST** = Use the VEGCOVER values listed in the veg_hist forcing files. Note: for this to work, VEGCOVER must be supplied in the veg_hist files and listd in the global parameter file as one of the variables in the files. <br><br>Default = FROM_DEFAULT. |
| SNOW_BAND             | integer <br> [string] | N/A <br> [path/filename] | Maximum number of snow elevation bands to use, and the name (with path) of the snow elevation band file. For example: `SNOW_BAND 5 path/filename`.  To turn off this feature, set the number of snow bands to 1 and do not follow this with a snow elevation band file name.  <br><br>Default = 1. |
| CONSTANTS             | string    | path/filename     | Constants / Parameters file name |


# Lake Parameters

The following options only take effect when the lake model is running.

| Name          | Type      | Units                                         | Description                                                                                                                                                                                                                                               |
|-------------- |--------   |---------------------------------------------  |--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   |
| LAKES         | string    | FALSE or path/filename                        | Options for handling lakes: <li>**FALSE** = do not simulate lakes lake parameter path/filename = simulate lakes and read the given file for lake model parameters <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| LAKE_PROFILE  | string    | TRUE or FALSE                                 | Options for describing lake profile: <li>**FALSE** = VIC computes a parabolic depth-area profile for the lake basin <li>**TRUE** = VIC reads user-specified depth-area profile from the lake parameter file <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| EQUAL_AREA    | string    | TRUE or FALSE                                 | Options for computing grid cell areas: <li>**FALSE** = Grid cell boundaries and centers fall on a regular lat-lon grid, i.e. grid cells appear as squares when plotted in geographic projection; this means that the grid cells do not have equal areas. <li>**TRUE** = Grid cells have equal area, i.e. they appear as squares when plotted in an equal-area projection; this means that their boundaries do not fall on a regular lat-lon grid and the cell centers are not equally-spaced in latitude and longitude. <br><br>Default = FALSE.                                                                                                                                                                                                                                          |
| RESOLUTION    | float     | decimal degrees of latitude or area in km<sup>2</sup> | Options for grid cell resolution: <li>If **EQUAL_AREA = FALSE**: width of grid cells, in decimal degrees latitude or longitude. <li>If **EQUAL_AREA = TRUE**: area of grid cells, in km<sup>2</sup>. <br><br>Default = none; this MUST be set by the user to match the grid cell size if the lake model is running.                                                                                                                                                    |

# Define Output Files

The following options describe the output files. Click [here](OutputFormatting.md) for more information.

| Name                  | Type      | Units             | Description                                                                                                                                                                               |
|---------------------- |---------  |---------------    |----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   |
| LOG_DIR            | string    | path name         | Name of directory where log files should be written (optional, default is stdout) |
| RESULT_DIR            | string    | path name         | Name of directory where model results are written                                                                                                                                         |
| OUT_STEP              | integer   | hours             | Output time step length                                                                                                                                                                   |
| SKIPYEAR              | integer   | years             | Number of years to skip before starting to write output file. Used to reduce output by not including spin-up years.                                                                       |
| COMPRESS              | string    | TRUE or FALSE     | if TRUE compress input and output files when done (uses gzip)                                                                                                                             |
| BINARY_OUTPUT         | string    | TRUE or FALSE     | If TRUE write output files in binary (default is ASCII).                                                                                                                                  |
| ALMA_OUTPUT           | string    | TRUE or FALSE     | Options for output units: <li>**FALSE** = standard VIC units. Moisture fluxes are in cumulative mm over the time step; temperatures are in degrees C <li>**TRUE** = units follow the ALMA convention. Moisture fluxes are in average mm/s (kg/m<sup>2</sup>s) over the time step; temperatures are in degrees K <br><br>Default = FALSE. [Click here for more information.](OutputFormatting.md)                                                                                                                                                                         |
| MOISTFRACT            | string    | TRUE or FALSE     | Options for output soil moisture units (default is FALSE): <li>**FALSE** = Standard VIC units. Soil moisture is in mm over the grid cell area <li>**TRUE** = Soil moisture is volume fraction                                                                                                                                                   |
| PRT_HEADER            | string    | TRUE or FALSE     | Options for output file headers (default is FALSE): <li>**FALSE** = output files contain no headers <li>**TRUE** = headers are inserted into the beginning of each output file, listing the names of the variables in each field of the file (if ASCII) and/or the variable data types (if BINARY) <br><br>[Click here for more information.](OutputFormatting.md)                                                                                                                                                          |
| PRT_SNOW_BAND         | string    | TRUE or FALSE     | if TRUE then print snow variables for each snow band in a separate output file (`snow_band_*`). <br><br>*NOTE*: this option is ignored if output file contents are specified. |
| N_OUTFILES\*            | integer   | N/A               | Number of output files per grid cell. [Click here for more information](OutputFormatting.md).                                                                                                                    |
| OUTFILE\*               | <br> string <br>| <br>prefix <br>| Information about this output file: <br>Prefix of the output file (to which the lat and lon will be appended) <br> This should be specified once for each output file. [Click here for more information.](OutputFormatting.md) |
| OUTVAR\*                | <br> string <br> string <br> string <br> integer <br> | <br> name <br> format <br> type <br> multiplier <br> | Information about this output variable:<br>Name (must match a name listed in vic_driver_shared.h) <br> Output format (C fprintf-style format code) <br>Data type (one of: OUT_TYPE_DEFAULT, OUT_TYPE_CHAR, OUT_TYPE_SINT, OUT_TYPE_USINT, OUT_TYPE_INT, OUT_TYPE_FLOAT,OUT_TYPE_DOUBLE) <br> Multiplier - number to multiply the data with in order to recover the original values (only valid with BINARY_OUTPUT=TRUE) <br><br> This should be specified once for each output variable. [Click here for more information.](OutputFormatting.md)|

\* *Note: `N_OUTFILES`, `OUTFILE`, and `OUTVAR` are optional; if omitted, traditional output files are produced. [Click here for details on using these instructions](OutputFormatting.md).*

# Obsolete Options from Earlier Versions

The following options are no longer supported.

| Name          | Type      | Units             | Description                                                                                                                                                                                                                                   |
|-------------  |--------   |-----------------  |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| TIME_STEP     | integer   | hour      | Simulation time step length        |
| SNOW_STEP     | integer   | hour      | Snow model time step length        |
| STARTHOUR     | integer   | hour      | Hour model simulation starts        |
| GRND_FLUX     | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE, compute ground heat flux and energy balance; if FALSE, do not compute ground heat flux. Default: If FULL_ENERGY or FROZEN_SOIL are TRUE, GRND_FLUX is automatically set to TRUE; otherwise GRND_FLUX is automatically set to FALSE.  |
| MIN_LIQ       | string    | TRUE or FALSE     | Version 4.1.1 only. Options for handling minimum soil moisture in presence of ice (default is FALSE): <li>**FALSE** = Use residual moisture as lower bound on soil moisture in Brooks-Corey/Campbell and other relationships involving liquid water. <li>**TRUE** = Use (`residual moisture * unfrozen water fraction` as function of temperature) as lower bound on soil moisture in Brooks-Corey/Campbell and other relationships involving liquid water. |
| GLOBAL_LAI    | string    | TRUE or FALSE     | If TRUE the vegetation parameter file contains an extra line for each vegetation type that defines monthly LAI values for each vegetation type for each grid cell. <br><br>*NOTE*: This option has been replaced by the two options LAI_SRC and VEGPARAM_LAI. |
| OUTPUT_FORCE  | string    | TRUE or FALSE     | If TRUE, perform disaggregation of forcings, skip the simulation, and output the disaggregated forcings. |
| FORCEHOUR    | integer   | hour                      | Hour meteorological forcing files start                             |
| MEASURE_H    | decimal   | m      | Height of humidity measurement        |
| PRT_FLUX      | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print energy fluxes debugging files .  |
| PRT_BALANCE   | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print water balance debugging files .  |
| PRT_SOIL      | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print soil parameter debugging files .  |
| PRT_VEGE      | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print vegetation parameter debugging files .  |
| PRT_GLOBAL    | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print global parameter debugging files .  |
| PRT_ATMOS     | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print forcing data debugging files .  |
| PRT_SNOW      | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print snow debugging files .  |
| PRT_MOIST     | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print soil moisture debugging files .  |
| PRT_TEMP      | string    | TRUE or FALSE     | Versions 4.1.1 and earlier. If TRUE print soil thermal debugging files .  |
| DEBUG_DIR     | string    | char * pathname   | Versions 4.1.1 and earlier.  Debugging files output directory (default directory is the current directory, '.'). |

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
MODEL_STEPS_PER_DAY  8   # number of model time steps per day (set to 1 if FULL_ENERGY = FALSE, set to > 4 if FULL_ENERGY = TRUE)
SNOW_STEPS_PER_DAY  8   # number of time steps per day for which to solve the snow model (should = MODEL_STEPS_PER_DAY if MODEL_STEPS_PER_DAY > 1)
RUNOFF_STEPS_PER_DAY  8   # number of time steps per day for which to solve the runoff model (should be >= MODEL_STEPS_PER_DAY)
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
#           # other options listed in vic_driver_shared.h
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
#BINARY_STATE_FILE       FALSE  # TRUE if state file should be binary format; FALSE if ascii

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
FORCING1    (put the forcing path/prefix here)  # Forcing file path and prefix, ending in "_"
FORCE_FORMAT    BINARY  # BINARY or ASCII
FORCE_ENDIAN    LITTLE  # LITTLE (PC/Linux) or BIG (SUN)
N_TYPES     4   # Number of variables (columns)
FORCE_TYPE  PREC    UNSIGNED    40
FORCE_TYPE  TMAX    SIGNED  100
FORCE_TYPE  TMIN    SIGNED  100
FORCE_TYPE  WIND    SIGNED  100
FORCE_STEPS_PER_DAY    24  # Forcing time step length (hours)
FORCEYEAR   2000    # Year of first forcing record
FORCEMONTH  01  # Month of first forcing record
FORCEDAY    01  # Day of first forcing record
GRID_DECIMAL    4   # Number of digits after decimal point in forcing file names
WIND_H          10.0    # height of wind speed measurement (m)
ALMA_INPUT  FALSE   # TRUE = ALMA-compliant input variable units; FALSE = standard VIC units

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
#VEGLIB_VEGCOVER    FALSE   # TRUE = veg lib file contains 12 monthly values of partial vegcover fraction for each veg class, between the LAI and albedo values
#VEGPARAM_LAI   TRUE    # TRUE = veg param file contains LAI information; FALSE = veg param file does NOT contain LAI information
#VEGPARAM_ALB   FALSE    # TRUE = veg param file contains albedo information; FALSE = veg param file does NOT contain albedo information
#VEGPARAM_VEGCOVER  FALSE    # TRUE = veg param file contains veg_cover information; FALSE = veg param file does NOT contain veg_cover information
#LAI_SRC    FROM_VEGLIB    # FROM_VEGPARAM = read LAI from veg param file; FROM_VEGLIB = read LAI from veg library file
#ALB_SRC    FROM_VEGLIB    # FROM_VEGPARAM = read albedo from veg param file; FROM_VEGLIB = read albedo from veg library file
#VEGCOVER_SRC   FROM_VEGLIB    # FROM_VEGPARAM = read veg_cover from veg param file; FROM_VEGLIB = read veg_cover from veg library file
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
RESULT_DIR      (put the result directory path here)    # Results directory path
OUTPUT_STEPS_PER_DAY   0      # Output interval (hours); if 0, OUT_STEP = MODEL_STEPS_PER_DAY
SKIPYEAR    0   # Number of years of output to omit from the output files
COMPRESS    FALSE   # TRUE = compress input and output files when done
BINARY_OUTPUT   FALSE   # TRUE = binary output files
ALMA_OUTPUT FALSE   # TRUE = ALMA-format output files; FALSE = standard VIC units
MOISTFRACT  FALSE   # TRUE = output soil moisture as volumetric fraction; FALSE = standard VIC units
PRT_SNOW_BAND   FALSE   # TRUE = write a "snowband" output file, containing band-specific values of snow variables; NOTE: this is ignored if N_OUTFILES is specified below.

#######################################################################
#
# Output File Contents
#
# As of VIC 4.0.6 and 4.1.0, you can specify your output file names and
# contents # in the global param file (see the README.txt file for more
# information).
#
# If you do not specify file names and contents in the global param
# file, VIC will produce the same set of output files that it has
# produced in earlier versions, namely "fluxes" and "snow" files, plus
# "fdepth" files if FROZEN_SOIL is TRUE and "snowband" files if
# PRT_SNOW_BAND is TRUE.  These files will have the same contents and
# format as in earlier versions.
#
# The OPTIMIZE and LDAS_OUTPUT options have been removed.  These
# output configurations can be selected with the proper set of
# instructions in the global param file.  (see the output.*.template
# files included in this distribution for more information.)
#
# If you do specify the file names and contents in the global param file,
# PRT_SNOW_BAND will have no effect.
#
# Format:
#
#   N_OUTFILES    <n_outfiles>
#
#   OUTFILE       <prefix>
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#
#   OUTFILE       <prefix>
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#   OUTVAR        <varname>       [<format>        <type>  <multiplier>]
#
#
# where
#   <n_outfiles> = number of output files
#   <prefix>     = name of the output file, NOT including latitude
#                  and longitude
#   <varname>    = name of the variable (this must be one of the
#                  output variable names listed in vic_driver_shared.h.)
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
#
#######################################################################
```
