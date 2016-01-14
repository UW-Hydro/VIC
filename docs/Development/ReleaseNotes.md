# Release Notes

This is the list of changes to VIC between each release.  For full details, see the commit logs at [https://github.com/UW-Hydro/VIC](https://github.com/UW-Hydro/VIC).

### Known Issues in Current Release

For a list of known issues and their fixes (in bug-fix updates), visit the VIC GitHub [Issues](https://github.com/UW-Hydro/VIC/issues) page.

### Version Check

To check which release of VIC you are running:

Type `vicNl -v` or for VIC 5 and later, `vic_{driver}.exe -v`.

------------------------------

## VIC 5.0.0

**Release date: (Unreleased)**

This is a major update from VIC 4. The VIC 5.0.0 release aims to have nearly identical physics as VIC 4.2 while providing a clean, refactored code base supporting multiple drivers. There are a number of backward incompatible changes. See the VIC Github page for more details on the changes included in this release.

#### New Features:

1. "vic_run"

	Although the physics and model behavior of VIC 5.0.0 should be nearly identical to VIC 4.2, the source code has undergone a major cleanup and reorganization. We have separated the physical core ("vic_run") from the driver source code. This work has improved the extensibility and readability of the model.

2. Classic Driver

	The Classic Driver provides similar functionality as VIC 4, including ASCII and binary I/O, and a time-before-space evaluation loop order. The classic driver is maintained for two main reasons: 1) to provide some level of backward compatibility for existing VIC users that wish to continue using VIC using a traditional approach, and 2) to allow VIC to be run at individual grid cells, without requiring the infrastructure needed by the Image Driver. Classic Driver specific documentation can be found here.

3. Image Driver

	The Image Driver adds a number of features to the user interface of the VIC model. Most notably, it uses a space-before-time evaluation loop order, netCDF I/O, and parallelization using Open-MPI.  Image Driver specific documentation can be found here.

4. Constants File ([GH#192](https://github.com/UW-Hydro/VIC/pull/173))

	Earlier versions of VIC included many hard-coded parameters and constants.  We have consolidated these constants into a single structure and developed a input file that allows users to modify parameters at run-time.  See here for more information.

5. Logging ([GH#173](https://github.com/UW-Hydro/VIC/pull/173))

	A set of logging Macros have been added to all drivers and `vic_run`. The logging level can be set in the driver `Makefile` via the `LOG_LVL` variable. The logging Macros provide the filename and line number in the source code to aid in debugging.  Additionally, when compiler support is available, a traceback is printed when VIC exits during runtime.  

6. Sub-hourly Timestep ([GH#188](https://github.com/UW-Hydro/VIC/pull/188))

	Previous versions of VIC were limited to a minimum timestep of one hour.  The units of the VIC timestep have been changed from hours to seconds and the minimum timestep is now one second. If you intend on running VIC at a timestep less that one hour, we suggest significant testing.

7. Calendar Support ([GH#188](https://github.com/UW-Hydro/VIC/pull/188))

	Earlier versions of VIC used the standard Gregorian calendar.  Because many modern climate models use non-standard calendars, we have implemented all [CF compliant calendars](http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-20010629.htm#cal). The standard calendar remains the VIC default.  See the documentation for individual drivers for how to set the calendar option.

8. Sample Datasets

9. Tests Datasets

#### Backwards Incompatible Changes:

1.  Classic Driver I/O ([GH#227](https://github.com/UW-Hydro/VIC/pull/227) and ?)

	The format of ASCII forcing and output files has changed in VIC 5. These changes were motivated by the desire to improve simulation metadata tracking and reproducibility of VIC simulations.

	- Forcing files now require date stamps for each timestep and a header specifies the names of the forcing variables.   
	- Output files now include a header with simulation metadata and variable names. The `PRT_HEADER` option has been depreciated.

2.  Classic Driver Global Parameter Options:

	A number of global parameter options have changed for the Classic Driver, relative to VIC 4.

	- `TIME_STEP` (int, units: hours) has been changed to `MODEL_STEPS_PER_DAY` (int)
	- `SNOW_STEP` (int, units: hours) has been changed to `SNOW_STEPS_PER_DAY` (int)
	- `OUT_DT` (int, units: hours) has been changed to `OUTPUT_STEPS_PER_DAY` (int)
	- `FORCE_DT` (int, units: hours) has been changed to `FORCE_STEPS_PER_DAY` (int)

#### Depreciated Features:

1.  Removed unused global parameter option `MEASURE_H` ([GH#284](https://github.com/UW-Hydro/VIC/pull/284).)
2.  Removed MTCLIM ([GH#288](https://github.com/UW-Hydro/VIC/pull/288)).

	Previous versions of VIC used MTCLIM to generate missing forcing variables required to run VIC.  This led to confusion by many users and considerably more complex code in the Classic Driver. VIC forcings are now required to be provided at the same time frequency as the model will be run at (`SNOW_STEPS_PER_DAY` or `MODEL_STEPS_PER_DAY`). The following options have been removed from the Classic Driver:

	- `LW_TYPE`
	- `LW_CLOUD`
	- `MTCLIM_SWE_CORR`
	- `VP_INTERP`
	- `VP_ITER`
	- `OUTPUT_FORCE`

	We are providing a stand-alone version of MTCLIM that produces subdaily VIC meteorological forcings.  That tool is available [here](http://mtclim.readthedocs.org).

#### Bug Fixes:

------------------------------

## VIC 4.2.b [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22307.svg)](http://dx.doi.org/10.5281/zenodo.22307)

Source code is available here: [![VIC.4.2.b](https://img.shields.io/badge/VIC-4.2.b-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.b)

**Release date: (January 22, 2015)**

This is a bugfix update from 4.2.a.

!!! Note "Note: Final Release of VIC 4 Development Track"
	This is the last release of the VIC Version 4 development track.  The next release will be VIC.5.0 and will include backward incompatible changes.

#### Bug Fixes:

1.  Fixed memory error in `initialize_atmos` when OUTPUT_FORCE = TRUE. ([GH#201](https://github.com/UW-Hydro/VIC/issues/201))

	Previously, access to unitialized elements of the veg_con and veg_hist structure was attempted when OUTPUT_FORCE = TRUE, causing a memory error and the model to crash.  This fix sets these elements inside a `if (!options.OUTPUT_FORCE)` block allowing the OUTPUT_FORCE option to work as expected.

------------------------------

## VIC 4.2.a

Source code is available here: [![VIC.4.2.a](https://img.shields.io/badge/VIC-4.2.a-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.a)

**Release date: (December 21, 2014)**

This is a bugfix update from 4.2.

#### Bug Fixes:

1.  Fixed uninitialized bare soil albedo. ([GH#178](https://github.com/UW-Hydro/VIC/issues/178))

	Previously, bare_albedo was unset for the bare soil case (`iveg!=Nveg`). This fix sets the bare_albedo to the global variable value of `BARE_SOIL_ALBEDO`.

2.  Cleanup of frozen soil option constraints.

	Removed hardcoded, behind the scenes checks for the `EXP_TRANS` and `NO_FLUX` global parameter values for case of `QUICK_SOLVE=TRUE` in `calc_surf_energy_bal`.

------------------------------

## VIC 4.2

Source code is available here: [![VIC.4.2](https://img.shields.io/badge/VIC-4.2-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2)

**Release date: (November 19, 2014)**

This is a minor release from 4.1.2.  It includes several new features, bug fixes, and a few backwards incompatible changes.

#### New Features:

1.  Added partial vegetation cover (within each tile).

	Added a time-varying (non-climatological) partial vegcover fraction within vegetated tiles.  Previously, vegetation was assumed to cover 100% of the land surface in a vegetated tile (i.e., "big leaf" scheme).  This assumption is not valid in general, and in particular becomes very inaccurate in arid environments (e.g., open shrublands), or as LAI decreases to near 0.  In such cases, evaporation from bare soil between the plants becomes a major (or dominant) component of total evapotranspiration.

	This partial veg cover fraction ("vegcover") is treated the same way as LAI and albedo: new options (`VEGCOVER_SRC` and `VEGPARAM_VEGCOVER`) tell VIC whether the veg param file contains 12 climatological vegcover values and whether to use those or the ones in the veg library.  A new forcing variable ("VEGCOVER") can be included in forcing files along with LAI and albedo, and must be specified in the global parameter file in the same way.  There is an additional option controlling whether `VEGCOVER` values appear in the veg library file: VEGLIB_`VEGCOVER`.  `FALSE` by default, if TRUE it tells VIC to expect 12 monthly `VEGCOVER` values in each veg class, after the 12 LAI and before the 12 albedo values.

	Internally, VIC uses the partial vegcover fraction to divide each veg tile into the area covered by plants and the area in between the plants.  LAI and canopy moisture and snow storage are rescaled by 1/vegcover to get plant-specific values before canopy evap, transpiration, and canopy snow dynamics are computed (vegcover is not allowed to go below the value `MIN_VEGCOVER` in `vicNl_def.h`).  Bare soil evap is computed for the bare soil component of the tile.  Total evapotranspiration is computed as the area-weighted sum of canopy evap and transpiration from the vegetated fraction and bare soil evap from the bare soil fraction.  Finally, LAI and canopy moisture storage are rescaled back to the tile-area-average values before output.

2.  Added non-climatological time-varying veg parameters.

	Added ability to read timeseries of LAI and albedo as forcing variables.  In addition, climatological albedo values can now be given in the veg parameter file in a similar manner to LAI.

	These changes involved adding a new veg_hist data structure to contain the timeseries of LAI, Wdmax, and albedo, as well as adding LAI, Wdmax, and albedo variables to the veg_var structure to store the current values of these variables.

	To allow for climatological albedo values to be specified in the veg parameter file (similar to LAI), a new option has been introduced: ALB_SRC.  The constants LAI_FROM_VEGPARAM and LAI_FROM_VEGLIB have changed name to FROM_VEGPARAM and FROM_VEGLIB so that they can be used to specify the source of either variable.

	VIC now handles LAI and albedo as follows:
		1. If specified as a variable in one of the forcing files, the values from the forcing file will be used in the simulation instead of the values in the veg parameter or veg library files.
		2. If not specified as a variable in one of the forcing files, values will be taken from the veg parameter file if a) they are listed there and b) `LAI_SRC` (for LAI) and/or `ALB_SRC` (for albedo) are set to `FROM_VEGPARAM`.
		3. If not supplied as a forcing and [EITHER not listed in the veg parameter file OR listed there but LAI_SRC or ALB_SRC is set to FROM_VEGLIB], values will be taken from the veg library.


3.  New forcing variables for carbon cycle.

	Added the following new input forcing variables (used for simulations of carbon cycle processes):

	- `CATM`: Atmospheric CO2 mixing ratio [ppm]
	- `FDIR`: Fraction of incoming shortwave that is direct [fraction]
	- `PAR`:  Photosynthetically active radiation [W/m2]

	These variables are optional; if not supplied as forcings, VIC will use default values for them, as follows:

	- `CATM`: Value of CatmCurrent defined in `vicNl_def.h`
	- `FDIR`: Value computed by MTCLIM module
	- `PAR`:  `SHORTWAVE*SW2PAR`, with `SW2PAR` defined in `vicNl_def.h`

	Similarly, added the following new output variables:

	- `OUT_CATM`:   (equals CATM)
	- `OUT_COSZEN`: Cosine of the solar zenith angle, computed by MTCLIM module
	- `OUT_FDIR`:   (equals FDIR)
	- `OUT_PAR`:    (equals PAR)

4.  Added simulation of photosynthesis.

	Added simulation of photosynthesis.  The photosynthesis formulation was taken from the BETHY model (Knorr, 2000), which in turn used the Farquhar model for C3 plants and the Collatz model for C4 plants.  In addition, inhibition of photosynthesis under saturated conditions (as described by Frolking et al, 2002) is allowed for.

	This feature requires several new veg parameters to be in the veg library file:

	- `Ctype`:          Photosynthetic pathway; can be C3 or C4
	- `MaxCarboxRate`:  Maximum carboxlyation rate at 25 deg C (mol(CO2)/m2s)
	- `MaxETransport`:  Maximum electron transport rate at 25 deg C (mol(CO2)/m2s) (C3 plants)
	- `CO2Specificity`: CO2 specificity at 25 deg C (mol(CO2)/m2s) (C4 plants)
	- `LightUseEff`:    Light-use efficiency (mol(CO2)/mol(photons))
	- `NscaleFlag`:     TRUE = nitrogen-scaling factors are applicable to this veg class
	- `Wnpp_inhib`:     Moisture level (fraction of maximum moisture) above which photosynthesis experiencing saturation inhibition, i.e. too wet for optimal photosynthesis; only applies to top soil layer
	- `NPPfactor_sat`:  Photosynthesis multiplier (fraction of maximum) when top soil layer is saturated

	There are several new output variables associated with this feature:

	- `OUT_GPP`:  Gross primary productivity [g C/m2d]
	- `OUT_RAUT`: Autotrophic respiration [g C/m2d]
	- `OUT_NPP`:  Net primary productivity [g C/m2d]
	- `OUT_APAR`: Absorbed PAR [W/m2]

	By default, this feature is turned off.  To turn this feature on, set `CARBON` to `TRUE` in the global parameter file.

	When this feature is turned on, you can choose to compute stomatal resistance via the Jarvis formulation (the formulation used by all previous versions of VIC) or as a function of photosynthetic demand. This is determined by the setting of `RC_MODE` in the global parameter file.  A value of `RC_JARVIS` (which is the default) selects the Jarvis formulation.  A value of `RC_PHOTO` selects the photosynthetic demand formulation.

5. Added simulation of soil carbon storage and fluxes.

	Added simulation of soil carbon storage and fluxes.  This formulation was taken mostly from the LPJ model (Sitch, 2003), which in turn used a Lloyd-Taylor model for the dependence of soil respiration on soil temperature.  The dependence of soil respiration on soil moisture was based on the formulation of Yi et al (2012) but modified to allow a small respiration rate under saturated conditions.

	At this point, we do not simulate the storage of carbon in living biomass. Therefore, the flux of carbon into the soil (litterfall) is set equal to the total NPP of the previous calendar year, spread evenly over the current year.  As in the LPJ model, soil carbon is stored in 3 pools: litter (fast), intermediate, and slow; with associated turnover times of 2.86 y, 33.3 y, and 1,000 y, respectively.  Litterfall enters the litter pool.  Carbon exits the litter pool through respiration (RhLitter).  A fraction (fAir) of this respired carbon is in the form of CO2 and is vented directly to the atmosphere (RhLitter2Atm).  The remainder is sent to the intermediate and slow pools in the proportions fInter and (1-fInter), respectively.  These pools also respire carbon, which is assumed to be in the form of CO2 and vented directly to the atmosphere.

	There are several new output variables associated with this feature:
	- `OUT_RHET`: Total heterotrophic respiration vented to the atmosphere (= RhLitter2Atm+RhInter+RhSlow)  [g C/m2d]
	- `OUT_NEE`:  Net Ecosystem Exchange (= NPP-RHET) [g C/m2d]
	- `OUT_LITTERFALL`: Flux of carbon from living biomass into litter pool [g C/m2d]
	- `OUT_CLITTER`: Carbon density in the litter pool [g C/m2]
	- `OUT_CINTER`: Carbon density in the intermediate pool [g C/m2]
	- `OUT_CSLOW`: Carbon density in the slow pool [g C/m2]

	This feature is part of the carbon cycle, controlled by the setting of the `CARBON` option in the global parameter file.

6. Added soil moisture content for half-space below bottom soil layer

 	VIC's soil thermal profile can extend well below its soil hydrologic layers.  Previously, the moisture content of these soil thermal nodes was set to that of the bottom soil layer.  Now, the moisture content can be set to a user-specified value, SLAB_MOIST_FRACT, defined in `vicNl_def.h`.

#### Depreciated Features:

1.  Removed the `DIST_PRCP` option.

	Removed the `DIST_PRCP` option and all functions and variables associated with it.  Removed the `dist` array dimension from the cell and veg_var data structures.  Renamed the `dist_prcp` data structure to `all_vars`.

#### Bug Fixes:

1.  Fixed division by 0 and nans in output when there is no liquid water available to satisfy evaporative demand

	Previously, runoff() scaled estimated evaporation for each frost subarea by that subarea's portion of available liquid moisture, via summing available liquid moisture over all subareas and computing the ratio of each subarea's moisture to the sum.  There was no check on whether the sum > 0, resulting in the possibility of division by 0 when no liquid moisture is available.  This has been fixed (a check was added).  In addition, the apportionment of evaporation to subareas originally included a check on whether all of the evaporative demand was met, with a warning statement if not.  This check and warning have been removed, since the evaporation values are subsequently modified to reflect what actually evaporated (i.e. it's ok for the initial estimate to exceed available moisture without affecting the water balance).

2.  Fixed negative liquid soil moisture for bare soil conditions

	Previously, runoff() only checked whether total (liquid+ice) soilmoisture was > residual moisture, but not whether liquid soil moisturewas positive.  In some cases, in the bare soil tile, liquid soilmoisture could occasionally go negative.  This has been fixed by addinga check on liquid soil moisture to runoff().

3.  Fixed incorrect handling of case of a mix of cells with and without lakes.

	VIC was neither reading the lake parameter file correctly nor initializing the lake data structures correctly for the case of a mix of cells with and without lakes within a single lake parameter file.  This has been fixed.

4.  Fixed use of tmp_moist array without initialization.

	Fixed use of tmp_moist array without initialization.

5. Fix for crash when FROZEN_SOIL, EXP_TRANS and IMPLICIT all == TRUE

	Extended the "cold nose" hack to the "warm nose" condition, and also extended to cover the IMPLICIT scheme.

	This will be superceded by a more bug-free soil temperature scheme in the next major release of the model.

6.  Fixed incorrect assignment of input forcing variables that are moisture fluxes (all forms of precipitation and channel inflow) when ALMA_INPUT is TRUE.

	When ALMA_INPUT was TRUE, VIC was not rescaling moisture fluxes such as precipitation to an hourly time step correctly in initialize_atmos. This led to incorrect assignment of these fluxes to the atmos array. This has been fixed.

7.  Fixed selection of starting point in forcing file when starting in the middle of a day

	For the case of STARTHOUR not equal to 0, VIC was not finding the correct starting record in the forcing file, due to its missing a check on the hour of the forcing record.  This has been fixed.

8. Fixed incorrect handling of case of a mix of cells with and without lakes.

	VIC was neither reading the lake parameter file correctly nor initializing the lake data structures correctly for the case of a mix of cells with and without lakes within a single lake parameter file.  This has been fixed.

9.  Fixed use of tmp_moist array without initialization.

10. Fixed use of uninitialized soil moisture values on first time step.

	The tmp_moist array, used in initialize_model_state() as an input to compute_runoff_and_asat(), was initialized within an if statement that caused it to be sent to compute_runoff_and_asat() without initialization in some cases.  This has been fixed by moving the initialization of tmp_moist outside the if statement.

11. Fixed errors in forcing disaggregation under certain input cases.

	Fixed bugs in the following cases:

	1. User supplied daily incoming shortwave (not sub-daily)
	2. User supplied daily specific or relative humidity without supplying average daily pressure or temperature, respectively (with which to convert these to daily vapor pressure).

12. Fixed bug in root zone calculation.

	Fixed infinite loop that was occurring when the total  of root zone depths exceeded the total soil depth and one of the root zone boundaries coincided with a soil layer boundary.

13. Fixed error in passing `SensibleHeat` to `func_atmos_energy_bal`.

	Replaced (`*SensibleHeat`) with `SensibleHeat` in argument lists of `root_brent`, `error_print_atmos_energy_bal` and `solve_atmos_energy_bal`.

14. Fixed use of uninitialized variable in cold nose fix for frozen soil

	Fixed use of uninitialized variable in cold nose fix for frozen soil.  Code was attempting to check all nodes for a cold nose, but this check requires accessing the value of the next node, which is undefined when we check the bottom node. Now the code does not check the bottom node (which is unlikely to experience a cold nose anyway).

15. Fixed bug in converting from ALMA_INPUT moisture flux units

	Fixed bug in converting from ALMA_INPUT moisture flux units to traditional units (was multiplying by number of seconds in model step when should have been multiplying by number of seconds in forcing step).

16. Fixed incorrect reporting of canopy energy balance terms.

	VIC was not summing under- and over-story energy fluxes correctly for the case of a forest canopy with snow.  This only affected the reporting of energy balance terms in VIC's output; internal calculations were fine.

	In addition, a bug (sign error in flux summation) in VIC's calculation of the surface energy balance (for output only) was fixed.  This reduces the vast majority of the energy balance errors reported during the course of a typical energy balance simulation.

17. Fixed incorrect summing of rain and snow components of precipitation over grid cell

	The amounts of rainfall and snowfall over the lake (or inundated wetland) were being omitted from the grid cell totals.  This has been fixed.

18. Vapor pressure incorrect if user supplies (QAIR or REL_HUMID) + PRESSURE as input forcings instead of vapor pressure.

	For the cases of the combination of (QAIR or REL_HUMID) plus PRESSURE supplied as input forcings instead of VP, the logic distinguishing between daily and sub-daily supplied PRESSURE was flawed, resulting in incorrect values in both cases.  This has been fixed.

19. Incorrect handling of user-supplied `tskc` (cloud fraction) for `LW_CLOUD==LW_CLOUD_DEARDORFF`

	Previous versions of VIC (before 4.1.2) used a full-sky longwave formulation taken from two formulas in the Bras hydrology text.  For the new Deardorff full-sky longwave formulation, the dependence on cloud fraction is different from the old Bras formulation.   In 4.1.2 (and 4.1.2.a-b), the new Deardorff formulation did not account for the possibility of user-supplied cloud fraction; if the user supplied cloud fraction as an input forcing, the resulting longwave was wrong.  This has been fixed.

20. Changed default settings of `MTCLIM_SWE_CORR` and `LW_TYPE` to reflect best general settings

	In light of the findings of Bohn et al. (2012), we have changed the default setting of `MTCLIM_SWE_CORR` to `FALSE` and of `LW_TYPE` to `LW_PRATA`.  These settings give forcing estimates that are less biased in general.

21. Vapor pressure set to 0 if user supplies (`QAIR` or `REL_HUMID`) + `PRESSURE` as input forcings instead of vapor pressure.

	For the cases of the combination of (QAIR or REL_HUMID) plus PRESSURE supplied as input forcings instead of VP, VIC was supposed to compute VP from (QAIR or REL_HUMID) and PRESSURE, then transfer the computed VP to the atmos data structure.  This transfer was being skipped, and vapor pressure was consequently set to 0 during the simulation.  This has been fixed.

22. Computed longwave sometimes is extremely large at high latitudes.

	Previously (VIC 4.1.2, 4.1.2.a, and 4.1.2.b only), when SHORTWAVE and VP were supplied to VIC as input forcings (and LONGWAVE was NOT supplied as a forcing), the incoming longwave radiation computed by VIC would in rare cases become extremely large.  This happens only at high latitudes in winter when the theoretical clear-sky potential solar radiation is very small.  If the supplied VP was large enough, it could cause the internal variable t_tmax (clear-sky transmittance) to go negative.  This in turn would cause the internal variable t_fmax (cloud transmittance) to go negative as well.  This, finally, would cause computed LONGWAVE values to become extremely large, if the LW_CLOUD method was set to DEARDORFF.  This has been fixed.

23. VIC was unnecessarily requiring WIND to be supplied as an input forcing.

	VIC was requiring `WIND` (or the zonal and meridional components of wind, `WIND_E` and `WIND_N`) to be supplied as an input forcing.  Now VIC allows WIND to be omitted.  If `WIND` is not supplied as an input forcing, VIC will supply a default wind speed, defined in `vicNl_def.h` as `DEFAULT_WIND_SPEED` and currently set to 3.0 m/s.

24. Cloud fraction tskc was not accounting for the case in which observed incoming shortwave is supplied as a forcing.

	In the absence of observations, VIC's estimate of cloud fraction, tskc, is a function of some intermediate quantities that are computed within the MTCLIM algorithm (in mtclim_vic.c).  These intermediate terms can be computed from either observed daily shortwave radiation (if available) or simulated daily shortwave radiation.  The computation of tskc was previously taking place in a part of the code where only the simulated daily shortwave radiation was available.  Thus, tskc would not reflect the actual amount of incident shortwave, even if observed incident shortwave was supplied as a forcing.

	This has been fixed.  The tskc computation has been moved to another location in the code where the observed daily shortwave can be accessed (if supplied by the user as a forcing).

25. Incorrect timing of disaggregated radiation and air temperature when daily forcings are supplied and off_gmt is 0.

	VIC was double-shifting the timeseries of subdaily radiation and temperature in cases in which VIC was *not* given sub-daily incoming shortwave as an input forcing *and* soil parameter "off_gmt" was not set to the local time zone offset (i.e. not set to `longitude*24/360`). This caused incorrect timing of the diurnal cycle of radiation and air temperature.

26. Disaggregated radiation is constant throughout the day when daily incoming shortwave radiation is supplied as an input forcing.

	When VIC was supplied with *daily* incoming shortwave as an input forcing variable, VIC would fail to disaggregate this correctly to sub-daily time steps; it would simply repeat the daily average for every sub-daily time step.  This has been fixed.  Now VIC will compute a diurnal cycle whose daily average matches the supplied daily average.

#### Miscellaneous:

1. Better out-of-box behavior for soil temperature scheme

	Added constraints to help ensure efficient, physically reasonable simulation of the soil temperature profile:

	1. Set default values of IMPLICIT and EXP_TRANS to TRUE.

	2. Made "cold" (no-spinup) initial soil temperatures more consistent with air temperature and bottom boundary temperature.

	3. Added validation of option.Nnodes for EXP_TRANS=TRUE to guarantee that, for the given soil temperature bottom boundary depth "dp" (also known as the damping depth), there are at least 3 nodes within the top 50 cm of the soil column.  This is to constrain errors to a reasonable size.  To satisfy this condition, the following relationship must hold: `Nnodes >= 5*ln(dp+1)+1`

	Some examples:

	| dp(m) | minimum Nnodes |
	|-------|----------------|
	| 4     | 9              |
	| 7     | 12             |
	| 10    | 14             |
	| 25    | 18             |
	| 50    | 21             |

	VIC will exit with an error message to this effect if Nnodes is too small for the given value of dp.

2. Cleanup of compile-time options

	Cleanup of `user_def.h,` either by removing options/settings or moving them to `vicNl_def.h` (becoming run-time options if appropriate). user_def.h has been removed.

	The following options were removed:

	- `NO_REWIND`.  All parameter files must contain the same grid cells in
	  the same order.
	- `QUICK_FS`.
	- `EXCESS_ICE`.
	- `OUTPUT_FORCE_STATS`.
	- `LWAVE_COR`.  If the user wishes to correct any of the forcing
	  variables, the user can do this externally to VIC.

	The following options were moved from compile-time to run-time (in the `options_struct` in `vicNl_def.h`):

	- `LOW_RES_MOIST` (renamed to `LOG_MATRIC`)
	- `OUTPUT_FORCE`
	- `CLOSE_ENERGY`
	- `SPATIAL_SNOW`
	- `SPATIAL_FROST`

	The following terms were moved to `vicNl_def.h` without modification:

	- `VERBOSE` (still a compile-time option)
	- Max array dimensions such as `MAX_VEG`, `MAX_LAYERS`, etc.
	- Constants such as `MAXIT_FE` and `LAI_WATER_FACTOR`

------------------------------

## VIC 4.1.2 and earlier

Changes prior to VIC 4.2 are archived in the VIC Git repository. The change log, prior to VIC.5 is archived [here](https://github.com/UW-Hydro/VIC/blob/VIC.4.2.b/src/ChangeLog).
