# Release Notes

This is the list of changes to VIC between each release.  For full details, see the commit logs at [https://github.com/UW-Hydro/VIC](https://github.com/UW-Hydro/VIC).

### Known Issues in Current Release

For a list of known issues and their fixes (in bug-fix updates), visit the VIC GitHub [Issues](https://github.com/UW-Hydro/VIC/issues) page.

### Version Check

To check which release of VIC you are running:

- For VIC 4, type `vicNl -v`
- For VIC 5 and later, type `vic_{classic,image}.exe -v`

------------------------------

## VIC 5.1.0 rc1

<!-- TODO -->
<!-- [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61422.svg)](http://dx.doi.org/10.5281/zenodo.61422) -->

**Release date: (April 27, 2018)**

Source code is available here: [![VIC.5.1.0](https://img.shields.io/badge/VIC-5.1.0-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.5.1.0)

This is a minor update from VIC 5.0.1. The VIC 5.1.0 includes new features, such as a new streamflow routing extension and extended parallelization using OpenMP. The release also includes a number of bug fixes for the CESM driver. See the VIC Github page for more details on the changes included in this release.

#### Model enhancement:

1. Improved calculation of drainage between soil layers ([GH#656](https://github.com/UW-Hydro/VIC/pull/656))

   Drainage from upper layer to adjacent lower layer is calculated according to Brook & Corey curve (where drainage rate is a function of upper-layer soil moisture). In previous versions, a simple numerical solution is applied which uses the timestep-beginning upper-layer soil moisture to calculate drainage rate, and assume this constant rate over the entire timestep. This can cause unreasonably large drainage if the curve has a steep shape and when soil moisture is high. Now, the current version uses exact integral (instead of numerical solution) for layer drainage calculation.

2. Fixes for the CESM driver

   1. [GH#642](https://github.com/UW-Hydro/VIC/pull/642)

    - Using correct fill value datatypes in MPI Gather steps
    - Updated state file name time step to be period-ending rather than period-beginning
    - Set the state file name to the RASM case ID
    - Removed decimal point for missing values for unsigned integers
    - Create dummy forcings when initializing the model (so that there is forcing data for the first time step)
    - Changed pressure units from kPa to Pa
    - Fixed bug that prevented using the correct local domain grid cells in `cesm_put_data.c`
    - Changed reference temperature units from Celsius to Kelvin in `cesm_put_data.c`

   1. [GH#695](https://github.com/UW-Hydro/VIC/pull/695)

    - Fix sign for latent heat fluxes passed from VIC to the coupler
    - Fix sign for longwave radiation passed from VIC to the coupler

   1. [GH#696](https://github.com/UW-Hydro/VIC/pull/696)

    - Changes names of CESM driver functions `trim` and `advance_time` to `trimstr` and `advance_vic_time`, respectively, to avoid conflicts with WRF functions with the same names when compiling RFR case.

   1. [GH#702](https://github.com/UW-Hydro/VIC/pull/702)

    - Fixes Julian day for the first timestep in the dmy struct for the CESM driver.

   1. [GH#710](https://github.com/UW-Hydro/VIC/pull/710)

    - Refactor the cesm_put_data.c routine in the CESM driver to use values from out_data directly, rather than computing them separately in cesm_put_data.c.

   1. [GH#716](https://github.com/UW-Hydro/VIC/pull/716)

    - Fixes initialization of coupler fields and calculates temperature and upwelling longwave to pass to WRF during initialization.

   1. [GH#718](https://github.com/UW-Hydro/VIC/pull/718)

    - Updates the cesm_put_data.c routine in the CESM driver to pass gridcell-averaged albedo to the coupler.

   1. [GH#726](https://github.com/UW-Hydro/VIC/pull/726)

    - Updates the cesm_put_data.c routine in the CESM driver to include the correct units for evap passed to the coupler.

   1. [GH#732](https://github.com/UW-Hydro/VIC/pull/732)

     - Updates the cesm_put_data.c routine in the CESM driver to include the correct units for sensible heat flux and updates the rofliq calculation to be correct (previously only OUT_BASEFLOW was being divided by global_param.dt).

   1. [GH#734](https://github.com/UW-Hydro/VIC/pull/734)

     - Updates the cesm_put_data.c routine in the CESM driver to include the correct signs for turbulent heat fluxes and evaporation. Previously we had switched the signs to agree with the image driver and they should instead be in accordance with the sign conventions for coupled models, which differ from those of land surface models. Also, eliminate populating the `l2x_Sl_ram1` field with aero_resist to agree with the VIC 4 implementation in RASM.

   1. [GH#739](https://github.com/UW-Hydro/VIC/pull/739)

    - Updates the cesm_put_data.c routine in the CESM driver to include the correct signs for the wind stresses and fixes a bug in calculating friction velocity (previously it was missing a square root).

   1. [GH#744](https://github.com/UW-Hydro/VIC/pull/744)

    - Updates the cesm_interface_c.c routine to include missing timers and the VIC RUN timer in the CESM driver.

   1. [GH#746](https://github.com/UW-Hydro/VIC/pull/746)

    - Updates the cesm_interface_c.c routine in the CESM driver to populate the nrecs, endyear, endmonth and endday fields in the global_param struct to make them available to vic_finalize for timing tables (specifically the secs/day columns).  

3. Speed up NetCDF operations in the image/CESM drivers ([GH#684](https://github.com/UW-Hydro/VIC/pull/684))

    These changes speed up image driver initialization, forcing reads, and history writes by only opening and closing each input netCDF file once.

4. Added two new timers to measure time in I/O operations ([GH#703](https://github.com/UW-Hydro/VIC/pull/703))

    These two timers count the CPU and WALL time spent in ``vic_force`` and ``vic_write``. The accumulated time from these timers is printed out at the end of each simulation in the timing table. See also [GH#442](https://github.com/UW-Hydro/VIC/pull/442).

5. Added gridcell-averaged albedo (STATE_AVG_ALBEDO) as a state file variable ([GH#712](https://github.com/UW-Hydro/VIC/pull/712))

    This is for use in the CESM driver for VIC to pass to WRF, but has been implemented in the core structure of VIC (in vic_run) for consistency with the classic and image drivers. Running VIC from a cold start now also includes calculation of gridcell-averaged albedo.

6. Cleanup of the initialization sections of the ``image`` and ``cesm`` drivers ([GH#701](https://github.com/UW-Hydro/VIC/pull/701))

    Codified behavior in the initialization of the ``image`` and `cesm` drivers that requires the parameter variables `AreaFract`, `Pfactor`, `zone_fract`, and `Cv` must sum exactly to 1.0. If using the `SNOW_BAND` option, the area weighted `elevation` must match the mean grid cell elevation (`elev`). VIC will print *warnings* if any of these criteria are violated.  

7. Added thread parallelization using OPENMP ([GH#712](https://github.com/UW-Hydro/VIC/pull/712))

    The VIC image and CESM drivers now may be optionally compiled with OPENMP to enable shared memory thread parallelization. This option should improve the parallel scaling of these drivers by reducing the number of MPI messages and increasing message size.

8. Added streamflow routing extensions ROUT_STUB and ROUT_RVIC for the VIC image driver ([GH#231](https://github.com/UW-Hydro/VIC/pull/231))

    The VIC image driver can be optionally compiled with ROUT_RVIC to enable routing in image mode (ROUT_STUB is the default extension which means no routing). With ROUT_RVIC enabled, the output variable ``OUT_DISCHARGE`` is available, and there will also be an extra state variable ``STATE_ROUT_RING`` stored in the state file.

9. Moved MAX_ITER_GRND_CANOPY, which controls the maximum number of ground-canopy iterations in CLOSE_ENERGY mode for vegetation types with an overstory, to the parameters struct ([GH#771](https://github.com/UW-Hydro/VIC/pull/771))

    Previously this was set in the surface_fluxes.c numerics routine for ground-canopy iterations, which meant that that routine had to be altered to change the maximum number of iterations. It has now been moved to the parameters struct so that it can be overriden in the constants file.

10. Updated new snow density function by adding a cap to new snow density that is set in the parameters struct by the parameter SNOW_NEW_SNOW_DENS_MAX ([GH#776](https://github.com/UW-Hydro/VIC/pull/776))

    Previously the change in cold content of the snowpack term (deltaCC in the snow_data_struct) would get unreasonably large if the Hedstrom and Pomeroy 1998 equation used to calculate snow density, which depends only on air temperature, was calculated with air temperatures above about 2 deg C. We use this term to calculate the ground flux from the snowpack and snow depth, which resulted in extremely small snow depths and unreasonably large ground fluxes from the snowpack (and thus changes in snowpack cold content). Now there is a cap on new snow density with the new parameter SNOW_NEW_SNOW_DENS_MAX as well as a snow depth below which we disregard the ground flux from the snowpack (1.e-8).

10. Miscellaneous clean-up:

    [GH#723](https://github.com/UW-Hydro/VIC/pull/723)

     - Added support for veg_hist forcings (non-climatological) in image mode
     - Fixed erroneous allocation of extra veg tile in image mode
     - Simplified looping over veg tiles and bands in vic_run() and prepare_full_energy()
     - Replaced lengthy data structures with local pointers in vic_run()
     - Simplified out_prec, out_rain, and Melt arrays
     - Updated names of variables and options for LAI and FCANOPY in documentation to match their new names in the code
     - Removed constants MAX_VEG and MAX_BANDS from code; all arrays that were declared with those lengths were replaced with dynamic allocations.  This allowed for specification of veg libraries containing more classes without recompiling the code, and more efficient memory usage.

    [GH#766](https://github.com/UW-Hydro/VIC/pull/766)

     - Improved logic in computing soil evaporation (esoil), primarily in func_surf_energy_bal(), by creating explicit terms for transpiration (transp) and esoil in the layer data structure.

#### Bug Fixes:

1. NetCDF forcing files are now closed at the last timestep in stead of after the last timestep. ([GH#774](https://github.com/UW-Hydro/VIC/pull/774))

2. Renamed "fcov" to "fcan" in image driver to better match variable code name ([GH#673](https://github.com/UW-Hydro/VIC/pull/673))

------------------------------

## VIC 5.0.1

**Release date: (February 1, 2017)**

#### Bug Fixes:

1. Fixed image driver history file name timestamp ([GH#635](https://github.com/UW-Hydro/VIC/pull/635))

	After the fix, the timestamp appeared in the image driver output history filename is the beginning time of the time period in the file.

2. Fixed forceskip rounding bug ([GH#639](https://github.com/UW-Hydro/VIC/pull/639))

	After the fix, the `forceskip` variable in the global parameter structure (i.e., the number of timesteps to skip in the forcing data for the simulation period) is rounded correctly (before the fix, rounding error might cause 1-timestep offset in the simulation results).

3. Fixed a problem with image restarts when using multiple processors ([GH#638](https://github.com/UW-Hydro/VIC/pull/638))

	After the fix, only the master node is assigned the task of validating state file dimensions and coordinate variables. Multiprocessing was also added to the VIC testing framework.

4. Ensured that the mask variable in the input domain file must be integer type; otherwise an error is raised. ([GH#645](https://github.com/UW-Hydro/VIC/pull/645))

5. Fixed a bug related to `make_lastday` function ([GH#647](https://github.com/UW-Hydro/VIC/pull/647))

	Before the fix, the input arguments to function `make_lastday` are sometimes in a wrong order. The bug caused error when trying to write state file on a leap day.

6. Fixed a bug related to writing two-dimensional lat/lon variables to a state file ([GH#652](https://github.com/UW-Hydro/VIC/pull/652))

	Before the bug fix, two-dimensional lat/lon variables were not populated correctly and were written as fill values to a state file. Now two-dimensional lat/lon variables are correctly populated and written.

7. Fixed a bug related to `dz_node` and `node_depth` variables in image driver output state file ([GH#657](https://github.com/UW-Hydro/VIC/pull/657))

	Before the fix, `dz_node` and `node_depth` in image driver output state file were not spatially distributed, which was wrong. Now these two variables are spatially distributed in the output state file.

8. Fixed a bug related to `run_cell` and `mask` variables in image driver inputs ([GH#662](https://github.com/UW-Hydro/VIC/pull/662))

	Before the fix, active cell was controlled by `mask` variable in the domain file in image driver, and `run_cell` variable in the parameter file was not actually used. Now `run_cell` variable in the parameter file controls active cells (`run_cell` must be within the mask defined by the domain file).

9. Fixed a time precision bug for long simulations ([GH#668](https://github.com/UW-Hydro/VIC/pull/668))

	Before the fix, the timestamps of long VIC runs were incorrect in some cases due to precision issue in timestamp generation. This resulted in incorrect output timestamps after running for a long period of time, or output termination. Please refer to [GH#668](https://github.com/UW-Hydro/VIC/pull/668) for details on this bug fix.

10. Fixed a bug related to forcing and simulation start time ([GH#671](https://github.com/UW-Hydro/VIC/pull/671))

	Before the fix, there would be an error if the simulation start time is later than the forcing start time that year AND the simulation spans multiple years. Fixed this bug.

------------------------------

## VIC 5.0.0 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.61422.svg)](http://dx.doi.org/10.5281/zenodo.61422)

**Release date: (September 2, 2016)**

Source code is available here: [![VIC.5.0.0](https://img.shields.io/badge/VIC-5.0.0-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.5.0.0)

This is a major update from VIC 4. The VIC 5.0.0 release aims to have nearly identical physics as VIC 4.2 while providing a clean, refactored code base supporting multiple drivers. There are a number of new features, bug fixes, and backward incompatible changes. See the VIC Github page for more details on the changes included in this release.

#### New Features:

1. "vic_run" ([GH#7](https://github.com/UW-Hydro/VIC/issues/7))

	Although the physics and model behavior of VIC 5.0.0 should be nearly identical to VIC 4.2, the source code has undergone a major cleanup and reorganization. We have separated the physical core ("vic_run") from the driver source code. This work has improved the extensibility and readability of the model.


2. Classic Driver ([GH#7](https://github.com/UW-Hydro/VIC/issues/7))

	The Classic Driver provides similar functionality as VIC 4, including ASCII and binary I/O, and a time-before-space evaluation loop order. The Classic Driver is maintained for two main reasons:

	1. to provide some level of backward compatibility for existing VIC users that wish to continue using VIC using a traditional approach, and,
	2. to allow VIC to be run at individual grid cells, without requiring the infrastructure needed by the Image Driver. Documentation for the Classic Driver can be found [here](../Documentation/Drivers/Classic/ClassicDriver/).


3. Image Driver ([GH#7](https://github.com/UW-Hydro/VIC/issues/7))

	The Image Driver adds a number of features to the user interface of the VIC model. Most notably, it uses a space-before-time evaluation loop order, netCDF I/O, and parallelization using MPI.  Image Driver specific documentation can be found [here](../Documentation/Drivers/Image/ImageDriver/).


4. Constants File ([GH#192](https://github.com/UW-Hydro/VIC/pull/173))

	Earlier versions of VIC included many hard-coded parameters and constants.  We have consolidated these constants into a single structure and developed an input file that allows users to modify parameters at run-time.  See [here](../Documentation/Constants/) for more information.


5. Logging ([GH#173](https://github.com/UW-Hydro/VIC/pull/173))

	A set of logging Macros have been added to all drivers and `vic_run`. The logging level can be set in the driver `Makefile` via the `LOG_LVL` variable. The logging Macros provide the filename and line number in the source code to aid in debugging.  Additionally, when compiler support is available, a traceback is printed when VIC exits during runtime. When the `LOG_DIR` variable is provided in the global parameter file, VIC will write its log(s) to log files instead of printing to stdout.


6. Sub-hourly Timestep ([GH#188](https://github.com/UW-Hydro/VIC/pull/188))

	Previous versions of VIC were limited to a minimum timestep of one hour. The units of the VIC timestep have been changed from hours to seconds and the minimum timestep is now one second. If you intend on running VIC at a timestep of less than one hour, we suggest extensive testing.


7. Calendar Support ([GH#188](https://github.com/UW-Hydro/VIC/pull/188))

	Earlier versions of VIC used the standard Gregorian calendar.  Because many modern climate models use non-standard calendars, we have implemented all [CF compliant calendars](http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-20010629.htm#cal). The standard Gregorian calendar remains the VIC default. See the documentation for individual drivers for how to set the calendar option (e.g. [classic](../Documentation/Drivers/Classic/GlobalParam/#main-simulation-parameters)).


8. Sample Datasets ([GH#387](https://github.com/UW-Hydro/VIC/pull/387))

	The [VIC_sample_data](https://github.com/UW-hydro/VIC_sample_data) repository contains the necessary input datasets (forcings and parameters) to run short simulations of the VIC model for both the classic and image driver.


9. Tests Datasets ([GH#79](https://github.com/UW-Hydro/VIC/issues/79))

	See https://github.com/UW-Hydro/VIC/issues/79 for more information. A temporary location of the test data is here: ftp://ftp.hydro.washington.edu/pub/gergel/VIC5_test_data/


10. Testing and Continuous Integration ([GH#190](https://github.com/UW-Hydro/VIC/pull/190))

	A comprehensive testing platform has been implemented and is available for public use along with the VIC model. A small subset of the test platform is run on [Travis-CI](https://travis-ci.org/UW-Hydro/VIC), which facilitates continuous integration of the VIC test platform. More information on the test platform is [here](Testing.md).


11. Run-time profiling and timing ([GH#442](https://github.com/UW-Hydro/VIC/pull/442))

	A timing module has been added to VIC in order to assess the computational cost and throughput of the VIC model. New output variables (`OUT_TIME_VICRUN_WALL` and `OUT_TIME_VICRUN_CPU`) document the time spent in `vic_run` for each variable. Additionally, a timing table is printed to `LOG_DEST` at the end of each simulation.

#### Backwards Incompatible Changes:

1.  Classic Driver I/O Formatting ([GH#18](https://github.com/UW-Hydro/VIC/issues/18), [GH#204](https://github.com/UW-Hydro/VIC/issues/204), [GH#227](https://github.com/UW-Hydro/VIC/pull/227))

	The format of ASCII forcing and output files has changed in VIC 5. These changes were motivated by the desire to improve simulation metadata tracking and reproducibility of VIC simulations.

	- Output files now include a header with simulation metadata and variable names. The `PRT_HEADER` option has been deprecated.

2.  Classic Driver Global Parameter Options

	A number of global parameter options have changed for the Classic Driver, relative to VIC 4.

	- `TIME_STEP` (int, units: hours) has been changed to `MODEL_STEPS_PER_DAY` (int)
	- `SNOW_STEP` (int, units: hours) has been changed to `SNOW_STEPS_PER_DAY` (int)
	- `OUT_DT` (int, units: hours) has been changed to `OUTPUT_STEPS_PER_DAY` (int)
	- `FORCE_DT` (int, units: hours) has been changed to `FORCE_STEPS_PER_DAY` (int)
	- `BINARY_STATE_FILE` (TRUE or FALSE) has been changed to `STATE_FORMAT` (BINARY or ASCII)
	- `BINARY_OUTPUT` (TRUE or FALSE) has been changed to `OUT_FORMAT` (BINARY or ASCII)

3.  State files now include seconds ([GH#464](https://github.com/UW-Hydro/VIC/pull/464))

	- There is a new global parameter option, `STATESEC`.  This specifies the time step at the end of which state will be saved, in units of seconds.  In other words, if you have an hourly time step (3600 sec) and you want to save state at the end of the final time step of the day (which is 86400 seconds long), subtract 3600 from 86400 to get a STATESEC of 82800.  This corresponds to the first second of the final time step.  State will be saved at the end of that time step.  
	- When the state save date is appended to state filenames, STATESEC will be included so that the date will have the format YYYYMMDD_SSSSS.

4.  Classic Driver Output Variables ([GH#352](https://github.com/UW-Hydro/VIC/pull/352))

	Computation of potential evapotranspiration (PET) has been simplified, reducing the number of output variables from 6 to 1.  The following output variables have been removed:

    - `OUT_PET_SATSOIL` (potential evapotranspiration from saturated bare soil)
    - `OUT_PET_H2OSURF` (potential evapotranspiration from open water)
    - `OUT_PET_SHORT` (potential evapotranspiration (transpiration only) from short reference crop (grass))
    - `OUT_PET_TALL` (potential evapotranspiration (transpiration only) from tall reference crop (alfalfa))
    - `OUT_PET_NATVEG` (potential evapotranspiration (transpiration only) from current vegetation and current canopy resistance)
    - `OUT_PET_VEGNOCR` (potential evapotranspiration (transpiration only) from current vegetation and 0 canopy resistance)

    These have been replaced by:

    - `OUT_PET` (potential evapotranspiration, which = area-weighted sum of potential transpiration and potential soil evaporation; potential transpiration is computed using the Penman-Monteith equation with architectural resistance and LAI of the current veg cover)

#### Deprecated Features:

1.  Removed unused global parameter option `MEASURE_H` ([GH#284](https://github.com/UW-Hydro/VIC/pull/284))
2.  Removed MTCLIM ([GH#288](https://github.com/UW-Hydro/VIC/pull/288))

	Previous versions of VIC used MTCLIM to generate missing forcing variables required to run VIC. This led to confusion by many users and considerably more complex code in the Classic Driver. VIC forcings are now required to be provided at the same time frequency as the model will be run at (`SNOW_STEPS_PER_DAY`).

	As part of this change, the following options have been removed from the Classic Driver:

	- `LW_TYPE`
	- `LW_CLOUD`
	- `MTCLIM_SWE_CORR`
	- `VP_INTERP`
	- `VP_ITER`
	- `OUTPUT_FORCE`

	As part of this change, the following output variables have been removed from the Classic Driver:

	- `OUT_COSZEN`
	- `OUT_TSKC`

	In the future, we would like to provide a stand-alone version of MTCLIM that produces subdaily meteorological forcings. We are looking for community support for this feature ([GH#17](https://github.com/UW-Hydro/VIC/issues/17))

3. Removed `LONGWAVE` and `SHORTWAVE` forcing types ([GH#379](https://github.com/UW-Hydro/VIC/pull/379)).

	Previous versions of VIC allowed users to specify either `LONGWAVE` or `LWDOWN` to denote the incoming longwave radiation flux and `SHORTWAVE` or `SWDOWN` to denote the incoming shortwave radiation flux. We have removed these duplicate options, standardizing on the more descriptive `LWDOWN` and `SWDOWN`.

	Similarly, output variables `OUT_NET_LONG` and `OUT_NET_SHORT` have been replaced with `OUT_LWNET` and `OUT_SWNET`, respectively.

4.  Changed the name of the variable `VEGCOVER` to `FCANOPY`, since this more accurately captures the meaning of the term (i.e., the fractional area of the plant canopy within the veg tile). Similarly changed `OUT_VEGCOVER` to `OUT_FCANOPY`.

    Similarly, changed the names of the following global parameter file options:

    - `VEGLIB_VEGCOVER` --> `VEGLIB_FCAN`
    - `VEGPARAM_VEGCOVER` --> `VEGPARAM_FCAN`
    - `VEGCOVER_SRC` --> `FCAN_SRC`

#### Bug Fixes:

1. Miscellaneous fixes to lake module ([GH#425](https://github.com/UW-Hydro/VIC/pull/425))

	Several lake processes (aerodynamic resistance, albedo, latent/sensible heat fluxes, net radiation, etc) were reported incorrectly or not at all in output files. This has been fixed. In addition, in the absence of an initial state file, lake temperatures were initialized to unrealistic temperatures (the air temperature of the first simulation time step). To fix this, we now initialize the lake temperature to annual average soil temperature.

2. Fix for computation of soil layer temperatures when soil thermal nodes do not reach the bottom of the soil column. ([GH#467](https://github.com/UW-Hydro/VIC/pull/467))

	Previously, if the soil thermal damping depth was shallower than the bottom of the deepest soil layer, and `FROZEN_SOIL==TRUE`, VIC would abort when estimating layer ice contents because it could not estimate a layer temperature if the thermal nodes did not completely span the layer.  Now, a layer temperature is estimated even when thermal nodes do not completely span the layer, and the error no longer occurs.

3. Fix related to exact restart ([GH#481](https://github.com/UW-Hydro/VIC/pull/481), [GH#507](https://github.com/UW-Hydro/VIC/pull/507), [GH#509](https://github.com/UW-Hydro/VIC/pull/509))

	Previously, VIC did not produce the same results (fluxes and states) if a simulation was separated into multiple shorter-period runs by saving the state variables and restarting. This was due to:

	1. The MTCLIM algorithm resulted in slightly different sub-daily meteorological variable values for different lengths of forcings (MTCLIM is deprecated in the current version)

	2. A few bugs resulting in inexact restart. The following bugs have been fixed:

		- The prognostic state variable `energy.Tfoliage` (foliage temperature) is now saved to the state file
		- Two flux variables `energy.LongUnderOut` and `energy.snow_flux` are now saved to the state file.

			!!!Note
					This is a temporary solution to ensure exact restart. A better way of handling these two flux variables needs to be done in the future (see [GH#479](https://github.com/UW-Hydro/VIC/issues/479))

4. Fix for binary state file I/O ([GH#487](https://github.com/UW-Hydro/VIC/pull/487))

	Fixed a bug so that the binary format state file I/O works correctly.

5. Fix for a physical constant (water heat capacity) ([GH#574](https://github.com/UW-Hydro/VIC/pull/574))

	Fixed a bug where volumetric heat capacity of water should be used in `func_canopy_energy_bal` (previously specific heat capacity was used).

------------------------------

## VIC 4.2.d [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.56058.svg)](http://dx.doi.org/10.5281/zenodo.56058)

!!! Note "Note: Final Release of VIC 4 Development Track"
	This is the last release of the VIC Version 4 development track.  The next release will be VIC.5.0 and will include backward incompatible changes. Future updates the VIC 4 development track will be made on the `support/VIC.4.2.d`.

Source code is available here: [![VIC.4.2.d](https://img.shields.io/badge/VIC-4.2.d-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.d)

**Release date: (June 20, 2015)**

This is a bugfix update from 4.2.c.

#### Bug Fixes

1. Fixed uninitialized `dmy_struct` when `OUTPUT_FORCE==TRUE` and `BINARY_OUTPUT==TRUE` ([GH#393](https://github.com/UW-Hydro/VIC/issues/393))

1. Fixed uninitialized vegetation parameters when `VEGPARAM_LAI==FALSE` ([GH#455](https://github.com/UW-Hydro/VIC/issues/455))

------------------------------

## VIC 4.2.c [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.35302.svg)](http://dx.doi.org/10.5281/zenodo.35302)

Source code is available here: [![VIC.4.2.c](https://img.shields.io/badge/VIC-4.2.c-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.c)

**Release date: (December 12, 2015)**

This is a bugfix update from 4.2.b.

#### Bug Fixes:

1.  Documented how VIC 4.2 needs user to specify veg_lib and veg_param files when OUTPUT_FORCE = TRUE. ([GH#305](https://github.com/UW-Hydro/VIC/issues/305))

	Prior to release 4.2, a user could run VIC in OUTPUT_FORCE mode with only a soil parameter file and forcing files.  This functionality is now broken as of release 4.2 and will not be fixed.  Users must either supply veg_lib and veg_parameter files (which the user is likely to have anyway) or use the standalone forcing disaggregator under cevelopment for use with release 5.0.  The documentation was updated to describe this issue as of release 4.2.c.

2. Added architectural resistance of 100 s/m for soil evaporation. ([GH#306](https://github.com/UW-Hydro/VIC/issues/306))

	Testing at approx. 60 eddy covariance towers ([Bohn and Vivoni, 2016](../Documentation/References.md#other-historical-references)) has indicated that soil evaporation is too high with the prior architectural resistance of 0 s/m and too low with a value of 200 s/m.  Further refinement would be ideal but this is a good ballpark figure.

3. Compute aerodynamic conductance of each veg tile as area-weighted average of conductances of vegetated and exposed soil fractions of the tile. ([GH#306](https://github.com/UW-Hydro/VIC/issues/306))

	The prior formulation was not the final version used in [Bohn and Vivoni (2016)](../Documentation/References.md#other-historical-references), but was mistakenly added to the codebase instead of the formulation used here.  This fixes the mistake.

4. Fix overwriting of veg_lib structure with values of current cell in veg_param file. ([GH#319](https://github.com/UW-Hydro/VIC/issues/319))

	Previously, VIC overwrote the LAI, albedo, and vegcover values in the copy of the veg library stored in memory (which is supposed to be constant reference values that apply to all grid cells) with those from the veg_parameter file pertaining to the current grid cell.  Values for veg classes not present in the current grid cell therefore were those of the last grid cell that contained those veg classes.  This did not affect performance but interfered with diagnostics while debugging.

5. Lake parameter validation. ([GH#308](https://github.com/UW-Hydro/VIC/issues/308))

	Previously, there were minimal checks performed on the values of the depth-area relationship.  This allowed unphysical values to be specified, leading to all manner of unphysical behaviors.  This has been fixed.

6. Fix lake water balance errors. ([GH#308](https://github.com/UW-Hydro/VIC/issues/308) and [GH#316](https://github.com/UW-Hydro/VIC/issues/316))

	Previously, precipitation over the lake was scaled by the lake area fraction twice, resulting in water balance errors.  This has been fixed.

------------------------------

## VIC 4.2.b [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22307.svg)](http://dx.doi.org/10.5281/zenodo.22307)

Source code is available here: [![VIC.4.2.b](https://img.shields.io/badge/VIC-4.2.b-blue.svg)](https://github.com/UW-Hydro/VIC/releases/tag/VIC.4.2.b)

**Release date: (January 22, 2015)**

This is a bugfix update from 4.2.a.

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

#### Deprecated Features:

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
