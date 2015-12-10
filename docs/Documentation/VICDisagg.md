# How to Use VIC as a Meteorological Forcing Disaggregator

In normal operation, VIC disaggregates daily forcings to a sub-daily time step, either the model time step (TIME_STEP) or the snow time step (SNOW_STEP) if this is different from the model time step. Normally, these disaggregated forcings are used internally but never written to the output files. However, there are many cases in which we want to store these disaggregated forcings; for example, preparing sub-daily forcings for other hydrologic models, or examining the inputs to VIC's physics functions.

For this reason, VIC has an option (OUTPUT_FORCE) to write its internal sub-daily forcings to output files. When OUTPUT_FORCE is TRUE in the [global parameter file](GlobalParam.md), VIC will run as before, but after disaggregating each grid cell's forcings, it will move to the next grid cell without running a simulation.

Here is a detailed procedure for running VIC in OUTPUT_FORCE mode:

## Global Parameter File

Create a [global parameter file](GlobalParam.md). Because VIC will not actually run a simulation, it will ignore some of the file's contents. The only settings that VIC will pay attention to in this mode are:

*   `TIME_STEP`: Set this to the sub-daily interval that you want your forcings disaggregated to
*   `SNOW_STEP`: Set this equal to `TIME_STEP`
*   `STARTYEAR`,` STARTMONTH`, `STARTDAY`, `STARTHOUR`: Set these to the start date for your disaggregated forcings
*   `ENDYEAR`, `ENDMONTH`, `ENDDAY`: Set these to the end date for your disaggregated forcings
*   All the variables in the forcing section should be the same as before, i.e. these describe the **input** (daily) forcings that you are **reading**: `FORCING1`, `FORCING2` (if applicable), `FORCE_FORMAT`, `FORCE_ENDIAN`, `N_TYPES`, `FORCE_TYPE` (there must be one of these for each input variable, e.g. PREC, TMAX, TMIN, WIND), `FORCE_DT`, `FORCEYEAR`, `FORCEMONTH`, `FORCEDAY`, `FORCEHOUR`, `GRID_DECIMAL`, `WIND_H`, `MEASURE_H`, and `ALMA_INPUT`
*   `VEGLIB`,`VEGLIB_VEGCOVER`,`VEGPARAM`,`ROOT_ZONES`,`VEGPARAM_LAI`,`VEGPARAM_VEGCOVER`,`VEGPARAM_ALB`: Set to the same values as you would for a full VIC simulation
*   `RESULT_DIR`: Set to the name of the directory where the disaggregated forcings should be written
*   `OUT_STEP`: Set to 0
*   `ALMA_OUTPUT`: For standard VIC forcings, set to FALSE; for ALMA-compliant forcings (often required by other models) set to TRUE
*   `BINARY_OUTPUT`: Set this to FALSE to produce ASCII forcings, TRUE to produce BINARY forcings
*   `SKIPYEAR`: We recommend setting this to 0
*   `OUTPUT_FORCE`: This must be set to TRUE
*   `N_OUTFILES`, `OUTFILE`, `OUTVAR`: These can be omitted; by default, VIC will produce 1 output file per grid cell, names "full_data_lat_lon_," where lat, lon = latitude and longitude of te grid cell's center. These default output files will contain the following variables:

| Name          | Units (ALMA_OUTPUT FALSE)     | Units (ALMA_OUTPUT TRUE)  |
|-----------    |---------------------------    |-------------------------- |
| PREC          | [mm/timestep]                 | [mm/sec]                  |
| AIR_TEMP      | [deg C]                       | [K]                       |
| SHORTWAVE     | [W/m<sup>2</sup>]                        | [W/m<sup>2</sup>]                    |
| LONGWAVE      | [W/m<sup>2</sup>]                        | [W/m<sup>2</sup>]                    |
| DENSITY       | [kg/m<sup>3</sup>]                       | [kg/m<sup>3</sup>]                   |
| PRESSURE      | [kPa]                         | [Pa]                      |
| VP            | [kPa]                         | [Pa]                      |
| WIND          | [m/s]                         | [m/s]                     |

If you wish to output different forcing variables from the ones wirtten by default, you will need to [specify these in the global parameter file](OutputFormatting.md).

## Run VIC

`vicNl -g global_param_file`

Disaggregated forcings should now reside in the directory you specified in the global parameter file.
