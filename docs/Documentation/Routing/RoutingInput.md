# Routing Model Input Files

The routing model was developed by Dag Lohmann, please refer to the references given below for the methodology.  The model transports grid cell surface runoff and baseflow produced by VIC-Nl within each grid cell to the outlet of that grid cell then into the river system. The within cell routing uses a Unit Hydrograph approach and the channel routing uses the linearized Saint-Venant equation. The river routing model assumes all runoff exits a cell in a single flow direction. For a comprehensive example (for the Stehekin Basin) that contains example routing input files see the file [vic.sample.stehekin.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Datasets/vic.sample.stehekin.tgz), also available under the "Sample Data Sets" section of the [download page](../../SourceCode/Code.md).

* * *

## Input File Setup

The input file contains the parameter file names and parameter values required by the routing model. This file is passed as a command line argument to the routing model (i.e. rout `infile`). Note, some of the parameter files are optional, if these are not provided a constant value must be specified instead.

The format of the input file; and an example are given here:

| MAIN TITLE                                        | # INPUT FILE FOR THE COLUMBIA BASIN. |
|---------------------------------------------------|--------------------------------------|
| TEXT                                              | # NAME OF FLOW DIRECTION FILE        |
| `flow direction file`                             | direc.cmb                            |
| TEXT                                              | # NAME OF VELOCITY FILE              |
| boolean (.TRUE. or .FALSE.)                       | .false.                              |
| `flow velocity file` or  float                    | 1.5                                  |
| TEXT                                              | # NAME OF DIFF FILE                  |
| boolean (.TRUE. or .FALSE.)                       | .false.                              |
| `diffusion file` or  float                        | 800                                  |
| TEXT                                              | # NAME OF XMASK FILE                 |
| boolean (.TRUE. or .FALSE.)                       | .false.                              |
| `xmask file` or  float                            | 25000                                |
| TEXT                                              | # NAME OF FRACTION FILE              |
| boolean (.TRUE. or .FALSE.)                       | .true.                               |
| `contributing fraction file` or  float            | ./rout_input/fraction.cmb            |
| TEXT                                              | # NAME OF STATION FILE               |
| station location file                             | stations.cmb                         |
| TEXT                                              | # PATH OF INPUT FILES AND PRECISION  |
| location of vic input files and prefix            | ./vic/vic_out/fluxes_                |
| No. of decimal places used in VIC input filenames | 3                                    |
| TEXT                                              | # PATH OF OUTPUT FILES               |
| output directory                                  | rout_out/                            |
| TEXT                                              | # MONTHS TO PROCESS                  |
| start and stop year/month of the VIC simulation   | 1969 1 1979 12                       |
| first and last year/month to write output         | 1969 1 1979 12                       |
| TEXT                                              | # NAME OF UNIT HYDROGRAPH FILE       |
| `unit hydrograph file`                            | uh_all                               |

* * *

## Routing Parameter Files

The following files may be required by the model, depending on the flags set in the above input files:

1. Fraction File
2. **Flow Direction File**
3. Flow Velocity File
4. Flow Diffusion File
5. Xmask File
6. Station Location File
7. **UH File**

*Items in **bold** are always required*

The parameter values in the [Flow Velocity File](#flowveloc), [Flow Diffusion File](#difffile), and [UH File](#UH_ALLfile) can all be [calibrated](../Calibration.md#routcal).

The start and stop year and month refer to the period over which the VIC simulation, which is the input to the routing model, was run.
The first and last year refer to the period for which the results of the routing are to be written to output.

## Fraction File

The fraction file is gridded information about the fraction of each grid cell that flows into the basin being routed. This allows the user to more accurately define the basin area, since edge cell can contribute more or less than 100% of their runoff and baseflow components to thebasin.

The format is an arc/info ascii grid. It contains a 6-line header that tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution.

More information on the [Fraction File](Fraction.md) and [creating the file](PrepRoutingParams.md#fractionfile).

## Flow Direction File

The flow direction file tells the routing model how all of the grid cells are connected in the routing net.

The format is an arc/info ascii grid. It contains a 6-line header that tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution.

More information on the [Flow Direction File](FlowDirection.md) and [creating the file](PrepRoutingParams.md#flowdirectionfile).

## Flow Velocity File

This file contains information about the velocity (m/s) for the river routing component of the model.

The format is an arc/info ascii grid. It contains a 6-line header that tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution.

More information on the [Flow Velocity File](FlowVelocity.md), [creating the file](PrepRoutingParams.md#flowvelocityfile), and [calibrating the file](../Calibration.md#routcal).

## Flow Diffusion File

This file contains the flow diffusion (m<sup>2</sup>/s) parameter used in river routing component of the model.

The format is an arc/info ascii grid. It contains a 6-line header that tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution.

More information on the [Flow Diffusion File](FlowDiffusion.md), [creating the file](PrepRoutingParams.md#flowdiffusionfile), and [calibrating the file](../Calibration.md#routcal).

## Fmask File

The values in the xmask file are related to the size (in metres) of a cell.

The format is an arc/info ascii grid. It contains a 6-line header that tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution.

More information on the [Xmask File](Xmask.md) and [creating the file](PrepRoutingParams.md#xmaskfile).

## Ftation Location File

The station location file tells the routing model from which grid cells to produce output flow data. Any number of stations may be defined within the basin, as well as a single basin outlet, where the routing network leaves the defined basin. Each line defining a station is followed by another that tells the routing model whether or not a uh_s file has been generated for the current station location. If set to NONE the routing model generates a new uh_s file in the current directory, otherwise it will read the defined uh_s file.

More information on the [Station Location File](StationLocation.md).

## FH File

This file contains the grid cell impulse response function.

More information on the [UH File](UH.md) and [calibrating the file](../Calibration.md#routcal).

* * *

## Runoff Time Series Files (typically from VIC)

These files provide the time series of daily grid cell runoff and baseflow (one file per grid cell). These are based on the traditional [VIC output flux files](../FluxOutputFiles.md), and are expected to be in ASCII column format, with the following columns:

`YYYY MM DD SKIP SKIP RUNOFF BASEFLOW ...`

Where

*   `YYYY, MM, DD` = 4-digit year, 2-digit month, and 2-digit day
*   `SKIP` = a data column containing any data (in VIC fluxes files, the SKIP columns are typically Precip and Evap); these are ignored by the routing model.
*   `RUNOFF, BASEFLOW` = daily runoff and baseflow from the grid cell, in units of [mm/day]; these will be summed inside th e routing model to arrive at total channel inflow.
*   ... = any number of other data columns, all of which will be skipped by the routing model.

* * *

## Fditing the Source Code for Specific Basins

The [routing code](RunRouting.md) uses hard coded array dimensions. Before [compilation](RunRouting.md) the user should check in rout.f that **NROW** and **NCOL** are greater than or equal to the number of rows and columns specified in the direction file header, and that **NYR** is greater than the number of years over which flows are to be routed.  Also, make sure **PMAX** exceeds the total number of grid cells to be routed. If the dimensions are insufficient a warning will be generated and the program will terminate.
