# Preparation of Meteorological Forcing Files

Meteorological Forcing Files may be prepared using a variety of methods. Below you will find the standard method of preparation along with additional processing methods that have been used in our lab.

Obtain precipitation and temperature data (as a minimum) from the [NCDC](http://hurricane.ncdc.noaa.gov/CDO/cdo) or our [web site](MetData.md) and grid it to the final model resolution. Additional meteorological data can also be compiled as input to the VIC model (see e.g. [Meteorological Forcing Files](ForcingData.md)), otherwise these values are parameterized internally to the model. If desired, disaggregate the precipitation into sub-daily increments. This produces a second set of "forcing" files for precipitation.

*   Additional Information:

    *   [Gridding meteorological data from the National Climatic Data Center](GridNCDC.md)

    *   [Forcing files for 2 degree global runs](ftp://ftp.hydro.washington.edu/pub/CE/HYDRO/nijssen/vic_global/index.html)

    *   [Forcing files for 1/8th degree runsfor the continental U.S.](http://www.hydro.washington.edu/Lettenmaier/Data/VIC_retrospective/index.html)

    *   [Disaggregation of daily precipitation](p_disag.md)

# Additional Processing

## Preprocessing of Monthly Wind Fields

## Introduction

This section presents a methodology on how to generate a gridded dataset of monthly windspeeds for the Variable Infiltration Capacity (VIC) model. Programs discussed in this section may be found in the file [met_prep_other.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/met_prep_other.tgz), also available under the "Preparing VIC Model Inputs" section of the [download page](../SourceCode/Code.md). The first step is to extract data from the Surface Airways CDs and then use the data to generate spatially gridded datasets. The interpolation method used is a simple distance averaging method. The final data are then added to your data files as a fourth column. Basic skills in UNIX is required.

## Extracting wind data from Surface Airways CDs

*   First the raw wind data must be extracted from the Surface Airways CD's. Put in the CD and follow these steps.
*   Mark the states you are interested in by pressing SPACE
*   Step to station ENTER on one of the states
*   Filter stations by using F2, FILTER, PARAMETER, CTRL-I, WIND DIRECTION AND SPEED, ENTER
*   Export using F4 into ASCII, ALL STATIONS. This file will be called ^M windstation file.
*   To extract wind data press ENTER and then SUMMARY.
*   Use CTRL-U to get km/h units.
*   Use F4 to export data as ASCII. All records should be extracted.
*   FTP the ASCII files to the UNIX system.
*   You should now have one or several wind station files and one or several wind datafiles.

## Concatenating datafiles

*   Concatenate the wind station file and the wind datafiles. Remember to concatenate the wind station files in the same rder as the wind datafiles.
*   UNIX: cat nvut1.dta wind1.dta wyom1.dta > ! wind_cbr8.dta
*   You should now have only one wind stationfile and one wind datafile.
*   wind_cbr8.dta wind_cbr8.asc

## Format_wind.scr

*   The wind stationfile and the wind datafile is now inputfiles to the next processing program. Format_wind.scr formates the data before the data can be gridded.
*   format_wind.scr windstationfile winddatafile outputfile

## Implement_wind.c

*   The outputfile from format_wind.scr is now input file to a C program that does the final gridding and adds winddata to your VIC metdata. Even the wind data are monthly, daily values are added to the VIC metdata. For each month the same value is added. `Implement_wind.c` also needs a filelist of the metdatafiles.
*   See the #define statements in the top of the program. Filenames and directories are hardcoded. The program have been used for the Ark-Red and the Columbia River at 1/8th degree.

* * *

## Rescaling of Daily Air Temperature

## Introduction

This section presents a short description on how to rescale the air temperature data that was processed using the "REGRID" programs. Several C-programs and AML-programs are used. Programs discussed in this section may be found in the file [met_prep_other.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/met_prep_other.tgz), also available under the "Preparing VIC Model Inputs" section of the [download page](../SourceCode/Code.md).

## Ascending mask

One of the programs that will be used demands a special file called an ascending mask. This maskfile is almost the same as a regular maskfile, but instead of float gridcell values there are integers in ascending order as cell values. The program `ascend.c` will read a maskfile and put ascending integers into the valid cells.

Use: `ascend inmaskfile > outmaskfile`

## Calculate mean monthly airtemperatures

`mk_monthly_airtemp.c` calculates the longterm monthly means from the output file from the "regrid" program (cmb_tmin.grd). 12 monthly asciigrid files are generated and will be used later on.

Look in `run_mk_monthly_airtemp` to set in and out files, etc.

Remember to make the output directories (`./month/`) before you start the program.

## Get PRISM monthly fields

The PRISM fields for your basin have to be extracted. This is done in ArcInfo. get_prism_airtemp.aml is an AML script that calculates the PRISM fields for each active gridcell in your maskfile. Two inputs are needed. The ascending maskfile and the directory where the PRISM asciigrids are stored.

Also from this program there is 12 asciigridded files as output. These files together with the monthly asciigrid files will be used to scale the timeseries.

## Rescale airtemperature

`rescale_airtemp.c` is a program that scales the `*.grd` file from the "regrid" programs. All the timesteps are scaled with the value of `New_Tmin = Tmin + PRISM - Monthly`

The program generates a new outputfile that should have same number of lines as total number of timesteps in your data.

Look in run_rescale_airtemp.scr to see IN and OUT directories, etc..

The same programs are used for TMAX and TMIN.

After rescaling both precipitation and airtemperature `vicinput.c` can be used as before
