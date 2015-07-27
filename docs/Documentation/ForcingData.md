# VIC Meteorological Forcings File

## Forcing Data Files

VIC can take either daily (Precip, Tmax, Tmin, Windspeed) or sub-daily meteorological data as inputs (or a combination of daily and sub-daily inputs). VIC has some flexibility in which variables and combinations of variables it can use. The contents and formats of the meteorological data files must be described in the [global parameter file](GlobalParam.md#DefineForcingFiles).

Here is a [List of possible input variables](InputVarList.md) (these can always be found in the file vicNl_def.h).

Click here for a description of [one method of preparing met forcing files](PrepMetInputs.md).

Because VIC allows the user to supply any combination of the necessary met variables, there is no single format for these files. For this reason, the user must give VIC all of the relevant details necessary to read the forcing files (e.g., file type - ascii or binary, file "endian-ness", which variables are in the file and in what order, the units of the variables, and start date of the file). This must be done by inserting this information into the forcings section of the [global parameter file](GlobalParam.md).

VIC needs the following meteorological variables, at the model timestep, to run: precipitation, air temperature, wind speed, atmospheric pressure and density, vapor pressure (or vapor pressure deficit or relative humidity or specific humidity), incoming shortwave (solar) radiation, and incoming longwave (or thermal) radiation. However, VIC can estimate some of these quantities internally, so that the user need not supply all of them.

The minimum set of variables that VIC requires the user to supply are: daily total precipitation (rain and/or snow), daily max and min temperature, and daily average wind speed. In this case, VIC uses the MTCLIM algorithms ([Kimball et al., 1997](References.md#Other); [Thornton and Running, 1999](References.md#Other)) to convert daily min and max temperature to humidity and incoming shortwave radiation. VIC then uses the Tennesee Valley Authority algorithm ([Bras, 1990](References.md#Other)) to deduce incoming longwave radiation from humidity and temperature. VIC also computes atmospheric pressure and density from grid cell elevation and global mean pressure lapse rate. Finally, VIC converts these daily quantities into sub-daily by:

1.  making some assumptions about what time of day the min and max air temperatures occur
2.  interpolating the time series of min and max air temperature with a spline
3.  distributing the shortwave (solar) radiation throughout the day according to solar zenith angle
4.  assuming vapor pressure, atmospheric pressure and density, and wind speed are constant throughout the day
5.  computing sub-daily longwave from sub-daily temperature and constant vapor pressure
6.  computing sub-daily vapor pressure deficit as sub-daily saturation pressure (at sub-daily air temperature) - sub-daily (constant) vapor pressure
7.  dividing daily total precipitation into equal sub-daily amounts

VIC also accepts more than these minimum 4 variables. If the user wishes to supply more than the minimum set of variables, for example sub-daily {precip, air temp, pressure, humidity, wind speed, and sw/lw radiation} terms from a reanalysis product, or sub-daily observations from a meteorological station, this is also OK. In fact, most combinations of variables are acceptable, as long as the user supplies VIC with at least the minimum 4 variables mentioned above. For example, the user could supply daily {precip, tmin, tmax, wind} in one forcing file and sub-daily sw and lw radiation in another file (i.e. the two files need not have the same time step - as long as the sub-daily file is at the same timestep as the model timestep).

**NOTE: Important Information Regarding the Time Zone of Sub-Daily Forcing Inputs**

As of release 4.1.2, we have fixed inconsistencies (present in all previous releases) in timing of sub-daily forcings computed by VIC and sub-daily forcings supplied by the user. Now, we have adopted the following convention:

*   The parameter **off_gmt** (from the [soil parameter file](SoilParam.md)) determines how VIC interprets sub-daily time steps relative to the model start date and time.
*   An **off_gmt** value of **0** indicates that the model start date/time is relative to **Greenwich Mean Time (GMT)**.
*   An **off_gmt** value of **(grid_cell_longitude*24/360)** indicates that the model start date/time is relative to **local time**.
*   When outputting sub-daily results, VIC's output files are referenced to the model start date/time; therefore they are controlled by **off_gmt** (off_gmt=0 means VIC results are referenced to GMT; off_gmt=(grid_cell_longitude*24/360) means VIC results are referenced to local time).
*   **Daily** supplied forcings are assumed to start/end at midnight in **local time**; the forcing start date/time is thus in **local time**. When VIC disaggregates daily forcings into sub-daily forcings, **off_gmt** will be used to determine the time lag between the start of the forcing's diurnal cycle and the start of the VIC simulation.
*   **Sub-daily** supplied forcings are assumed to occur relative to **the time zone indicated by off_gmt**. Therefore, if VIC outputs these sub-daily forcings, they will occur at the exact same time of day as in the input files.

Therefore, if mixing daily and sub-daily forcing inputs, it is important that any sub-daily forcing inputs be shifted as necessary to be in the time zone indicated by off_gmt.

In all cases, the user must list the variables contained in the forcing files, by their names (listed below) and in the exact order in which they appear in the forcing files. This is done in the [global parameter file](GlobalParam.md) using the [FORCE_TYPE](Info/FORCE_TYPE.md) parameter. Possible forcing file data types and units are:

| Variable Name     | Definition                                                                                        | Default Units     | ALMA units        |
|---------------    |-------------------------------------------------------------------------------------------------- |---------------    |-----------------  |
| AIR_TEMP          | sub-daily air temperature                                                                         | deg C             | K                 |
| ALBEDO            | surface albedo                                                                                    | fraction          | fraction          |
| CRAINF            | convective rainfall                                                                               | mm per step       | mm/s == kg/m2/s   |
| CSNOWF            | convective snowfall                                                                               | mm per step       | mm/s == kg/m2/s   |
| DENSITY           | atmospheric density                                                                               | kg/m3             | kg/m2             |
| FDIR              | Fraction of incoming shortwave that is direct                                                     | fraction          | fraction          |
| LAI_IN            | Leaf Area Index (LAI)                                                                             | m2/m2             | m2/m2             |
| LONGWAVE          | incoming longwave (thermal infrared) radiation                                                    | W/m2              | W/m2              |
| LSRAINF           | large-scale rainfall                                                                              | mm per step       | mm/s == kg/m2/s   |
| LSSNOWF           | large-scale snowfall                                                                              | mm per step       | mm/s == kg/m2/s   |
| PAR               | Photosynthetically active radiation                                                               | W/m2              | W/m2              |
| PREC              | total precipitation                                                                               | mm per step       | mm/s == kg/m2/s   |
| PRESSURE          | atmospheric pressure                                                                              | kPa               | Pa                |
| QAIR              | specific humidity                                                                                 | kg/kg             | kg/kg             |
| RAINF             | total rainfall                                                                                    | mm per step       | mm/s == kg/m2/s   |
| REL_HUMID         | relative humidity                                                                                 | fraction          | fraction          |
| SHORTWAVE         | incoming shortwave (solar) radiation                                                              | W/m2              | W/m2              |
| SNOWF             | total snowfall                                                                                    | mm per step       | mm/s == kg/m2/s   |
| TMAX              | daily maximum temperature                                                                         | deg C             | K                 |
| TMIN              | daily minimum temperature                                                                         | deg C             | K                 |
| TSKC              | cloud cover                                                                                       | fraction          | fraction          |
| VEGCOVER          | partial veg cover fraction                                                                        | fraction          | fraction          |
| VP                | atmospheric vapor pressure                                                                        | kPa               | Pa                |
| WIND              | wind speed                                                                                        | m/s               | m/s               |
| WIND_E            | East component of wind speed                                                                      | m/s               | m/s               |
| WIND_N            | North component of wind speed                                                                     | m/s               | m/s               |
| SKIP              | indicates a data column which is ignored, <br>e.g., year, month, and day columns if they are present  |               |                   |

By default, the variables are assumed to have units desribed in column 3 of the table above. To use the units in column 4 (ALMA), the user must set ALMA_INPUT TRUE in the [global parameter file](GlobalParam.md). ALMA units correspond more closely to the units used by reanalysis products or GCMs.

Forcing data files can be in short-int Binary or ASCII column formats. Details are below. Three examples of file type definitions are provided below using a standard daily input file containing precipitation, daily maximum and minimum air temperature and wind speed:

## ASCII Column Format

Below is an example of a 4 column daily forcing file:

    0.000 -10.00 -2.00 2.00
    0.000 -8.50 -1.10 2.00
    5.000 -7.90 0.10 2.00
    10.000 -8.20 -3.10 2.00
    ...


To read this file the global control file must have the following lines:

    FORCING1    FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     4
    FORCE_TYPE  PREC
    FORCE_TYPE  TMAX
    FORCE_TYPE  TMIN
    FORCE_TYPE  WIND
    FORCE_FORMAT    ASCII
    FORCE_DT    24
    FORCEYEAR   1950
    FORCEMONTH  1
    FORCEDAY    1
    FORCEHOUR   0
    FORCING2    FALSE


## Binary short int Format

Short integer binary files are significantly smaller than an equivalent ASCII column file. However, floating point data values must be converted by a multiplication factor in order to be stored in a single short int value. If a forcing value is always positive (such as precipitation), then using UNSIGNED short int values means that a larger range of positive values can be stored. To properly extract the data stored in the forcing file, both the multiplication factor, and whether the data is SIGNED or UNSIGNED must be included in the global control file.

One further piece of information must be provided to make sure that VIC reads the binary file correctly even if the file was created on a diffrent machine than the one which runs the simulation. There are two different archetectures present on our systems. Most of the workstations store data in big-endian format (most significant byte first), while the PC based systems (including all of the FreeBSD systems) store data in little-endian format (least significant byte first). VIC can determine the endian of the simulation machine, but the global control file must tell it the format of the machine on which the binary files were created.

Below is an example of what would be in the global control file in order to read the previous example if it had been encoded in short int binary form:

    FORCING1   FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     4
    FORCE_TYPE  PREC    UNSIGNED    40
    FORCE_TYPE  TMAX    SIGNED      100
    FORCE_TYPE  TMIN    SIGNED      100
    FORCE_TYPE  WIND    SIGNED      100
    FORCE_FORMAT    BINARY
    FORCE_ENDIAN    LITTLE
    FORCE_DT    24
    FORCEYEAR   1950
    FORCEMONTH  1
    FORCEDAY    1
    FORCEHOUR   0
    FORCING2    FALSE


## Using two forcing files

As the final example, assume that you have daily precipitation, and minimum and maximum air temperature, but you also have hourly wind speed. To save space, it makes the most sense to save the data sets separately (otherwise the daily data must be converted to hourly, or the hourly to daily). VIC is capable of reading forcing data from two sources with two different time steps. Below is an example of what the global control file would look like:

    FORCING1    FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     4
    FORCE_TYPE  PREC
    FORCE_TYPE  TMAX
    FORCE_TYPE  TMIN
    FORCE_TYPE  SKIP
    FORCE_FORMAT    ASCII
    FORCE_DT    24
    FORCEYEAR   1950
    FORCEMONTH  1
    FORCEDAY    1
    FORCEHOUR   0
    FORCING2    FORCING_DATA/LDAS_ONE_DEGREE/data_
    N_TYPES     2
    FORCE_TYPE  SKIP
    FORCE_TYPE  WIND
    FORCE_FORMAT    ASCII
    FORCE_DT    1
    FORCEYEAR   1950
    FORCEMONTH  1
    FORCEDAY    1
    FORCEHOUR   0
