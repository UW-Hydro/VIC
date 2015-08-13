# Disaggregation of Daily Precipitation into Sub-Daily Amounts

## General Description of Disaggregation Procedure

A set of three programs and scripts has been developed to do a relatively quick disaggregation of daily precipitation for VIC input files. The default in VIC, when daily data is provided, is to just divide daily total evenly throughout the day. As shown on the plot below, routed flows produced by VIC can be fairly sensitive to this. There are many disaggregation routines available, but the most common are not designed to preserve the hour-specific statistics common to each month or season for a region. These are important features, as seen in plots organized by hour of day for different precipitation bins (Bin 1=0-5; 2=5-10;3=10-15;4=15-20;5=>20 mm/d) and season (1=DJF; 2=MAM;3=JJA;4=SON) ([Rainfall Durations in Arkansas](Images/ARK_P_dur_stats.jpg); [Occurence Hour in Nebraska -- note season 3](Images/NE_P_hour_stats.jpg)). The simple approach used here is to use the hourly data from the nearest precipitation station for determining the CDF of hourly precipitiation occurence and the CDF of hours of duration for events in different rainfall categories. This disaggregation process presumes that the original daily forcing files have been created with the time-of-observation adjustment (included in the gridding process). More specifics are below.

This technique is summarized in the manuscript: Maurer, E.P., A.W. Wood, J.C. Adam, D.P. Lettenmaier, and B. Nijssen, 2002, A Long-Term Hydrologically-Based Data Set of Land Surface Fluxes and States for the Conterminous United States, _J. Climate_ 15(22), 3237-3251\.

## Program Descriptions

There are three programs, bundled into the compressed tarfile [temporal_disagg_precip.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_disagg_precip.tgz):

1.  Script (reformat_prcp.scr) to take the output from the data downloaded from the Earthinfo Hourly Precipitation Cds and place it into a useable format

2.  Program (make_lookup_table.c) to create the statistics look-up tables for each station

3.  Program (disagg_prec.c) to do the final disaggregation

## Downloading the programs and data

C programs and scripts: [temporal_disagg_precip.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_disagg_precip.tgz), also available under the "Preparing Meteorological Data" section of the [download page](../SourceCode/Code.md).

Lookup tables for the continental U.S. (by season): [US_LOOKUP_TABLES.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/US_LOOKUP_TABLES.TAR.gz), also available under the "Preparing Meteorological Data" section of the [download page](../SourceCode/Code.md).

## More Specifics

The first step is to select a list of stations from the Earthinfo hourly precipitation CD. Filtering should look for relatively long record years, and high % coverage. For some reason, downloading the station list does not result in a file including latitudes and longitudes, so select all stations (ctrl-a) and step to the summary view. Select all records again, and download this into a station information file (use default ascii format). You might consider saving this with a .sta extension.

Second, step back to the station list, and then step to the "hourly" data. Select all records and download the data in default ascii format, with a file extension of .hly (These file extensions can be easily changed, but for these are the ones that match this documentation and the scripts). Gzip the .hly files.

Edit the top of reformat_prcp.scr to make sure EXT_DAT and EXT_INF match the hly and sta extensions given above. Also, edit INFO and DAT to define output names for the station information and hourly data files. Create a text file with just the roots of the .hly file names. For example, if you have four data files downloaded from the Cds: ../data/cent_hp.hly.gz ../data/west1_hp.hly.gz ../data/east_hp.hly.gz ../data/west2_hp.hly.gz the file name list would look like:

```
../data/cent_hp
../data/east_hp
../data/west1_hp
../data/west2_hp
```

and now run the script with the command line argument of the name of the file name list.

The data are analyzed separately for 5 daily precipitation bins (0-5mm, 5-10, 10-15, 15-20, >20). The data are also separated by either month or season (DJF, MAM, JJA, SON). Two lookup tables are created from the hourly station data:

1.  a file consisting of the probabilities (hourly CDF) for daily rainfall duration for 1 to 24 hours. Actually, the program counts the number of hours of rain, regardless of whether the rain fell in more than one event during the day.

2.  a file of the probabilities (also an hourly CDF) of rainfall occurence for each of 24 hours. The programs were originally designed for determining rainfall start hour, so do not be confused by the occasional reference to "start time."

To create these lookup tables, edit the run_make_lookup.scr, which oddly enough is used to run the make_lookup.c program. The things to change are the names of the station information and hourly data files created by the reformatting script above, the names of the output files for the two lookup tables, and a flag for whether you want to separate the statistics by season (SEAS_FLG = 1) or by month (SEAS_FLG =anything else):

```
set STA_INF = hly_prcp.inf
set STA_DAT = hly_prcp.dat
set OUT_HR  = hp_hour.lookup
set OUT_DUR = hp_dur.lookup
set SEAS_FLG = 1
```

Lookup tables for the continental U.S. (and the associated station information file), with a seasonal breakdown of statistics, can be downloaded above so it is not necessary to access the CDs for work in the U.S.

Finally, use these lookup tables to disaggregate the daily data to any sub-daily increment (1, 2, 3, 4, 6, 8, or 12 hours). By default the disagg_prec.c program looks for input forcing files with P, Tmax, Tmin, and Wind stored in 2-byte binary format. For ascii format, change the flag at the top of the code. Edit the run_disagg_prec.scr script to have the correct filenames. These should be obvious, but are also described in the script.

```
set STA_INF = hly_prcp.inf
set S_LOOKUP = hp_hour.lookup
set D_LOOKUP = hp_dur.lookup
set IN_DIR = /home/edm/mo_met/forcing/
set OUT_DIR = /home/edm/mo_met/forcing/disag_p/
set DISAG_HR = 3
set SEAS_FLG = 1
set PMULT = 40
set ST_YR = 1949
set ST_MO = 1
```

Now the VIC model can be run with the disaggregated precipitation data. The global file can be easily changed to include the new sub-daily data with the original daily forcing data. An example of a portion of the original global control file, which specifies a daily forcing file with precipitation, Tmax, Tmin, and Wind in 2-byte binary format:

```
FORCING1   /home/edm/mo_met/forcing/data_
N_TYPES         4
FORCE_TYPE      PREC    UNSIGNED      40
FORCE_TYPE      TMAX    SIGNED       100
FORCE_TYPE      TMIN    SIGNED       100
FORCE_TYPE      WIND    SIGNED       100
FORCE_FORMAT    BINARY
FORCE_ENDIAN    LITTLE
FORCE_DT        24
FORCEYEAR       1949  
FORCEMONTH        1
FORCEDAY          1
FORCEHOUR         0
```

which can be modified to include the separate disaggregated (in this case, 3 hourly) precipitation data:

```
FORCING1   /home/edm/mo_met/forcing/data_
N_TYPES         4
FORCE_TYPE      SKIP    UNSIGNED    40
FORCE_TYPE      TMAX    SIGNED     100
FORCE_TYPE      TMIN    SIGNED     100
FORCE_TYPE      WIND    SIGNED     100
FORCE_FORMAT    BINARY
FORCE_ENDIAN    LITTLE
FORCE_DT        24
FORCEYEAR       1949  
FORCEMONTH      1
FORCEDAY        1
FORCEHOUR       0
FORCING2        /home/edm/mo_met/forcing/disag_p/data_
N_TYPES         1
FORCE_TYPE      PREC    UNSIGNED    40
FORCE_FORMAT    BINARY
FORCE_ENDIAN    LITTLE
FORCE_DT        3
FORCEYEAR       1949
FORCEMONTH      1
FORCEDAY        1
FORCEHOUR       0
```
