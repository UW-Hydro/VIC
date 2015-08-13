# Temporal Aggregation of Hourly Output to Daily

When the model is run in full energy balance mode, the resulting output files contain information for each time step. This can produce large volumes of data even when using binary output formats, it also can produce more data then is really needed. This document describes scripts and programs used to temporally aggregate hourly output to daily quantities which are more manageable, and are used by the routing program, and plotting scripts described in the next few documents. All programs mentioned below may be found in the [temporal_agg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_agg.tgz), also available under the "Post-processing" section of the [download page](../SourceCode/Code.md).

* * *

## Building a File List

The file list contains the file extensions, latitude and longitude of all gridded data files that are to be processed. The file list should appear as shown below:

```
_45.25_-96.25 42.25 -96.25
_45.25_-96.75 42.25 -96.75
_45.75_-96.25 42.75 -96.25
```

Where the first column is the file extension, the second column is the latitude, and the last column is the longitude. To simplify the process of making a file list when you have a large number of grid cells, the script make_file_list.script has been supplied (see the RESULTS directory if you downloaded the sample directory structure). This script reads all of the file extensions from a given directory using a given file prefix (e.g. fluxes), and creates a file called file.list in the current directory. **NOTE:** replace "awk" with "gawk" or "nawk" when running on a SUN or other system where awk does not include all functions.

## Temporally Aggregating ASCII Data

The programs reform_vic_3L_asc.c, reform_vic_3L_snow_asc.c, and reform_vic_3L_fdepth_asc.c (all programs located in RESULTS/SOURCE of the example directory structure) will read through each respective type of [ASCII output data files](FluxOutputFiles.md) temporally aggregating them from hourly to daily values. Output from these programs will be in a similar format, except that the hour column will not appear. Variables like precipitation, runoff, evaporation and baseflow are summed to give daily totals. Taking the last daily value reduces storage terms like soil moisture to daily values that can be used in calculating the water balance. Energy flux terms, like net radiation and latent heat, are averaged for each day.

Each program requires the same information on the command line: the path where the hourly files are located; the file prefix (fluxes, snow, fdepth); the output path for the daily files; the [file list](#filelist); the number of files to process; and the number of days to process. An example command line is shown below:

```reform_vic_3L_asc RES_3L_EWB/ fluxes RES_3L_EWB/DAILY/ file.index 4.060```

## Temporally Aggregating Binary Data

The programs reform_vic_3L_bin.c, reform_vic_3L_snow_bin.c, and reform_vic_3L_fdepth_bin.c (all programs located in RESULTS/SOURCE of the example directory structure) will read through each respective type of [Binary output data files](FluxOutputFiles.md) temporally aggregating them from hourly to daily values. Output from these programs will be the same as that of the [ASCII reducing programs](#asciidata). Variables like precipitation, runoff, evaporation and baseflow are summed to give daily totals. Taking the last daily value reduces storage terms like soil moisture to daily values that can be used in calculating the water balance. Energy flux terms, like net radiation and latent heat, are averaged for each day.

Each program requires the same information on the command line: the path where the hourly files are located; the file prefix (fluxes, snow, fdepth); the output path for the daily files; the [file list](#filelist); the number of files to process; and the number of days to process. An example command line is shown below:

```reform_vic_3L_asc RES_3L_EWB/ fluxes RES_3L_EWB/DAILY/ file.index 4.060```

## Temporally Aggregating LDAS Binary Data

The program extract_daily_means_from_LDAS_binary.c first determines the time step of the LDAS binary file and whether or not it contains frozen soil information. Then it reads through the file and reduces variables from sub-daily to daily values. Precipitation, runoff, evaporation, and baseflow are summed to give daily totals. All other variables are averaged to produce daily average values for each day.

Usage: ```extract_daily_means_from_LDAS_binary <result dir> <result file prefix> <output directory> <file extension list> <number of files> <start date> <end date>```

This program reads the LDAS binary files and extracts daily average or summed parameters for the specified period of interest. The output file is in ASCII column format:

```<year> <month> <day> <prec> <evap> ...```

All dates are in the form MMDDYYYY (MM = month, DD = day, YYYY = year) and computations start at midnight the morning of that day.

This program is currently configured to read files with 3 soil layers.
