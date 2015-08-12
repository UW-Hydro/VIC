# Spatially Aggregating Data Files

This page describes a set of programs designed to spatially aggregate VIC model data and parameter files by factors of 2 (e.g. 1/8 -> 1/4 -> 1/2 -> 1 -> 2 degrees). Programs included with this page were developed by Bernt Viggo Matheussen and Keith Cherkauer.

All programs discussed on this page are included in the file [spatial_agg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/spatial_agg.tgz), also available under the "Post-processing" section of the [download page](../SourceCode/Code.md).

There are three programs that need to be run for each halving of the forcing data resolution:

*   [Meteorological Forcing Files](#A1)
    *   [make_new_fraction.c](#A1.1) - which creates a fraction file for the new resolution based on the routing fraction file of the current resolution
    *   [latlong.c](#A1.2) - which generates a three column file with latitude, longitude and fractional coverage for the current resolution.
    *   [aggregate_metdata.c and aggregate_binary_metdata.c](#A1.3) - which aggregates the forcing data.
*   [Routing Model Files](#A2)
    *   [make_accumulation.c](#A2.1) - computes the accumulation of runoff for each model pixel.
    *   [make_new_direction.c](#A2.2) - works out the flow direction file for the lower resolution file.
*   [Soil Parameter Files](#A3)
    *   [aggregate_soil_data.c](#A3.1) - aggregates soil parameters in ARCINFO grid format.
*   [Vegetation Parameter Files](#A4)
    *   [aggregate_veg_param.c](#A4.1) - aggregates the vegetation parameter file.

* * *

## Meteorological Forcing Files

## make_new_fraction

```Usage: make_new_fraction < fraction file >```

`< fraction file >` is the routing fraction file for the current model resolution WITH AN ARC/INFO HEADER.
Output is to stdout.

This program reduces the resolution of the provided routing fraction file by half (e.g. 1/8 -> 1/4 degree). The new fraction file is sent to stdout, and can be redirected to a new file.

WARNING: This program requires that the fraction file have a full ARC/INFO header. Older versions of the routing code did not include the header as part of the file, so check your fraction file and add the header from the direction file if necessary.

## latlong

Usage: ```latlong < inputfile >```

This program reads an ARC/INFO grid file and returns an 3 column ASCII file (latitude, longitude, and fraction). Run it on the fraction file for both the current resolution and the new resolution. Both column files are needed to run aggregate_forcing_data. Output from the program is to stdout, and can be redirected to a file.

## aggregate_metdata and aggregate_binary_metdata

Usage:  ```aggregate_metdata < old fraction xyz > < old fraction data dir > < new fraction xyz > < new fraction data dir > < old resolution > < number of parameters >```

This program aggregates the gridded metdata files by a factor of 2 (e.g. 1/8 degree is aggregated to 1/4 degree). It uses the three column ASCII files created by [latlong](#A2) to determine the file names at the old and new resolutions. The program can be run on data files with any number of parameters, as long as each is in a separate column, and is continuous through the data file.

The following variables need to be defined in the source code before compiling:

*   MAX - the maximum number of grid cells at any resolution,
*   STEPS - the number of records in the data file (e.g. for daily data STEPS = number of days),
*   MAX_PARAM - the maximum number of parameters in the data file,
*   OLD_DECIMAL_PLACES - the number of decimal places used for the latitude and longitude of the old resolution files, and
*   NEW_DECIMAL_PLACES - the number of decimal places used for the latitude and longitude of the new resolution files (the simpliest approach is to set this the same as OLD_DECIMAL_PLACES).

* * *

## Routing Model Files

## make_accumulation

Usage: ```make_accumulation < direction file >```
        This program computes the flow accumulation for each grid cell
        based on the given routing direction file.  The accumulation
        file is output to stdout.

This program reads the current routing model direction file and computes the flow accumulation for each grid cell. The resulting accumulation grid file is output to stdout. Accumulation at the basin outlet should be equal to the number of grid cells within the basin.

## make_new_direction

Usage: ```make_new_direction < direction file > < accumulation file >```

This program builds a new routing direction file at twice the
resolution (e.g. 1/8 to 1/4) of the provided direction file.  
The accumulation file for the current direction file must also
be provided.  The new direction file is output to stdout.

This program reads the high resolution direction and accumulation files, and determines the most appropriate drainage direction for each grid cell at the coarser resolution. The new drainage direction is based on the direction of flow from the cell with the highest accumulation (which must also leave the lower resolution cell). If the resulting direction is a diagonal it is checked to see if it should actually be a diagonal, or if it really flow into the neighboring low resolution cell. The new direction file is output to stdout.

* * *

## Soil Parameter Files

## aggregate_soil_data

Usage: ```aggregate_soil_data <high res="" grid=""> <low res="" mask=""> <low res="" output="" gird="">```

<high res="" grid=""> is the high resolution ARC/INFO grid.
<low res="" mask=""> is an ARC/INFO grid mask file at the lower resolution.
This program aggregates the high resolution ARCINFO file to
the resolution of the low resolution mask file (only aggregates
by a factor of 2, e.g. 1/8 to 1/4, 1/4 to 1/2) and outputs a
new ARCINFO grid file to <low res="" output="" grid="">.</low> </low></high></low></low></high>

This program aggregates ARCINFO ASCII grid files (typically containing soil parameters) by a factor of 2 (e.g. 1/8 to 1/4, 1/4 to 1/2). Either the direction file or the fraction file generated above can be used as the low resolution mask file.

* * *

## **Vegetation Parameter Files**

## aggregate_veg_param

Usage: ```aggregate_veg_param <high res="" vegparam="" file=""> <high res="" cellnum="" file=""> <low res="" vegparam="" file=""> <low res="" cellnum="" file=""> <# root zones>```

This program aggregates the given vegetation parameter file by a factor of 2 (e.g. 1/8 to 1/4, 1/4 to 1/2, etc.).
<high res="" vegparam="" file=""> is the vegetation parameter file for the high resolution model simulation.
<high res="" cellnum="" file=""> is an ARCINFO ASCII grid cell number file for the high resolution model basin.
<low res="" vegparam="" file=""> is the vegetation parameter file for the low resolution (output) model simulation.
<low res="" cellnum="" file="">is an ARCINFO ASCII grid cell number file for the low resolution model basin.</low> </low></high></high></low></low></high></high>

This program aggregates...
