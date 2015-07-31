# Gridding Meteorological Data from the National Climatic Data Center

This document provides a brief explanation of the methodology involved in developing a gridded meteorological data set (daily precipitation, maximum and minimum temperatures, and daily average wind speed) for the hydrology VIC model using data over the U.S. The process includes two general steps, based on two data sources: 1) taking raw data from the EarthInfo National Climate Data Center (NCDC) CDs and converting them into VIC input files, and 2) appending more recent data downloaded from the [NCDC On-Line Web Site](http://hurricane.ncdc.noaa.gov/CDO/cdo). The core of the gridding process is the interpolation routine called SYMAP (Shepard, D.S., Computer Mapping: the SYMAP Interpolation Algorithm, In: _Spatial Statistics and Models_, Gaile and Willmott, eds., 1984). For precipitation, all the interpolated data are scaled to match long-term monthly means from the [PRISM](http://www.prism.oregonstate.edu/) monthly precipitation dataset. For 10 meter daily wind data, gridded data are obtained from the [NCEP/NCAR Reanalysis](http://www.esrl.noaa.gov/psd/data/gridded/reanalysis/), and linearly interpolated to the VIC grid resolution.

The preprocessing steps are performed using Unix shell scripts and programs written both in C and in Fortran. A basic knowledge of Unix, C and Fortran is presumed, although the programs should not require alteration.

The source code of the programs used, and documentation (in MS Word and PDF format) descriping in detail the necessary steps to build a VIC meteorological input dataset, is downloadable from this website.

Download the file [GRID_2000.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/GRID_2000.TAR.gz), also available on the [download page](../SourceCode/Code.md) under the section "Programs for gridding meteorological data from the National Climatic Data Center". This file contains all the necessary source code and documentation for gridding meteorological data.

The main program involved in gridding the met data is called "regrid". As of 1/10/2001, the regrid program has been revised. The new version, included in GRID_2000.TAR.gz, can be compiled on FreeBSD machines. The previous version used dynamic memory allocation, but could only be compiled and run on HP-UX systems.

The file `GRID_2000.TAR.gz` is compressed together with TAR and gzip.

To uncompress use: `gzip -d GRID_2000.TAR.gz`

To extract the files use: `tar -xfv GRID_2000.TAR`

* * *

There is also daily meteorological data for Canada on CDs available from Environment Canada, entitled "Canada Daily Climate Data." Programs to process the downloaded data are stored in the file [CANADA_MET.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/CANADA_MET.TAR.gz), also available on the [download page](../SourceCode/Code.md) under the section "Programs for gridding meteorological data from the National Climatic Data Center". This includes a README file describing the steps to downloading and processing the data for input to VIC.

The University of Washington takes no responsibility for any damage or errors that these programs contain or may produce.
