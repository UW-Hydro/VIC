# Code

Note: The VIC model and the routing model have been developed for use on LINUX and UNIX platforms. To use VIC and/or the routing model on a Windows platform, we suggest downloading a free UNIX emulator such as [Cygwin](http://www.cygwin.com) and compiling/running these models within this emulator. Please do not ask us for help with installing or using Cygwin.

## VIC Model Source Code

The VIC model is written in C and typically compiled with the [GNU gcc](http://gcc.gnu.org/) compiler.

**Access:** The VIC model source code is archived in Git and is publicly available through [GitHub](https://github.com). To access the source code, visit [GitHub](https://github.com), create an account, and visit [github.com/UW-Hydro/VIC](https://github.com/UW-Hydro/VIC).

## Routing Model Source Code

The routing model is written in Fortran 77 and typically compiled with the [GNU g77](https://gcc.gnu.org/fortran/) compiler (not to be confused with GNU's gfortran compiler, which may also work but hasn't been tested here in our lab).

*   [Route 1.0](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Source_Code/route_code_1.0.tgz)

## Preparing VIC Model Inputs

Follow these instructions if you don't want to use our parameters.

*   Basin Delineation
    *   Arc/Info Utility Basin Preparation Programs: [basin_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/basin_prep.tgz)
*   Meteorological Data
    *   Programs for gridding meteorological data from the National Climatic Data Center: [GRID_2000.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/GRID_2000.TAR.gz)
    *   Programs for processing "Canada Daily Climate Data": [CANADA_MET.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/CANADA_MET.TAR.gz). Daily meteorological data for Canada is available on CDs from Environment Canada, entitled "Canada Daily Climate Data." This includes a README file describing the steps to downloading and processing the data for input to VIC.
    *   Additional Processing MethodPrograms: [met_prep_other.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/met_prep_other.tgz)
    *   Temporal Disaggregation of Precipitation Data
        *   Tools for TemporalDisaggregation of Precipitation Data: [temporal_disagg_precip.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_disagg_precip.tgz)
        *   Lookup tables for thecontinental U.S. (by season): [US_LOOKUP_TABLES.TAR.gz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/US_LOOKUP_TABLES.TAR.gz)
*   Parameter Files
    *   Vegetation Parameter Preparation: [veg_param_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/veg_param_prep.tgz)
    *   Soil Parameter Preparation: [soil_param_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/soil_param_prep.tgz)
    *   Snowband Preparation: [snowband_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/snowband_prep.tgz)
    *   LDAS File Preparation: [LDAS_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/LDAS_prep.tgz)

## Preparing Routing Model Inputs

*   Routing Model Input Preparation: [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz)
*   Jennifer Adam's Routing Model Preparation: [rout_prep_adam.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep_adam.tgz)

## Post-processing

*   Temporal Aggregation: [temporal_agg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_agg.tgz)
*   Spatial Aggregation: [spatial_agg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/spatial_agg.tgz)
*   Spatial Disaggregation: [spatial_disagg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/spatial_disagg.tgz)
*   Tool for converting between standardVIC input/output and NetCDF: [flux2nc.py](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/flux2nc.py)

## Calibration

*   Calibration Programs: [calibrate_other.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/calibrate_other.tgz) for method described on [Other Calibration Methods](../Documentation/CalibrateMethodOther.md) website.

## Plotting Scripts

*   Splus (or R) scripts for comparing 2 VIC runs (at one grid cell): [R_plot_scripts.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Plotting_Scripts/R_plot_scripts.tgz)
*   C-shell/perl/gmt script for plotting hydrographs: [plot_flow_STEHE.scr](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Plotting_Scripts/plot_flow_STEHE.scr)

## Scripts For Converting File Formats Between VIC Versions

*   Perl script for converting state files created by older versions of VIC to the format needed by the latest version of VIC: [convert_state_file.pl](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/convert_state_file.pl)
