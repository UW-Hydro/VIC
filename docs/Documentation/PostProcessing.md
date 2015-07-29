# Post-Processing VIC Model Output

VIC model output may require post-model run processing in order to make the output data more useful.

## Temporal Aggregation: Reducing Hourly Output to Daily

When VIC is run at a sub-daily time step (this is typical for FULL_ENERGY = TRUE or FROZEN_SOIL = TRUE), it will by default write its outputs at the same sub-daily time step. This can produce large volumes of data, much of which is not needed, even when the output is written in binary format. There are two options for solving this problem:

1.  VIC can aggregate the output to a daily time step before writing it to output files if the user specifies OUT_STEP = 24 in the [global parameter file](GlobalParam.md).
2.  If you wish to see the outputs at both sub-daily and daily time steps, you will need to have VIC write its outputs at sub-daily time steps and then you will need to reduce the data to daily as a post-processing script. For more information on this, [click here](TemporalAggregation.md).

## Spatial Aggregation

It is often desirable to aggregate model inputs or outputs to a coarser resolution. We have scripts for spatial aggregation of VIC input and output files, and routing input files, [here](SpatialAggregation.md).

## Temporal Disaggregation

VIC internally disaggregates daily meteorological forcings to sub-daily intervals using the Thornton and Running algorithm among others (reference to be posted soon). In this way, VIC can be used as a temporal disaggregator of forcings - simply recompile VIC with OUTPUT_FORCE set to TRUE in the [global parameter file](GlobalParam.md), set the TIMESTEP to a sub-daily interval in your [global parameter file](GlobalParam.md), and make sure there is no output file format specification in your [global parameter file](GlobalParam.md). By default, VIC will bypass its simulation and simply produce output files consisting of sub-daily PRECIP, AIR_TEMP, SHORTWAVE, LONGWAVE, PRESSURE, VP, DENSITY, and WIND. NOTE: the sub-daily precipitation is assumed to fall at the same average rate throughout the day.

To disaggregate precipitation so as to be non-uniform throughout the day, download the file [temporal_disagg_precip.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/temporal_disagg_precip.tgz), also available under "Preparing VIC Model Inputs" section on the [download page](../SourceCode/Download.md).

## Spatial Disaggregation: Sub-Sampling Data to a Finer Resolution

Scripts to subsample meteorological data to finer resolution are in the file [spatial_disagg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/spatial_disagg.tgz), also available under the "Post-processing" section on the [download page](../SourceCode/Download.md).

## VIC NetCDF Tools

NetCDF is a machine-independent format for representing scientific data. It is widely used among the atmospheric sciences community and is quickly becoming a standard format for use in comparing outputs between different models.

For more information on NetCDF Tools, [click here](VICNetCDF.md).
