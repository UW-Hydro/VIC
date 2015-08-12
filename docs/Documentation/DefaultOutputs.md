# Default Output Files

If the user does not give explicit output file instructions in the global parameter file, then VIC will produce a set of default output files, one per grid cell. The exact set of output files, and which variables they contain, is determined by model settings, as follows:

*   If OUTPUT_FORCE is TRUE in the [global parameter file](GlobalParam.md), VIC will produce 1 output file per grid cell, containing: PREC, AIR_TEMP, SHORTWAVE, LONGWAVE, DENSITY, PRESSURE, VP, and WIND. Units of these variables will depend on whether ALMA_OUTPUT is set to TRUE or FALSE in the [global parameter file](GlobalParam.md). These files will be named "full_data_lat_lon", where "lat" and "lon" are the latitude and longitude of the grid cell.
*   If OUTPUT_FORCE is FALSE in the [global parameter file](GlobalParam.md), the set of output files, and the variables they contain, will depend on which operation mode VIC is being run in (e.g. FULL_ENERGY, FROZEN_SOIL, etc. in the [global parameter file](GlobalParam.md)).

VIC can also be explicitly instructed to output these default files by pasting the contents of the file output.TRADITIONAL.410.template, which is included with the [VIC source code](../SourceCode/Code.md), into the bottom of your [global parameter file](GlobalParam.md). Of course, this is not necessary, since omitting any output file format information will cause VIC to make these files by default. But if you wish to have output that is similar to the default output files, but with fewer or additional variables, it is easiest to paste the contents into the global parameter file and then edit them to your liking. NOTE: to have output files exactly similar to those of 4.0.6 and earlier, use output.TRADITIONAL.template, instead of output.TRADITIONAL.410.template.

## Flux Files

The primary output file type for the VIC-NL model is the flux file, which contains information about moisture and energy fluxes for each time step. These output files are based on the model output files used for the PILPS-2C project. The number of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the [global parameter file](GlobalParam.md); when either FULL_ENERGY or FROZEN_SOIL are true, the file will contain the same number of columns as found in the PILPS-2C files for comparison purposes. Flux files are always written, regardless of the mode of operation. Flux files begin with the prefix "fluxes_".

For more information on the Flux Output Files, [click here](FluxOutputFiles.md).

## Snow Files

The snow file contains information about the snowpack, averaged across all elevation bands and vegetation tiles. The set of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the [global parameter file](GlobalParam.md). Snow files are always written, regardless of the mode of operation. Snow files begin with the prefix "snow_".

For more information on the Snow Output File, [click here](SnowOutputFile.md).

## Frozen Soil Files

When the model is run with FROZEN_SOIL set to TRUE in the [global parameter file](GlobalParam.md), a third output file is produced which contains soil thermal output parameters. Frozen soil files begin with the prefix "fdepth_".

For more information on the Frozen Soil Output File, [click here](FrozenSoilOutputFile.md).

## Snow Band Files

When the model is run with snow elevation bands, and PRT_SNOW_BAND is turned on in the [global parameter file](GlobalParam.md), snow pack information for each elevation band will be output to files in the results directory with the prefix "snow_band_". Energy fluxes are output only for the full energy balance model, so there are file descriptions for both full energy and water balance model output files.

For more information on the Snow Band Output Files, [click here](SnowBandOutputFiles.md).

## Lake Files

The lake file contains information about the lake fraction of the grid cell. Lake files are only written when LAKES is set equal to a valid lake parameter file, in the [global parameter file](GlobalParam.md). Lake files begin with the prefix "lake_".

For more information on the Lake Output File, [click here](LakeOutputFile.md).
