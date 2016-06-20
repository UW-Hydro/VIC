# Default Output Files - Image Driver

If the user does not give explicit output file instructions in the global parameter file, then VIC will produce a set of default output files. The exact set of output files, and which variables they contain, is determined by model settings, as follows:

## Flux Files

The primary output file type for the VIC-NL model is the flux file, which contains information about moisture and energy fluxes for each time step. These output files are based on the model output files used for the PILPS-2C project. The number of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the [global parameter file](GlobalParam.md); when either FULL_ENERGY or FROZEN_SOIL are true, the file will contain the same number of variables as found in the PILPS-2C files for comparison purposes. Flux file name is `fluxes.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

For more information on the Default Flux Output Files, [click here](FluxOutputFiles.md).

## Snow Files

The snow file contains information about the snowpack, averaged across all elevation bands and vegetation tiles. The set of variables in this file depends on the values of FULL_ENERGY and FROZEN_SOIL in the [global parameter file](GlobalParam.md). Snow files are always written, regardless of the mode of operation. Snow file name is `snow.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

For more information on the Snow Output File, [click here](SnowOutputFile.md).

## Frozen Soil Files

When the model is run with FROZEN_SOIL set to TRUE in the [global parameter file](GlobalParam.md), a third output file is produced which contains soil thermal output parameters. Frozen soil file name is `fdepth.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

For more information on the Frozen Soil Output File, [click here](FrozenSoilOutputFile.md).

## Lake Files

The lake file contains information about the lake fraction of the grid cell. Lake files are only written when LAKES is set equal to a valid lake parameter file, in the [global parameter file](GlobalParam.md). Lake file name is `lake.**.nc`, where `**` is determined by AGGFREQ variable specified in the [global parameter file](GlobalParam.md).

For more information on the Lake Output File, [click here](LakeOutputFile.md).
