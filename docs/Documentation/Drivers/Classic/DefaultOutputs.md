# Default Output Files

If the user does not give explicit output file instructions in the global parameter file, then VIC will produce a set of default output files, one per grid cell. The exact set of output files, and which variables they contain, is determined by model settings, as follows:

## Flux Files

The primary output file type for the VIC model is the flux file, which contains information about moisture and energy fluxes for each time step. These output files are based on the model output files used for the PILPS-2C project. The number of variables in this file depends on the values of `FULL_ENERGY` and `FROZEN_SOIL` in the [global parameter file](GlobalParam.md). When either `FULL_ENERGY` or `FROZEN_SOIL` are true, flux files are always written, regardless of the mode of operation. Flux files begin with the prefix `fluxes_`.

For more information on the Flux Output Files, [click here](FluxOutputFiles.md).

## Snow Files

The snow file contains information about the snowpack, averaged across all elevation bands and vegetation tiles. The set of variables in this file depends on the values of `FULL_ENERGY` and `FROZEN_SOIL` in the [global parameter file](GlobalParam.md). Snow files are always written, regardless of the mode of operation. Snow files begin with the prefix `snow_`.

For more information on the Snow Output File, [click here](SnowOutputFile.md).

## Frozen Soil Files

When the model is run with `FROZEN_SOIL` set to TRUE in the [global parameter file](GlobalParam.md), a third output file is produced which contains soil thermal output parameters. Frozen soil files begin with the prefix `fdepth_`.

For more information on the Frozen Soil Output File, [click here](FrozenSoilOutputFile.md).

## Snow Band Files

When the model is run with snow elevation bands, and default outputs are specified in the [global parameter file](GlobalParam.md), snow pack information for each elevation band will be output to files in the results directory with the prefix `snow_band_`. Energy fluxes are output only for the full energy balance model, so there are file descriptions for both full energy and water balance model output files.

For more information on the Snow Band Output Files, [click here](SnowBandOutputFiles.md).

## Lake Files

The lake file contains information about the lake fraction of each grid cell. Lake files are only written when LAKES is set equal to a valid lake parameter file, in the [global parameter file](GlobalParam.md). Lake files begin with the prefix `lake_`.

For more information on the Lake Output File, [click here](LakeOutputFile.md).
