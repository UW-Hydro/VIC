# VIC Model Inputs

VIC input files may be constructed using various programs and datasets. Below are general descriptions of each input file along with links to its structure and the methods that may be used to build it. For a comprehensive example, including an example of each input file, see the Stehekin Basin file in the Sample Data Sets section of the [downloads page](../../../Datasets/Datasets.md).

To run VIC, several sets of input data are necessary:

*   [Global Parameter File](GlobalParam.md): This is the main input file for VIC. It points VIC to the locations of the other input/output files and sets parameters that govern the simulation (e.g., start/end dates, modes of operation).
*   [Meteorological Forcing Files](ForcingData.md): Gridded, sub-daily timeseries of meteorological variables as inputs.
*   [Parameters File](Params.md): Spatially distributed parameters describing the land surface.
*   [Domain File](Domain.md): Domain information of VIC run.

And a few more are optional:

*   [Constants File](../../Constants.md): Model parameters that are constant in time and space.
*   [Initial State File](StateFile.md): Moisture storages (soil moisture, snow pack, etc), energy storages (soil temperatures, etc) and other information describing the current state of the system. A state file saved from a previous VIC simulation may be used as the initial state for another run.
*   [Lake/Wetland Parameter File](LakeParam.md): File containing lake model parameters. By default, VIC does not simulate lakes or other impoundment of surface water.
