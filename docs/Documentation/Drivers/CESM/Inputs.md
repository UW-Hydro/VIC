# CESM Driver Inputs

To run VIC, several sets of input data are necessary:

*   Input Configuration File: This file lists the individual input files used by the CESM driver. The heading (e.g. `wr50a`) describes the domain grid used. The individual input files described are described below.

```
[wr50a]
vic_domain = ${DIN_LOC_ROOT}/share/domains/domain.lnd.wr50a_ar9v4.100920.nc
vic_params = ${VIC5TESTDIR}/vic_params_wr50a_vic5.0.dev.nc
vic_global_param = vic.globalconfig.txt
vic_constants = vic.parameters.txt
```

*   [Global Parameter File](../Image/GlobalParam.md): This is the main input file for VIC. It points VIC to the locations of the other input/output files and sets parameters that govern the simulation (e.g., start/end dates, modes of operation).  TODO: Write custom Global Parameter page for CESM driver.
*   [Parameters File](../Image/Params.md): Spatially distributed parameters describing the land surface. The CESM driver uses the same Parameters file format as the Image Driver.

And a few more are optional:
*   [Constants File](../../Constants.md): Model parameters that are constant in time and space. The CESM driver uses the same constants file format as the Image Driver.
*   [Initial State File](../Image/StateFile.md): Moisture storages (soil moisture, snow pack, etc), energy storages (soil temperatures, etc) and other information describing the current state of the system. A state file saved from a previous VIC simulation may be used as the initial state for another run. The CESM driver uses the same state file format as the Image Driver.
