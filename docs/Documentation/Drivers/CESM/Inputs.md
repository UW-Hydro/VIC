# CESM Driver Inputs

!!! Warning
    The VIC CESM Driver is experimental and is still under development. No documentation or support is provided for running CESM in any way.

To run CESM or RASM with VIC, several sets of input data are necessary:

*   Input Configuration File: This file lists the individual input files used by the CESM driver. The heading (e.g. `wr50a`) describes the domain grid used and is used by `build_vic_namelist` via the `LND_GRID` environment variable set by the CESM build system. The individual input files described are described below. Currently, the `vic_global_param` and `vic_constants` are included in the `VIC/vic/drivers/cesm/bld/` directory. Additional default configuration files may be added in the future.

```
    [wr50a]
    vic_domain = ${DIN_LOC_ROOT}/share/domains/domain.lnd.wr50a_ar9v4.100920.nc
    vic_params = ${DIN_LOC_ROOT}/lnd/vic/vic_params_wr50a_vic5.0.dev.nc
    vic_global_param = vic.globalconfig.txt
    vic_constants = vic.parameters.txt
```

*   [Global Parameter File](../Image/GlobalParam.md): This is the main input file for VIC. It points VIC to the locations of the other input/output files and sets parameters that govern the simulation (e.g. input parameter files, modes of operation, and output variables).

    The CESM driver's global parameter file differs from the [Image Driver's global parameter file](../Image/GlobalParam.md) in the following ways:

    1.  Options specific to the coupling with CPL7 are not specified (e.g. `MODEL_STEPS_PER_DAY`).
    2.  Some options will be filled in by the CESM build system (e.g. `build_vic_namelist`).

*   [Parameters File](../Image/Params.md): Spatially distributed parameters describing the land surface. The CESM driver uses the same parameter file format as the Image Driver.

And a few more are optional:

*   [Constants File](../../Constants.md): Model parameters that are constant in time and space. The CESM driver uses the same constants file format as the Image Driver.
*   [Initial State File](../Image/StateFile.md): Moisture storages (soil moisture, snow pack, etc), energy storages (soil temperatures, etc) and other information describing the current state of the system. A state file saved from a previous VIC simulation may be used as the initial state for another run. The CESM driver uses the same state file format as the Image Driver.
