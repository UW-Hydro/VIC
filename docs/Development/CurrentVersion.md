# Current Version of VIC Model

## Current Version: 4.2 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22307.svg)](http://dx.doi.org/10.5281/zenodo.22307)

The current version of VIC is 4.2.  This is a minor release in the VIC 4 development track.  The previous release was VIC 4.1.2.m.

Here we provide a brief description of the features of 4.2 and differences between 4.2 and previous versions.

## Summary of new features

VIC release 4.2 contains the following new features relative to 4.1.2:

*   Ability to simulation partial vegetation cover
*   Ability to input non-climatological time-varying vegetation parameters
*   Simulation of photosynthesis
*   Simulation of soil carbon storage and fluxes
*   New options for handling frozen soils
*   Removal of some obsolete options and general clean-up

## Known Issues in Current Release

For a list of known issues and their fixes (in bug-fix updates), visit the VIC GitHub [Issues](https://github.com/UW-Hydro/VIC/issues) page.

## Differences between 4.2 and previous versions

For more information about VIC 4.2 features and the differences VIC 4.2 has with earlier versions, [click here](VersionSummaries.md).

## Backwards-compatibility

Existing model setup from VIC 4.1.2 will still run with VIC 4.2.  The only exception to this may be when using an old state file to restart from a previous simulation

Most of the changes in VIC 4.2 are new features that are turned off by default.  Therefore, similar, but not exact, results should be obtained by both VIC 4.1.2 and VIC 4.2.

## Earlier Versions

For information on archival versions of VIC, [click here](ArchivedVersions.md).

## Version Check

To check which release of VIC you are running:

Type `"vicNl -v"`

# Current Version of Routing Model

## Current Version: 1.0

This is the routing model of [Lohmann, et al. (1996)](../Documentation/References.md), as implemented at UW, as of 2001.

## How to Obtain the Code

For information on obtaining the routing model source code, [click here](../SourceCode/Code.md).
