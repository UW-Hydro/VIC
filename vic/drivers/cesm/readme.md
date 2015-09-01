CESM Driver for VIC Model
================

This driver was developed for coupling VIC within the Regional Arctic System
Model (RASM).

# Directory Structure

    | - cesm
    | --- bld
    | --- cpl_esmf
    | --- cpl_mct
    | --- include
    | --- src

# Building the Model

The CESM driver for VIC can be built in two ways.

1. The `bld/` directory contains the build scripts used by the CESM  `$case.build` script.

1. The `Makefile` is in the root driver directory is configured to build the CESM driver as a shared object. This makefile is provided for testing purposes only.
