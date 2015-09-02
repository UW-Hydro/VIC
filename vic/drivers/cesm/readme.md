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

1. The `Makefile` is in the root driver directory is configured to build the CESM driver as a shared object. This makefile is provided for testing purposes only.

1. The `bld/` directory contains the build scripts used by the CESM  `$case.build` script. To build the model in this configuration, follow these temporary steps on Garnet, Spirit, or Lightning:

  ```bash
  svn checkout https://svn.nps.edu/repos/racm/rasm/branches/vic5/ rasm_vic5
  cd rasm_vic5/models/lnd/
  mv vic/ vic4
  git clone git@github.com:UW-Hydro/VIC.git vic
  cd vic
  git checkout driver/cesm
  cd ../../../../scripts
  # follow typical steps to build RASM
  ```
