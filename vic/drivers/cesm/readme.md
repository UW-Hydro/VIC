CESM Driver for VIC Model
================

This driver was developed for coupling VIC within the Regional Arctic System
Model (RASM).

# Directory Structure

    | - cesm
    | --- bld
    | --- cpl_mct
    | --- include
    | --- src

# Building the Model

The CESM driver for VIC can be built in two ways.

1. The `Makefile` is in the root driver directory is configured to build the CESM driver as a shared object. This makefile is provided for testing purposes only.

1. The `bld/` directory contains the build scripts used by the CESM  `$case.build` script. To build the model in this configuration, follow these temporary steps on Garnet, Spirit, or Lightning:

  ```bash
  # Get vic5 input data
  cd $HOME
  export VIC5TESTDIR=$HOME/vic5_inputdata
  cd $VIC5TESTDIR
  /usr/bin/curl -O ftp://ftp.hydro.washington.edu/pub/jhamman/rasm_vic5_test_data/vic_params_wr50a_vic5.0.dev.nc

  # Get vic5 branch in rasm repo
  cd $HOME
  # Private repository, not publicly available
  svn co https://svn.nps.edu/repos/racm/rasm/branches/vic5/ rasm_vic5

  # move vic4 out of the way
  cd rasm_vic5/models/lnd/
  mv vic vic4

  # get vic5 rasm repo
  git clone git@github.com:UW-Hydro/VIC.git vic
  cd vic
  git checkout develop

  # follow typical steps to build RASM
  # NOTE: only set DEBUG flag to TRUE if running RI compset
  # (it does not work with WRF)
  cd $HOME/rasm_vic5/scripts
  today=$(date +'%Y%m%d')
  compset=RI # adjust for compset
  mach=spirit_intel # adjust for machine
  case_name=vic5.${compset}.test.${today}a
  create_newcase -case ${case_name} -res w5a_a94 -compset ${compset} -mach ${mach}
  cd ${case_name}
  ./cesm_setup
  ./xmlchange -file env_build.xml -id DEBUG -val TRUE
  ./${case_name}.build
  qsub ${case_name}.run
  ```

  ** Supported Machines **
  - [x] Thunder
  - [x] Lightning
  - [x] Garnet
  - [ ] Copper *(Not currently supported by RASM)*
