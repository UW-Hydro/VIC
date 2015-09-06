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
  # Get vic5 input data
  cd $HOME
  export VIC5TESTDIR=$HOME/vic5_inputdata
  cd $VIC5TESTDIR
  /usr/bin/curl --header 'Host: dl-web.dropbox.com' --header 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:40.0) Gecko/20100101 Firefox/40.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --header 'Referer: https://www.dropbox.com' --header 'Cookie: locale=en; t=_m_3OjJVG0C4_KFK1Cz1Vd_1; blid=AADO2Ivt6-gp5e0WvK8GH4c42fv31BSECHgI_d_XxhH2TQ; bjar=W3sidWlkIjogNDk1MDg2MjMsICJzZXNzX2lkIjogNzA0NjA1ODE4OTQxNjA3OTY3NzI2NjcxOTg1OTIwOTcxODcwMDEsICJleHBpcmVzIjogMTQ0MzcxNDY3OSwgInRlYW1faWQiOiAiIiwgInJvbGUiOiAicGVyc29uYWwifV0%3D' --header 'Connection: keep-alive' 'https://dl-web.dropbox.com/get/UW_RASM_shared/calibration/vic_params_wr50a_vic5.0.dev.nc?_subject_uid=49508623&w=AACXfsAfXEE1Fvjp04FKoS16Wu4jLEi-uT8Ad6dau8yYpA&dl=1' -o 'vic_params_wr50a_vic5.0.dev.nc' -L

  # Get vic5 branch in rasm repo
  cd $HOME
  svn co https://svn.nps.edu/repos/racm/rasm/branches/vic5/ rasm_vic5

  # move vic4 out of the way
  cd rasm_vic5/models/lnd/
  mv vic vic4

  # get vic5 rasm repo
  git clone git@github.com:UW-Hydro/VIC.git vic
  cd vic
  git checkout driver/cesm

  # follow typical steps to build RASM
  cd $HOME/rasm_vic5/scripts
  ./create_case ...
  ```

  ** Supported Machines **
  - [x] Spirit
  - [x] Lightning
  - [ ] Garnet
  - [ ] Copper *(Not currently supported by RASM)*
