VIC CESM Driver
================

!!! Warning
    The VIC CESM Driver is experimental and is still under development. No documentation or support is provided for running CESM in any way.

The CESM driver was developed for coupling VIC within the [Regional Arctic System Model (RASM)](http://uw-hydro.github.io/current_project/rasm). Although it has only been tested as part of RASM, it should be possible to apply this coupling to CESM 1.

## Dependencies:
In addition to the dependencies listed in the [image driver documentation](../Image/RunVIC/#dependencies), the CESM driver has the following dependencies:

1. A Fortran compiler with Fortran 2003 or later.  We routinely test VIC using the following compilers:

    - [GNU](https://gcc.gnu.org/fortran/) (`gfortran` version 4.9+)
    - [PGI](http://clang.llvm.org/) (`pgfortran` version 15.5+)
    - [Intel](http://clang.llvm.org/) (`ifort` version 16.0+)

## Building the Model

The CESM driver for VIC can be built in two ways.

1. The `Makefile` is in the root driver directory is configured to build the CESM driver as a shared object. This makefile is provided for testing purposes only.

2. Using the CESM build system. Users should refer to the [CESM documentation](http://www.cesm.ucar.edu/models/cesm1.0/cesm/) for more information on how to build cases.

!!! Note
    Note that the CESM driver is built as a library and is intended to be linked to the CESM coupler by the CESM build system, it is not a stand-alone executable.

## Configuring VIC as the land surface component in CESM.
Most of the inputs and options for the CESM driver are identical to the Image driver. See the [Image Driver Documentation](../Image/ImageDriver.md) for more details.

1.  [VIC CESM Driver Input Files](Inputs.md)
