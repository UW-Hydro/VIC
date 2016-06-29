# Running the VIC Image Driver

## Dependencies:
The Image Driver's has three dependencies:

1. A C compiler.  We routinely test VIC using the following compilers:

    - [GNU](https://gcc.gnu.org/) (`gcc` version 4+)
    - [Clang](http://clang.llvm.org/) (`clang` version 3+)

    VIC has also been compiled using these compilers:

    - [Intel](https://software.intel.com/en-us/c-compilers) (`icc`)
    - [PGI](http://www.pgroup.com/) (`pgcc`)

2. MPI.  We have tested VIC with the following MPI implementations:
    - [Open MPI](http://www.open-mpi.org/) (version 1.5.4+)
    - [MPICH](http://www.mpich.org/) (version 1.2+)

3.  [netCDF4](http://www.unidata.ucar.edu/software/netcdf/)

## Compiling
In most cases, you will need to edit the `NETCDF_PATH` and `MPI_PATH` variables in the `Makefile`.

If you want to use a compiler other than `gcc`, either edit the Makefile or set the `CC` environment variable, e.g.

        export CC=icc

- Change directory, `cd`, to the "Image Driver" source code directory and type `make`

        cd vic/drivers/image
        make

- If this completes without errors, you will now see a file called `vic_image.exe` in this directory. `vic_image.exe` is the executable file for the model.

## Run VIC

At the command prompt, type:

`vic_image.exe -g global_parameter_filename`

where `global_parameter_filename` = name of the global parameter file corresponding to your project.

## Other Command Line Options

VIC has a few other command line options:

- `vic_image.exe -v`: says which version of VIC this is
- `vic_image.exe -h`: prints a list of all the VIC command-line options
- `vic_image.exe -o`: prints a list of all of the current compile-time settings in this executable; to change these settings, you must edit the appropriate header files (e.g. `vic_def.h` or `vic_driver_shared.h`) and recompile using `make full`.
