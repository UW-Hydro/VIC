# Running the VIC Image Driver

## Dependencies:
The Image Driver has three dependencies:

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

!!! Note
    Compiling the Image Driver may also be done with [OpenMP](http://www.openmp.org/). Nearly all modern C compilers include the [OpenMP standard](http://www.openmp.org/resources/openmp-compilers/) and users will need to ensure that the makefile has the appropriate compiler flag (usually `-fopenmp`). See the discussion below for how to control OpenMP parallelization.

## Compiling
In most cases, you will need to edit the `NETCDF_PATH` and `MPI_PATH` variables in the `Makefile`.

If you want to use a compiler other than `mpicc`, either edit the Makefile or set the `MPICC` environment variable, e.g.

        MPICC=/path/to/mpi_c_compiler

The flags and libraries required to compile VIC with netCDF are automatically determined in the `Makefile`.  They can be overwritten by setting the following two environment variables.  These variables can be determined by running `nc-config --all`.

        NC_LIBS="-L/path/to/libs ..."
        NC_CFLAGS="-I/path/to/includes -your_c_flags ..."

In some versions of the MPI library (e.g. OPEN-MPI with Intel), you may also need to set the environment variable `MX_RCACHE=2` prior to compiling.

To enable the river routing extension, you must set the ROUT option. This includes setting extension in the Makefile and adding routing-specific input parameter file. For more information on how to enable the routing extension, see the [routing extension documentation](Routing.md).

- Change directory, `cd`, to the "Image Driver" source code directory and type `make`

        cd vic/drivers/image
        make

- If this completes without errors, you will now see a file called `vic_image.exe` in this directory. `vic_image.exe` is the executable file for the model.

## Run VIC

At the command prompt, type:

        ./vic_image.exe -g global_parameter_filename.txt

where `global_parameter_filename` = name of the global parameter file corresponding to your project.

The VIC image driver can be run using parallel processing with MPI and/or OpenMP.

!!! Note
    Users are encouraged to consult their system administrator for assistance in configuring the VIC image driver for parallel processing applications.

To run VIC image driver using multiple processors using MPI, type the following instead:

        mpiexec -np $n_proc ./vic_image.exe -g global_parameter_filename.txt

where `n_proc` = number of processors to be used. *Note that different MPI implementations may use different names for the MPI executable such as: `mpirun`, `mpiexec_mpt`, or `mpiexec.hydra`*.

To run the VIC image driver using multiple processors with OpenMP (threads), set the environment variable `OMP_NUM_THREADS`:

        export OMP_NUM_THREADS=8
        ./vic_image.exe -g global_parameter_filename.txt

These two parallelization methods may also be combined using a Hybrid OpenMP/MPI approach. However, that configuration is usually machine, compiler, or scheduler dependent.

## Other Command Line Options

VIC has a few other command line options:

- `./vic_image.exe -v`: says which version of VIC this is
- `./vic_image.exe -h`: prints a list of all the VIC command-line options
- `./vic_image.exe -o`: prints a list of all of the current compile-time settings in this executable; to change these settings, you must edit the appropriate header files (e.g. `vic_def.h` or `vic_driver_shared.h`) and recompile using `make full`.
