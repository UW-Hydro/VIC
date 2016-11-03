# Running the VIC Classic Driver

## Compiling

- Dependencies:
    The Classic Driver's only dependency is a C compiler that supports the C-99 standard.  We routinely test VIC using the following compilers:

    - GNU (`gcc` version 4+)
    - Clang (`clang` version 3+)

    VIC has also been compiled using these compilers:

    - Intel (`icc`)
    - PGI (`pgcc`)

- If you want to use a compiler other than `gcc`, either edit the Makefile or set the `CC` environment variable, e.g.

        export CC=icc

- Change directory, `cd`, to the "Classic Driver" source code directory and type `make`

        cd vic/drivers/classic
        make

*   If this completes without errors, you will now see a file called `vic_classic.exe` in this directory. `vic_classic.exe` is the executable file for the model.

## Run VIC

At the command prompt, type:

        ./vic_classic.exe -g global_parameter_filename

where `global_parameter_filename` = name of the global parameter file corresponding to your project.

where `global_parameter_filename`  name of the global parameter file corresponding to your project.

## Other Command Line Options

VIC has a few other command line options:

*   `vic_classic.exe -v`: says which version of VIC this is
*   `vic_classic.exe -h`: prints a list of all the VIC command-line options
*   `vic_classic.exe -o`: prints a list of all of the current compile-time settings in this executable; to change these settings, you must edit `vic_def.h` and recompile using `make full`.
