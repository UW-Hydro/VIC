VIC Python Driver
================

The VIC Python driver exposes the `vic_run` and `vic_driver_shared` function APIs through a Python interface. It may be used to test individual algorithms or to wrap the `vic_run` function.

!!! Warning
    The VIC Python Driver is experimental and is still under development.

## Building the Python Driver

The Python driver uses the Python `setuptools` package to handle the compilation and install of VIC. It has the following dependencies:

- C compiler with C-99 (e.g. gcc 4+ or clang 3+).
- Python version 2.7+ or 3.4+, with the following packages installed:
    - [CFFI](http://cffi.readthedocs.org/en/latest/index.html) version 1.2 or greater

- To set the C compiler used in the Python build process,

        export CC=icc

- Change directory, `cd`, to the "Python Driver" source code directory and run the `setup.py` script:

        cd vic/drivers/classic
        python setup.py install

!!! Note
    Note that the Python driver is built as a library and is intended to be linked to the Python coupler by the Python build system, it is not a stand-alone executable.

- Use your python interpreter to interface with the installed VIC library:

```
$ ipython
Python 3.5.1 |Continuum Analytics, Inc.| (default, Dec  7 2015, 11:24:55)
Type "copyright", "credits" or "license" for more information.

IPython 4.2.0 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: from vic import lib as vic_lib

In [2]: vic_lib.print_license()

  Variable Infiltration Capacity (VIC) macroscale hydrologic
  model version 5.0.0, Copyright (C) 2016 Computational
  Hydrology Group, Dept. of Civil and Environmental Engineering,
  University of Washington.  VIC comes with ABSOLUTELY NO
  WARRANTY. This is free software, you may redistribute it
  under certain conditions; see LICENSE.txt for details.

  Report Bugs and Issues to : https://github.com/UW-Hydro/VIC/issues
  VIC Users Email Listserve : vic_users@u.washington.edu
```
