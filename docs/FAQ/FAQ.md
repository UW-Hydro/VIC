# Frequently Asked Questions

Here are some commonly-asked questions and their answers. Another source of information is the [VIC Users Listserve](http://mailman.u.washington.edu/mailman/listinfo/vic_users). If you join this listserve, you can post questions to the vic user community.

1.  **Does VIC and/or the routing model run on Windows or MacOS?**

    The VIC model and the routing model have been developed for use on linux and unix platforms and as such it will run perfectly fine under OS X. On OS X, you may need to install some additional libraries such as NetCDF. For that you can use any of the OS X package managers such as [MacPorts](https://www.macports.org/index.php), [Homebrew](http://brew.sh), or [Fink](http://www.finkproject.org).

    To use VIC and/or the routing model on a Windows platform, we suggest downloading a free UNIX emulator such as [Cygwin](http://www.cygwin.com) or [Docker](../Development/Docker.md) and compiling/running these models within this emulator. Please do not ask us for help with installing or using Cygwin.

    !!! note
        We do not have the resources to test or support running VIC or routing model on any non-\*nix platforms. If you encounter problems while running VIC and/or the routing model on Windows or any other non-\*nix environment, we will not be able to offer you support. If you do think you have found a bug in VIC while running on a non-\*nix platform, please replicate the error on a \*nix system before reporting it.

2.  **Is VIC available in Fortran?**

    No. VIC is written and maintained in C.  However, the [CESM driver](../Documentation/Drivers/CESM/CESMDriver.md) includes Fortran bindings using the [`ISO_C_BINDING` standard](https://gcc.gnu.org/onlinedocs/gfortran/Interoperability-with-C.html). If you are looking to use VIC in conjunction with Fortran, the CESM driver source code would be a good place to start.

3.  **What are the hardware requirements for VIC? Do they take advantage of parallel processing / can they run in a multi-threaded environment?**

    As described above, VIC is designed to run on any linux or unix system. VIC is capable of being run on a laptop or on large distributed super-computers.

    ***Parallelization:***

    - The *[Classic Driver](../Documentation/Drivers/Classic/ClassicDriver.md)* code is not explicitly written for parallel processing. Nevertheless, its design will allow you to implement a "poor man's parallelization". Because VIC simulates each grid cell independently of the others, you can break up your study domain into separate sub-domains (these could be sub-basins or simply any arbitrary grouping of grid cells) and simulate each of these sub-domains separately, in parallel (on different nodes of a cluster, for example). Dividing the domain into sub-domains and combining them again later should be fairly easy, since VIC's input forcings and output results are organized into one file per grid cell. All you need to do is create separate soil parameter files, one for each sub-domain (these can even be copies of the original whole-domain soil parameter file, with the appropriate grid cells given a "1" as the first value and all other grid cells given a "0"). You will also need to create separate global parameter files. These would all be identical save for the soil parameter file that they point to. Of course the routing model must access all grid cells in a single run, but it runs very quickly in comparison to the VIC run.
    - The *[Image](../Documentation/Drivers/Image/ImageDriver.md) and [CESM](../Documentation/Drivers/CESM/CESMDriver.md) Drivers* utilize a the MPI standard for parallelization. This allows VIC to be efficiently run on larger clusters and super-computers. For information on how to run the VIC Image Driver in parallel, see the ["how to run VIC page"](../Documentation/Drivers/Image/RunVIC.md).

    ***Memory Requirements***

    - The *[Classic Driver](../Documentation/Drivers/Classic/ClassicDriver.md)* is not very memory-intensive, since it only stores the state of a single grid cell in memory at a given time.
    - The memory usage of the *[Image](../Documentation/Drivers/Image/ImageDriver.md) and [CESM](../Documentation/Drivers/CESM/CESMDriver.md) Drivers* scales with the number of grid cells on each processor.

4.  **Are there any Matlab tools written for analyzing and/or plotting VIC results?**

    Not that we know of. We generally use python to do our analysis. Check the [VIC Users Listserve](http://mailman.u.washington.edu/mailman/listinfo/vic_users).

5.  **How can I save VIC's screen output in a file?**

    The Classic and Image drivers both include the option to write the runtime logs to file. To use this option, set the `LOG_DIR` variable in your *global parameter file*.

6.  **VIC reports that the energy balance failed to converge in _snow_intercept_, _calc_atmos_energy_bal_, _calc_surf_energy_bal_, or _solve_T_profile_ - what can I do?**

    Failures to converge can arise for a number of reasons. Causes can include bad input parameters or forcings, poor choice of model timestep length, and instability of the model equations themselves. Here are some steps to take:

    1.  Check your input parameters and meteorological forcings - are they reasonable for this grid cell? Are they in the correct format?
    2.  Check your model time step or snow step - does the problem disappear if you reduce the time step length?
    3.  If your inputs are reasonable, consider setting `TFALLBACK` to TRUE in the *global parameter file*. Now, the simulation will continue using the previous time step's temperature value. The number of instances in which the previous step's T value was used will be reported at the end of each grid cell's simulation. In addition, several output variables are available to monitor when these instances occurred: `OUT_TFOL_FBFLAG`, `OUT_TCAN_FBFLAG`, `OUT_SURFT_FBFLAG`, `OUT_SOILT_FBFLAG`.
    4.  If you use the `TFALLBACK` option, make sure to output the various temperatures simulated by the model (`OUT_TFOLIAGE`, `OUT_TCANOPY`, `OUT_SURFT`, `OUT_SOILT`) and check that these values look reasonable. In some cases, convergence failures occur long after the model has ventured into unphysical behavior, and using the previous time step's T value simply preserves these unphysical temperatures.
    5.  Check if these behaviors are a [known issue under investigation](../Development/ModelDevelopment.md). If you do not see the issue listed, please create an issue [on Github](https://github.com/UW-Hydro/VIC/issues) and tell us about the issue.

8.  **Why does actual evapotranspiration (ET) sometimes exceed potential evapotranspiration (PET)?**

    VIC's actual evapotranspiration is the sum of canopy evaporation, transpiration, soil evaporation, canopy snow sublimation, and ground snow sublimation, averaged over the grid cell.  VIC's `OUT_PET` (release 5.0 onwards; `OUT_PET_NATVEG` for releases 4.1-4.2) computes PET as the ET that the current landscape (with its current vegetation, architectural resistance, and LAI) would produce in the absence of limitations from soil moisture, vapor pressure deficit, temperature, or insolation.  But VIC's canopy evaporation and sublimation are computed using a canopy resistance of 0, which allows for much higher ET rates than when architectural resistance and LAI are taken into account.  Thus, when there is water stored in the canopy, actual ET can exceed `OUT_PET`.  This generally happens during a small fraction of the time of a simulation (immediately following rain events, for example).

9.  **I want to output a variable, but I don't see it listed in the output variable list. What should I do?**

    VIC allows users to select which output variables to write to their output files. Variables that can be selected are listed [here](../Documentation/OutputVarList.md) as well as in the VIC source code file `vic_driver_shared.h` (listed as `OUT_*`). However, not all variables that VIC simulates are available for output. And of course if you are changing VIC's physics, any new variables you add are also not available for output (yet). In either of these cases, please see the [instructions on how to add a new output variable to VIC](../Documentation/HowToAddNewOutputVars.md).
