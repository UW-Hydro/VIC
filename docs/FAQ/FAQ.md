# Frequently Asked Questions

Here are some commonly-asked questions and their answers. Another source of information is the [VIC Users Listserve](http://mailman.u.washington.edu/mailman/listinfo/vic_users). If you join this listserve, you can post questions to the vic user community.

1.  **Does VIC and/or the routing model run on Windows or MacOS?**

    The VIC model and the routing model have been developed for use on linux and unix platforms and as such it will run perfectly fine under OS X. On OS X, you may need to install some additional libraries such as NetCDF. For that you can use any of the OS X package managers such as [MacPorts](https://www.macports.org/index.php), [Homebrew](http://brew.sh), or [Fink](http://www.finkproject.org).

    To use VIC and/or the routing model on a Windows platform, we suggest downloading a free UNIX emulator such as [Cygwin](http://www.cygwin.com) and compiling/running these models within this emulator. Please do not ask us for help with installing or using Cygwin.

    We do not have the resources to test or support running VIC or routing model on any non-\*nix platforms. If you encounter problems while running VIC and/or the routing model on Windows or any other non-\*nix environment, we will not be able to offer you support. If you do think you have found a bug in VIC while running on a non-\*nix platform, please replicate the error on a \*nix system before reporting it.

2.  **Is VIC available in Fortran?**

    No. VIC is written and maintained in C.

3.  **Is the routing model available in C?**

    For the most part, no. We do know of a version of the model written in C that is being tested by a group here at UW, but there are some discrepancies between it and the official Fortran version. When these are resolved, we may post the C version. If you would like to try the C version, take a look [here](https://github.com/UW-Hydro/VIC_Routing).

4.  **What are the hardware requirements for VIC and the routing model? Do they take advantage of parallel processing / can they run in a multi-threaded environment?**

    VIC's code is not explicitly written for parallel processing. Nevertheless, VIC's design will allow you to implement a "poor man's parallelization". Because VIC simulates each grid cell independently of the others, you can break up your study domain into separate sub-domains (these could be sub-basins or simply any arbitrary grouping of grid cells) and simulate each of these sub-domains separately, in parallel (on different nodes of a cluster, for example). Dividing the domain into sub-domains and combining them again later should be fairly easy, since VIC's input forcings and output results are organized into one file per grid cell. All you need to do is create separate soil parameter files, one for each sub-domain (these can even be copies of the original whole-domain soil parameter file, with the appropriate grid cells given a "1" as the first value and all other grid cells given a "0"). You will also need to create separate global parameter files. These would all be identical save for the soil parameter file that they point to. Of course the routing model must access all grid cells in a single run, but it runs very quickly in comparison to the VIC run.

    VIC is very I/O intensive - by default it appends to the grid cell's output file at the end of every time step, although the more recent versions (4.0.6, 4.1.1) allow you to aggregate sub-daily results to a daily time step before writing them, thus cutting down the number of file writes. Still, if you are running on a cluster that has many nodes but a common disk array with a small number of disk controllers, a performance bottleneck will occur in the writing of output files. In this case (which we have here in our lab), the best thing to do is write the output locally (assuming the nodes have their own storage) during the simulation, and then at the end of each node's job, have the job copy all the output back to the common array in a single command, which will be much more efficient in terms of block sizes than VIC's single-record writing.

    On the other hand, VIC is not very memory-intensive, since it only stores the state of a single grid cell in memory at a given time. I don't know off hand how much memory is uses, though that would be a good statistic to measure.

5.  **Are there any Matlab tools written for analyzing and/or plotting VIC results?**

    No. We generally use python to do our analysis. Check the [VIC Users Listserve](http://mailman.u.washington.edu/mailman/listinfo/vic_users).

6.  **How can I save VIC's screen output in a file?**

    The Classic and Image drivers both include the option to write the runtime logs to file. To use this option, set the `LOG_DIR` variable in your *global parameter file*.

    This strategy also works for the routing model screen output.

7.  **VIC runs fine in water-balance mode, but when I set FULL_ENERGY to TRUE, it crashes – why?**

    If you are running an older version of VIC (pre-4.1.1), one explanation might be that your soil parameter file has “nodata” values (e.g. -9999) for bubbling pressure, which is not needed for water balance runs. Bubbling pressure is necessary for energy-balance runs, but older versions of VIC did not check whether the values of bubbling pressure in the soil parameter file were valid. 4.1.1 does check and exits with an error message if bubbling pressure is not valid. If you need a quick-and-dirty way to generate values of bubbling pressure so that you can start an energy-balance run, you can use an [approximate relationship between _bubble_ and values of _expt_](../Documentation/Definitions.md) in your soil parameter file, derived by applying a linear regression to the values from Table 5.3.2 in Handbook of Hydrology and converting lambda to Brooks-Corey's _n_ via Table 5.1.1 in Handbook of Hydrology.

8.  **VIC reports that the energy balance failed to converge in _snow_intercept_, _calc_atmos_energy_bal_, _calc_surf_energy_bal_, or _solve_T_profile_ - what can I do?**

    Failures to converge can arise for a number of reasons. Causes can include bad input parameters or forcings, poor choice of model timestep length, and instability of the model equations themselves. Here are some steps to take:

    1.  Check your input parameters and meteorological forcings - are they reasonable for this grid cell? Are they in the correct format?
    2.  Check your model time step or snow step - does the problem disappear if you reduce the time step length? (Note: currently VIC cannot support time step lengths < 1 hour)
    3.  If your inputs are reasonable, consider setting TFALLBACK to TRUE in the *global parameter file*. Now, the simulation will continue using the previous time step's temperature value. The number of instances in which the previous step's T value was used will be reported at the end of each grid cell's simulation. In addition, several output variables are available to monitor when these instances occurred: OUT_TFOL_FBFLAG, OUT_TCAN_FBFLAG, OUT_SURFT_FBFLAG, OUT_SOILT_FBFLAG.
    4.  If you use the TFALLBACK option, make sure to output the various temperatures simulated by the model (OUT_TFOLIAGE, OUT_TCANOPY, OUT_SURFT, OUT_SOILT) and check that these values look reasonable. In some cases, convergence failures occur long after the model has ventured into unphysical behavior, and using the previous time step's T value simply preserves these unphysical temperatures.
    5.  Check if these behaviors are a [known issue under investigation](../Development/ModelDevelopment.md). If you do not see the issue listed, please create an issue on Github[contact us](https://github.com/UW-Hydro/VIC/issues) and tell us about the issue.

9.  **Why does actual evapotranspiration (ET) sometimes exceed potential evapotranspiration (PET)?**

    VIC's actual evapotranspiration is the sum of canopy evaporation, transpiration, soil evaporation, canopy snow sublimation, and ground snow sublimation, averaged over the grid cell.  VIC's OUT_PET (release 5.0 onwards; OUT_PET_NATVEG for releases 4.1-4.2) computes PET as the ET that the current landscape (with its current vegetation, architectural resistance, and LAI) would produce in the absence of limitations from soil moisture, vapor pressure deficit, temperature, or insolation.  But VIC's canopy evaporation and sublimation are computed using a canopy resistance of 0, which allows for much higher ET rates than when architectural resistance and LAI are taken into account.  Thus, when there is water stored in the canopy, actual ET can exceed OUT_PET.  This generally happens during a small fraction of the time of a simulation (immediately following rain events, for example).

10.  **Relative humidity from VIC 4.0.x sometimes exceeds 1.0 (100%) - what can I do?**

    Upgrade to VIC 4.1.1 or later. This problem does not exist in later versions of the model.

11.  **VIC 4.0.x sometimes gives large water balance errors when COMPUTE_TREELINE is TRUE - what can I do?**

    Please ignore these water balance errors. They are accounting errors only, resulting from a bug in the water balance check, rather than any problems in any water balance terms. These errors only occur in the water balance of the first time step of the simulation. They only occur in grid cells that have elevation bands above the treeline. These errors do not occur if COMPUTE_TREELINE is FALSE. These errors do not occur in VIC 4.1.1 or later.

12.  **Why does VIC give different results when compiled with optimization flags or compiled 64-bit instead of 32-bit?**

    We have heard various reports of VIC's results changing when users re-compile the code with different optimization flags or move from 32-bit to 64-bit platforms. At present we do not know the reasons for this. What we do know is:

    *   These differences in results may result from memory issues. This is especially likely if VIC gives different results from one run to the next using the exact same inputs. But we are not 100% certain that memory issues are the reason.
    *   These differences in results may be less pronounced when running VIC 4.1.1 or later. VIC 4.1.1 and later have been subjected to more comprehensive memory testing than earlier versions, and many bugs having to do with uninitialized variables, buffer overruns, incorrect freeing of allocated memory, etc. have been found and fixed. However even VIC 4.1.1 and later yield different results under different compiler options and on different platforms.
    *   Some of the differences may result from using inappropriate compiler option combinations, such as "-g" (debugging symbols) and "-O3" simultaneously. You should never combine the "-g" option with any compiler optimization options. Admittedly, this can make it difficult to debug an optimized version of the code.
    *   Some of the differences in results may in fact be introduced by the compiler when it attempts to optimize the code.
    *   Simply turning off the "-g" option without adding an optimization option will make the code run faster - this may be sufficient for your purposes and may reduce any differences in results.
    *   We have not yet performed rigorous testing of any VIC versions on 64-bit platforms. At this point, users must assume some risk when running VIC in a 64-bit environment. We hope to look into this further in the future.
    *   VIC's numerical schemes may be partially responsible for differences. VIC uses simple first-order forward difference schemes for most of its computations, although VIC 4.1.1 and later have the option to use an implicit scheme to model soil temperatures.

    In summary, probably the best thing you can do at this point is to upgrade to the latest version of VIC and see if turning off the "-g" option is sufficient to speed up your code without invoking any compiler optimization options. If running with FROZEN_SOIL or FULL_ENERGY set to TRUE, setting IMPLICIT to TRUE may make VIC's results more stable as well.

    It is worth noting that users have reported similar problems with other hydrologic models. See Martyn Clark's "numerical daemons" paper.

13.  **I want to output a variable, but I don't see it listed in the output variable list. What should I do?**

    VIC allows users to select which output variables to write to their output files. Variables that can be selected are listed [here](../Documentation/OutputVarList.md) as well as in the VIC source code file `vic_def.h` (listed as `OUT_*`). However, not all variables that VIC simulates are available for output. And of course if you are changing VIC's physics, any new variables you add are also not available for output (yet). In either of these cases, please see the [instructions on how to add a new output variable to VIC](../Documentation/HowToAddNewOutputVars.md).
