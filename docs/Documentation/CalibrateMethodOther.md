# Other VIC Model Calibration Methods

Calibration of the VIC Hydrologic Model is completed using an optimization process.

*   [Calibration Algorithms](#Algorithm)
*   [Discussion of Calibrater Applications](#Discussion)
*   [Calibration Procedure](#Procedure)

## Calibration Algorithms

Recently, attempts have been made to use numerical optimization techniques to automate the calibration procedure for the VIC model. Many papers have been written on the subject of optimization of rainfall-runoff models, but fewer attempts have been made to optimize the parameters of more complex hydrologic models. Currently two methods have been tested on the VIC macroscale model: the random autostart simplex method and a genetic optimization scheme. Both optimization methods attempt to minimize the differences between simulated and observed discharge records. The magnitude of the differences was given a value by computing the Nash-Sutcliffe R<sup>2</sup> coefficient (where 1 represents a perfect match and smaller values represent worse matches). The optimization schemes attempted to minimize the negative of this coefficient (-1 would therefore be a perfect match).

The random autostart simplex method starts by selecting a large number of random parameter sets, solving the model at each one. It then selects the best set of parameters from the randomly generated sets and uses those to start the simplex minimization algorithm. A good description of the simplex method can be found in the Numerical Recipes books. A simplified description is that the algorithm tries to corral the minimum within a geometric shape, the simplex, with N+1 apexes (N being the number of parameters being optimized). Once a minimum has been contained the algorithm begins to minimize the volume of the simplex, until all of its apexes are within a specified tolerance of each other. The problem with the straight simplex method is that it is very likely to find a local minimum, and not the global minimum that truly optimizes the model. To combat this problem the random autostart process is run several times, each time producing a new set of initial parameters for the model. If run enough times the algorithm should eventually locate the global optimization parameters.

Genetic optimization was developed as a way to get around the limitations of the simplex method, by following up on all of the good solutions, and hopefully not getting trapped in local minimums. The laws of natural selection are applied to the parameter sets driving the survival of the fittest. First the original population is created by randomly generating a large set of parameters. These are then sorted and assigned a probability of reproducing; individuals with the lowest -R<sup>2</sup> values are given the highest probabilities. Then a new generation is created by randomly selecting parents from the previous generation (based on the probability of their reproducing). The parameter sets for both parents are encoded as bit strings (DNA), fragments are selected from both parents, and recombined to form the child's parameters. The model is solved for all of the child parameter sets, and then they are allowed to breed a new generation. Over the course of several generations the weaker (less optimized) parameters are bred out of the population, with enough generations the algorithm should focus in on the optimized parameters.

## Discussion of Calibrater Applications

All simulations used for these tests optimized the model using five parameters: b<sub>i</sub>, D<sub>s</sub>, W<sub>s</sub>, D<sub>2</sub>, and D<sub>3</sub>. Optimization was conducted on the calibration period of 1980-1984, though the first year was not used in computing the R<sup>2</sup> values since it was the ramp up period for the model. Most of the simulations were conducted using the full energy balance mode with a three hour time step. Though tests were conducted on various catchments within the Upper Mississippi Basin, the frozen soil algorithm was not used since it significantly increases the computational time, and makes multiple optimization passes prohibitively expensive.

The simplex method was applied using random autostart populations of 75-100 paremeter sets. The entire cycle was repeated from 5-10 times for each subbasin. For most autostarts the simplex method converged to an R<sup>2</sup> tolerance of 0.001 within 75 iterations. Each autostart yielded different R<sup>2</sup> values (usually within +/- 0.1), and different parameter sets. In some cases a single parameter could change significantly between optimized sets, and yet yield fairly similar R<sup>2</sup> values. Plotting simulated discharges emphasizes the differences between the parameter sets, and shows that some sets are less physically sound (_e.g._ flat baseflow between peaks) than others, so it is left for the user to determine the set of optimized parameters that best maintains the physical reponse of the basin.

Genetic optimization was applied to the Wisconsin River (14 1/2 degree grid cells), for the same time period and with the same model options used for the simplex tests. The population was held to 100 parameter sets, and the model was run for 20 generations. Typically the genetic optimizer must be run for significantly more generations to find the true global minimum, but after 20 generations the best solution it had found was still in the original population. After finding a minimum R<sup>2</sup> of 0.70 in the first population, the best solution in the following generations ranged from 0.64 to 0.69\. This is hardly an improvement and simulations had already consumed 6 days. Therefore, the initial assessment is that the VIC model does not lend itself to optimization by the genetic algorithm and that the simplex method should be employed until a better optimization algorithm is found.

## Calibration Procedure

The following is a discussion of how the optimization was implemented for the VIC macroscale hydrologic model. Included are the scripts and source code used to run the optimizer, plus suggested modifications to both the VIC model and routing model source codes to make them function more efficently. The scripts and source code mentioned may be found in the file [calibrate_other.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/calibrate_other.tgz), available in the Calibration section of the [download page](../SourceCode/Download.shtml).

*   Determine if optimization is feasible in your situation:
    *   To determine the time necessary to run the optimization start by running the model with /usr/bin/time. Estimate the optimization time using the formula:  
         `(user time) * (Ntry) * (Ninit) * (80)`  

        *   **Ntry** is the number of optimizations to run (5-25).
        *   **Ninit** is the number of initial random parameter sets used to find the best starting points for the optimizer (75-100).
        *   **80** is the approximate number of iterations needed for each optimization to find a solution (this is likely to increase if you decrease Ninit).

*   Random autostart simplex method program
    *   The program optimize_vic.c uses the random autostart simplex method to optimize the function given by the [model run script](#ModelRunScript). It asks for two such scripts, the first the _model start script_ is used for the autostart sequences, while the _model optimize script_ is actually optimized in the simplex. This allows the user to generate the startup parameter sets using a course time step, or fewer model options before optimizing parameters for the model setup that will actually be used.
    *   Output from all runs (parameters as well as -R<sup>2</sup> values) will be stored in the optimization file. When the optimizer has finished each optimization, it will also write the set of optimized parameters to the file.
    *   Before compiling the program, check the values of the following parameters:
        *   **Ntry** is the number of times the optimizer will be run. Each time it is run it will produce another set of optimized parameters. The larger Ntrys is the more likely the optimizer will find the global minimum, but it will also take longer to run. A value of 5-10 should produce enough parameter sets to be able to select a decent calibration.
        *   **Ninit** is the number of random parameter sets to be generated before starting the optimizer. It can be set as low as N+1 (where N is the number of parameters you are calibrating), but higher values will produce more parameter sets to choose from when starting the optimizer. A value of 75-100 will produce enough sets that the optimizer should be able converge on an acceptable solution within 75-80 iterations.
    *   Since this program is liable to run for many days depending on the size of the basin, and the number of iterations it makes, output from the program should be directed to /dev/null, or an error log file (for later viewing), and the program backgrounded. _example run statement:_  
         `optimize_vic run_wis_3hr_EWB.script run_wis_1hr_fs_ewb.script wis_opt_log_021899 >& wis_error_log_021899 &`
*   Model run script
    *   The script run_opt_wis_calib.script changes the parameter files, runs the model, runs the routing model and then computes the R<sup>2</sup> value. The script is given the current values of the optimization parameters, and writes the R<sup>2</sup> value to a file which is read by the [optimization model](#OptiModel).
    *   The following items need to be set in the script to run it for a different setup:
        *   [the parameter modification script](#ParamAdjust)
        *   [the model global parameter file](#ControlFile)
        *   the directory where the optimization results are written
        *   [the routine to compute the R<sup>2</sup> value](#ComputeR2)
*   Parameter ajustment routine
    *   The script calibrate_wis_opti.script takes the values for five VIC parameters (currently set to b<sub>i</sub>, D<sub>s</sub>, W<sub>s</sub>, D<sub>2</sub>, and D<sub>3</sub>) and outs the new parameter files in the optimization directory.
    *   NOTE: this script is designed to work with the ARC/INFO ASCII grid soil parameter files. But it should be relatively simple to create scripts or programs to do the same thing for the ASCII table files.
    *   The file wis_file_list_opti is the file list used by the example model set-up.
*   Computing the Nash-Sutcliffe R<sup>2</sup> coefficient
    *   The program compute_R2.c computes the Nash-Sutcliffe R<sup>2</sup> coefficient between the daily observed flow, and the simulated daily flow.  
         `compute_R2 RESULTS/CURRENT_RUN/DAILY/WISCR.day ROUTING/WISCR/wisc_daily.txt 366 1460 0.80 $OUTFILE > $R2FILE`
    *   $OUTFILE contains only the R<sup>2</sup> value and is read by the [model run script](#ModelRunScript) and returned to the [optimization program](#OptiModel).
*   Model control file
    *   wis_global_ewb_calib_opt is an example control file which runs the energy and water balane VIC model at a three hour time step. See the technical notes on [model time step](TechnicalNotes/Timestep.shtml) and [frozen soil time step](TechnicalNotes/TimeStopfrozsoil.shtml) issues for more about selecting a suitable time step.
*   Streamlining the VIC model code for optimization
    *   To reduce model run time, it helps to minimize VIC file writes. This can be done by:
        1.  Specifying daily output, i.e. OUTPUT_STEP 24 in the [global parameter file](GlobalParam.shtml).
        2.  Outputting the bare minimum for input into the routing model, i.e. specifying one output file per grid cell, with only the following variables: OUT_PRECIP, OUT_EVAP, OUT_RUNOFF, and OUT_BASEFLOW. This can be specified in the [global parameter file](GlobalParam.shtml) via [flexible output configuration](OutputFormatting.shtml).
