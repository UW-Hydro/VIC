# VIC Model Outputs

## Model Results

The contents of the results files can be controlled by the user, by options/instructions in the global parameter file. If no output file instructions are given in the [global parameter file](GlobalParam.md), [VIC will create the same 2 or 3 output files as in earlier versions, by default.](DefaultOutputs.md)

[How to Control the Contents of VIC Output Files](OutputFormatting.md)

[List of possible output variables](OutputVarList.md) (given in the file vicNl_def.h)

[Default Output Files](DefaultOutputs.md)

[Other Output File formats](OtherOutputFiles.md)

[How to Add New Output Variables (that aren't currently available)](HowToAddNewOutputVars.md)

## Output to the Screen

VIC sends messages to the screen regarding status of the simulation, warnings, and errors. The amount of verbage sent to the screen can be controlled by the appropriate setting of VERBOSE in the file vicNl_def.h. We recommend redirecting VIC's screen output to a log file for later reference; this can be useful in diagnosing model behavior.

If you'd like to capture VIC's screen output in a log file, you can do (in C-shell):

  vicNl -g global_parameter_filename >& log.txt

where

`global_parameter_filename` = name of the global parameter file corresponding to your project.

`log.txt` = name of a log file to contain VIC's screen output.

When VERBOSE is TRUE, VIC prints the following output to the screen:

*   List of parameter files that are opened for reading
*   List of the values of compile-time and run-time options for this simulation (this is important for diagnosing the causes of model behavior)
*   For each grid cell, lists:
    *   Cell ID, latitude, and longitude
    *   Forcing file opened for reading
    *   Output files opened for writing <limaximum and="" cumulative="" average="" moisture="" energy="" balance="" residual="" errors="" as="" they="" are="" encountered="" through="" time="" during="" the="" simulation="" <li="">Warnings and error messages encountered during the simulation</limaximum>
    *   Final total moisture and energy balance residual errors

## State File (optional)

VIC can save the hydrologic state from any point in the simulation (usually the final state) to a file for the purpose of re-starting the simulation later (as an initial state file). This is useful for simulations that require lengthy spin-up periods or ensemble methods.

[State File Structure](StateFile.md)

The point at which the hydrologic state is saved, and the name/location of the state file, can be specified by the user in the [global parameter file](GlobalParam.md).
