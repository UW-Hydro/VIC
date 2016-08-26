# VIC Model Outputs - Image Driver

## Model Results

The contents of the results files can be controlled by the user, by options/instructions in the global parameter file. If no output file instructions are given in the [global parameter file](GlobalParam.md), [VIC will create default output files with default variables.](DefaultOutputs.md)

[How to Control the Contents of VIC Output Files](OutputFormatting.md)

[List of possible output variables](../../OutputVarList.md) (given in the file `vic_driver_shared_all.h`)

[Default Output Files](DefaultOutputs.md)

[How to Add New Output Variables (that aren't currently available)](../../HowToAddNewOutputVars.md)

## The VIC Runtime Logs

If the `LOG_DIR` variable is provided in the global parameter file, VIC will output its logging to a log file (file name is determined at runtime). The default logging location is `stderr`. The verbosity of these logs can be controlled by setting the `LOG_LVL` variable in the `Makefile`.

## State File (optional)

VIC can save the hydrologic state from any point in the simulation (usually the final state) to a file for the purpose of re-starting the simulation later (as an initial state file). This is useful for simulations that require lengthy spin-up periods or ensemble methods.

[State File Structure](StateFile.md)

The point at which the hydrologic state is saved, and the name/location of the state file, can be specified by the user in the [global parameter file](GlobalParam.md).
