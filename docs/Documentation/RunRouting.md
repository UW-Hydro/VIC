# Run Routing Model

## Download Source Code

*   The current routing model code is available in the [Routing Model Source Code](../downloads/Code.md) section of the download page.
*   The code is packaged in a compressed tar archive. Extract the contents via:

    `tar -xvzf filenamee`

    where `filename` = name of the file you downloaded, e.g. `route_code_1.1.tar.gz`

## Compile Routing Model

*   Change directory, cd, to the source code directory (extracted from the compressed tar archive above)
*   There should be 2 subdirectories: "samp_inputs" and "src"; cd to "src"
*   Edit [`rout.f`](RoutingInput.md) and make sure that the parameters NROW, NCOL, NYR, and PMAX are all large enough to contain the dimensions of the basin and the length of the simulation.
*   At the command prompt, type:

    `make`

*   If this completes without errors, you will now see a file called `rout` in this directory. `rout` is the executable file for the routing model.

## Run Routing Model

At the command prompt, type:

    `rout input_filenam_`

where

    `input_filename` of the main input file corresponding to your project

If you'd like to capture the routing model's screen output in a log file, you can do (in C-shell):

    `rout input_filename >& log.txt`

where

`put_filename` = name of the main input file corresponding to your project.

`log.txt` = name of a log file to contain the routing model's screen output.

* * *

Note: This version of the code uses a default value of `UH_DAY=96`, defined in `rout.f`. This was set to ensure the routing code would work in large basins, where routed flows take longer to exit the system (the previous default was 48). All .uh_s files created with one value of UH_DAY will not work if a new value is used in the routing code, and so .uh_s files must be re-created if they were originally made with the previous versions of the code.
