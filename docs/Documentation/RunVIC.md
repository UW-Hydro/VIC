# Running the VIC Model

## Source Code

*   To download the latest VIC source code, proceed to the [download page](../SourceCode/Code.md).
*   The code is packaged in a compressed tar archive. Extract the contents via:

`tar -xvzf filename`

    where `filename` = name of the file you downloaded, e.g. VIC_code_4.2.b.tar.gz

## Compile VIC

*   Change directory, cd, to the source code directory (extracted from the compressed tar archive above)
*   At the command prompt type:

`make`

*   If this completes without errors, you will now see a file called `vicNl` in this directory. `vicNl` is the executable file for the model.

## Run VIC

At the command prompt, type:

`vicNl -g global_parameter_filename`

where `global_parameter_filename` = name of the global parameter file corresponding to your project.

If you'd like to capture VIC's screen output in a log file, you can do (in C-shell):

`vicNl -g global_parameter_filename >& log.txt`

where `global_parameter_filename`  name of the global parameter file corresponding to your project.

`log.txt` = name of a log file to contain VIC's screen output.

## Other Command Line Options

VIC has a few other command line options:

*   `vicNl -v`: says which version of VIC this is
*   `vicNl -h`: prints a list of all the VIC command-line options
*   `vicNl -o`: prints a list of all of the current compile-time settings in this executable; to change these settings, you must edit `vicNl_def.h` and recompile using `make clean; make`.
