# Routing Model Output

The routing model stores its output in several files, as well as printing informative messages to the screen.

## Output Files

For each station indicated in the [station location file](StationLocation.md), the routing model produces several ASCII output files, containing time series of discharge at the station:

*   `station.day`, `station.day_mm`: daily discharge at station `station`, in units of cfs or mm over the basin, respectively. Format:

    `YYYY MM DD discharge`

*   `station.month`, `station.month_mm`: monthly discharge at station `station`, in units of cfs or mm over the basin, respectively. Format:

    `YYYY MM discharge`

*   `station.year`, `station.year_mm`: annual discharge at station `station`, in units of cfs or mm over the basin, respectively. Format:

    `YYYY discharge`

## Screen Output

The routing model prints the following information to the screen:

*   Number of days and months in the simulation
*   For each station indicated in the [station location file](StationLocation.md):
    *   Number of grid cells upstream of the station
    *   Openings of the UH files
    *   For each upstream cell, name of runoff time series file being read
    *   Any warnings or errors encountered (e.g. missing files, file format errors, etc)

If you'd like to capture the routing model's screen output in a log file, you can do (in C-shell):

`rout input_filename >& log.txt`

where

`input_filename` = name of the main input file corresponding to your project.

`log.txt` = name of a log file to contain the routing model's screen output.
