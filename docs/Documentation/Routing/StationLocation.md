# Routing: Station Location File

The station location file tells the routing model from which grid cells to produce output flow data. Any number of stations may be defined within the basin, as well as a single basin outlet, where the routing network leaves the defined basin. Each line defining a station is followed by another that tells the routing model whether or not a uh_s file has been generated for the current station location. If set to NONE the routing model generates a new uh_s file in the current directory, otherwise it will read the defined uh_s file. An example station location file is shown below:

```
1 UPMIS               14  2  -9999
NONE
1 WISCR               14 11  -9998
WISCR.uh_s
1 ANOKA                8 15  -9998
store/ANOKA.uh_s
```

The file contains information in the following columns:

1.  The first column indicates whether the station is active (1 = active, 0 = inactive),
2.  the second column is the basin name which will be used for all of the output files (the first 5 characters are used to form the root of the output files),
3.  the third column is the column, from the left, in the gridded direction file where the station is located,
4.  the fourth column is the row, from the bottom, in the gridded direction file where the station is located, and
5.  the fifth column is the basin area (required, but not used at present).
