# Routing: Flow Direction File

The flow direction file tells the routing model how all of the grid cells are connected in the routing net. It is also the only input file with a header. The header tells the routing model the lower left latitude and longitude, the number of rows and columns, and the grid cell resolution. An example of a flow direction file is shown below:

```
ncols           22
nrows           20
xllcorner       -97.000
yllcorner       38.000
cellsize        0.50
NODATA_value    0
0 0 0 5 5 5 6 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 4.0 4.0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 5 3 5 5 6 7 5 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 5 4.0 5 5 5 5 6 7 5 5 5 7 5 0 0 0 0 0 0
0 0 5 4.0 5 4 5 5 7 7 5 6 7 5 6 0 0 0 0 0 0
4.0 5 7 4.0 3 5 5 7 5 5 7 4 5 6 0 0 0 0 0 0
4.0 3 4.0.4.0.4 3 5 6 7 7 4 5 6 0 0 0 0 0 0
0 4.0 4.0 3 1 2 4.0 4 5 7 5 5 0 0 0 0 0 0 0
0 0 3 4 3 1 8 4 5 4.0 5 5 3 5 6 5 0 0 0 0 0
0 0 0 3 4 5 5 4 5 4 3 5 6 7 7 5 5 7 0 0 0 0
0 0 0 4 3 5 5 3 4.0.4 3 4 3 5 5 7 7 0 0 0 0
0 0 0 4.0 5 4.0.4 3 4.0.4 5 3 5 7 7 0 0 0 0
0 0 0 0 4.0 4.0 4.0 4 5 3 6 7 7 5 5 6 0 5 7
0 0 0 0 0 4.0 3 4 3 5 3 6 7 6 5 7 7 7 7 7 7
0 0 0 0 0 0 0 0 4.0 4 5 7 5 6 7 1 7 1 7 0 0
0 0 0 0 0 0 0 0 3 4.0 5 5 6 7 7 7 6 0 0 0 0
0 0 0 0 0 0 0 0 0 4.0 4 5 7 1 7 7 0 0 0 0 0
0 0 0 0 0 0 0 0 0 3 2 4 5 7 7 1 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 4.0 5 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

```

Where:

*   ncols is the number of columns in the file,
*   nrows is the number of rows in the file,
*   xllcorner is the longitude of the lower left corner of the grid,
*   yllcorner is the latitude of the lower left corner,
*   cellsize is the resolution of the grid, and
*   NODATA_value is the value that represents missing or unused grid cells.

Flow from each grid cell is given by a by a number:

1.  = north
2.  = northeast
3.  = east
4.  = southeast
5.  = south
6.  = southwest
7.  = west
8.  = northwest
