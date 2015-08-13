# Preparation of Routing Model Parameter Files

This page presents a short description of how to derive all the necessary [input files for the Routing model](../Routing/RoutingInput.md). Some knowledge of basic Arc/Info commands and operations is required, but hints as to which commands to use are provided. Be sure to familiarize yourself with the operation and basic theory behind the routing model (Lohmann et. al. 1998a, b).

# Basin Delineation

Most of the routing model's parameter files are Arc/Info-style ascii grids and are based on some of the files you produced when delineating the basin in Arc/Info. Instructions for delineating a basin from a DEM in Arc/Info are provided [here](../BasinDelineation.md).

# Fraction File

The fraction file in the routing model is used on the boundaries of the basin delineation so that a gridcell can have less area than 1.0, or 100%. There are several ways of doing this and only one way is presented here. This section explains how to produce the fraction file using the fine-resolution DEM produced when you [delineated your basin](../BasinDelineation.md). All programs referenced are available in the file [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz), also available under the "Preparing Routing Model Inputs" section of the [download page](../../SourceCode/Code.md).

Procedure:

1.  First, create a basin mask at the resolution of the DEM file.
    1.  Set all the values in the DEM cover (for example, for the Columbia basin, let's call it "colum_dem") to ones (1)

        `GRID: colum_ones1 = con(colum_dem > 0,1,0)`

    2.  Now you have to make a grid that covers a larger area than your real basin boundaries and set all values in this grid to zeros (0).
        1.  `GRID: gridclip na_dem colum_box COVER box_polygon` (where "box_polygon" is a cover consisting of a box spanning a larger lat-lon range than the basin; note: it is sufficient to pad the basin's lat-lon boundaries with one VIC grid cell on all sides of the basin)
        2.  `GRID: colum_zeros1 = con(colum_box > 0,0,1)`
    3.  You now have two grids, one with only ones (1) in it, and one with only zeros (0) in it. Merge these two grids using the MERGE command.

        `GRID: colum_ones_zeros = merge(colum_ones1,colum_zeros1)`

        It is important to have the grid with ones first in the parentheses.

2.  Now you must make a grid that has the final resolution of your VIC model. The grid must cover a larger area than the real basin boundaries and have a unique value in each cell. This can be accomplished via the program create_unique.c:

    `create_unique colum_mask.asc colum_unique.asc`

    where colum_mask.asc = Arc/Info-style ascii grid file containing the VIC-resolution mask for the basin (in this case, the Columbia basin); and colum_unique.asc = an output file that is similar to colum_mask.asc but with a unique value in each grid cell.

    Another way to do this is to use a spreadsheet. For example for the Columbia River (1/8th degree) a grid with 94rows and 100 columns was made. Starting in the upper left corner the value 1 was put in, and going right and downwards increasing values for each gridcell. Row 1 col 1 have the value 2, row 1 col 2 have the value 2, etc. The area that this grid covers can be the same as your box_polygon. Put an Arc/Info header with rows,cols,xllcorner, etc. in the top of the file.

3.  Import the unique-value file into Arc/Info:

    `ARC: asciigrid colum_unique.asc colum_unique int`

4.  Project this grid to the correct projection and resample it to the same resolution as the DEM. For this example, call the resampled file colum_resampled1 (GTOPO30 resolution = 30 arc second = 0.0083333333333333 degrees)

    `GRID: colum_resampled1 = resample(colum_unique,0.0083333333333333)`

5.  Now you must use the ZONALMEAN function in GRID. This function calculates the mean of the values given in a zone. See HELP ZONALMEAN for more details.

    `GRID: zmean1 = zonalmean(colum_resampled1,colum_one_zeros)`

6.  Resample this file to the modeling resolution. (1/8 degree = 0.125 degrees)

    `GRID: fraction1 = resample(zmean1,0.125)`

7.  The grid fraction1 now has a value from zero to one (0-1) in all the pixels. This value represents the fraction of the area of each gridcell that is in the basin. Convert this file into ascii in Arc/Info using:

    `ARC: gridascii fraction1 fraction1.asc`

8.  This file (fraction1.asc) is your fraction file and can be used directly by the routing model. Take a look at it and check that all the gridcells in the middle of the basin have only ones (1) in it, and that the pixels on the boundaries have a number between zero (0) and one (1). The gridcells that are not part of the basin should all be zeros.

# Flow Direction File

This section describes in detail how to get the final flow network file for the routing model. You will need the fraction1 grid from the [Fraction File](#fraction-file) section as well as the DEM. All programs referenced are available in the file [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz), also available under the "Preparing Routing Model Inputs" section of the [download page](../../SourceCode/Code.md).

## Coarse Resolution Basin Mask

1.  If you do not already have a coarse resolution basin mask (created in the [basin delineation](../BasinDelineation.md) step), you need to create one now.

    Alternatively you can create the mask using the fraction file created [above](#fraction-file). Make a new grid of the fraction1 grid that has the same value in all gridcells that are in the basin. The fraction1 file should have the same resolution as for the VIC model that you are building.

    `GRID: coarse_grid = con(fraction1 > 0,1,0)2.  Now use the GRIDPOLY command in Arc/Info to make a polygon that is an outline of your basin.

    `GRID: coarse_outline = gridpoly(coarse_grid,0.00)    This polygon is the coarse-resolution outline of your basin. On the boundaries it will look like cubes or squares, which is because you are using a grid that has the final model resolution.

## DEM

1.  Now use the basin outline and cut out the basin from the DEM.

    `GRID: gridclip na_dem coarse_dem COVER coarse_outline    The coarse_dem grid now has approximately 15x15 1 km pixels in each 1/8th degree box. This is also true for the gridcells on the boundary of the basin. If you are using another model resolution you will have a different set of pixels in you gridcells. The reason for this is because some of the programs that will be used later need the same number of 1 km pixels in each model resolution gridcell. Now we will use only coarse_dem.

2.  Use the FILL command in GRID to fill all the sinks in the coarse_dem grid. If the DEM has sinks in it the flow accumulation routine will not work.
    1.  GRID: fill coarse_dem coarse_fill
    2.  GRID: flowdir1 = flowdirection(coarse_fill)
    3.  GRID: flowacc1 = flowaccumulation(flowdir1)
    4.  ARC: `gridascii flowacc1 flowacc1.asc`

    The file flowacc1.asc contains the accumulated flow values in each 1-km resolution cell. The outlet of your river has the highest value.

3.  In our example case, we want the flow direction file to be at 1/8 degree latitude-longitude resolution. The script make_rout.scr from [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz) will create a flow direction file at your desired resolution from the high-resolution (1 km) flow accumulation file. Also note, this set-up expects your flow accumulation file values to be integers.
    1.  Edit the two lines of the script to define the input flow accumulation file (set ACC = flowacc1.asc) and the output flow direction (at 1/8th degree or other, ex, set OUT = flowdir.asc).
    2.  This script is hard-coded for creating a 1/8th degree flow direction file from a 30 arc-second flow accumulation file. The “15”s in this script are for the 1/8th degree (120/8, where 120 comes from 30-arc second, i.e., 120 30-arc-seconds/degree). The “0.125” is the 1/8th degree desired coarse resolution. These values will need to be changed for other resolutions.
    3.  Compile the programs
        *   flowgen.c (`gcc –lm –o flowgen flowgen.c`)
        *   convert.c (`gcc –lm –o con convert.c`)

        Make sure that flowgen.h is in the same directory as these programs.

    4.  Run make_rout.scr
4.  The awk script read_rout.awk will convert the routing direction file into an ASCII 2 column xy file, which can be used to generate a plot of the routing network using GMT. The plotting script upmiss.rout.script.clr shows how to use the awk script to generate the plot. It also uses the file mask.xy to mask out the basin area.
5.  The program calc_area.c uses the subroutine new_get_dist.c to compute the area of the routed basin. The program uses the flow direction file to get the latitude and longitude of active cells (those with values not equal to the NODATA_value). The resulting cell area is multiplied by the drainage fraction and summed for the total basin area.

    Note that this file has no header. Values in the file are 1 for 100% contributing area. Values greater than 1 indicate more than 100% contribution, while values less than 1 indicate less than 100% of the grid cell contributes.

# Flow Velocity File

The awk file velo.awk on creates a velocity file from the flow direction file:

`awk -f velo.awk _direction_file_ > _velocity_file_`

Use `nawk` on the Sun systems.

# Flow Diffusion File

The awk file diff.awk creates a diff file from the flow direction file:

`awk -f diff.awk _direction_file_ > _diff_file_`

Use `nawk` on the Sun systems.

# Geometry (Xmask) file

The program create_xmask.c creates an xmask file from the flow direction file, and uses a calculation of the actual horizontal, vertical, and diagonal flow distances for each grid cell to produce an xmask file:

*   compile with: `gcc -o create_xmask create_xmask.c -lm`
*   run with: `create_xmask "direction file" "xmask file"``

# Station Location File

This file lists the points at which hydrographs will be created. These are usually the points for which there are observed streamflow data - but be careful because often the published lat/long are slightly incorrect. For basins with many outflow points special scripts may be used, they are semi-automated. All programs referenced are available in the file [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz), also available under the "Preparing Routing Model Inputs" section of the [download page](../../SourceCode/Code.md).

Recommendation of using these scripts is not given unless you have many stations because it may be easier to do it by hand.

compile: `make`

output: specified in input file

"input": hard-coded

two output files hard-coded as well

can be load into Excel: num_cells_upstream

Useful information when deciding what RN and CN to use in the station file:

`lat = yllcorner + RN*cellres - cellres/2`

`long = xllcorner + CN*cellres - cellres/2`

Another useful program that is not actually necessary is gen_submask.c.

`gen_submask gen_input.basin`

This programs creates a mask file for all cells draining to a certain location. It also gives an estimate of the area draining into that point. This program does need the vic output, be fore-warned. The code is in submask_code/.

# Unit Hydrograph (UH) File

The [Unit Hydrograph File](UH.md) contains a simple unit hydrograph describing the distribution of travel times of surface runoff over the land surface to the nearest channel. There is no specific information describing the creation of this file. However, an example unit hydrograph file is available in the sample dataset for the Stehekin Basin. This unit hydrograph can be used elsewhere with appropriate modifications. Information on the unit hydrograph may be found in [Lohmann, et al. (1996, 1998)](../References.md).

# Elevation Mask File

When processing meteorological data an elevation maskfile is needed. This is a file that have the same resolution that your final model will have. Each pixel have the mean gridcell elevation in it. The gridcells that are not in the basin have a defined VOID number in them. This is usually a zero (0) or the number -9999.

The elevation maskfile is produced the exact same way as the fraction file but instead of using zeros and ones, use real elevations and zeros. Use the real boundaries cover to clip the DEM. Calculate the zonalmean using the cube_resampled1 and the true_dem1\. Resample to your resolution and gridascii the file.

# Plotting and Control of Flow Network

## Generic Mapping Tools (GMT)

There is a GMT script available to plot the routing network.

MO_GMT.SCR is a GMT script that uses a fortran program called concoord.f. Compile this as: f77 -o concoord concoord.f. Then use: mo_gmt.scr flow_dir1 RESOLUTION (0.125) A postscript file with the routing network is then produced. The program is not perfect and some manual corrections must now be done. This can be done by use of maps, and or plot the rivers and do a rough estimate of the flowdirection. It shoudn't matter that much because its mostly the cells on the boundaries that have errors. These cells mostly point in the direction towards the center of the basin.


## ArcInfo and C - programs

There are other available programs which can be used to plot your routing network. These programs are written by BVM 1999 and are ment to be a substitution to the programs by Nijssen and O'Donnell. The programs are written in C and works together with ArcInfo.

## Gridnet.c

This program reads the direction file that is produced by the flowgen.c program. It then generates a textfile that can be imported into ArcInfo. Use GENERATE, INPUT filename, LINES The input file to this program have to have only integers in it. Void numbers have to be zeros (0), an example is given below.

```
0 0 0 5
3 3 3 5
1 1 3 5
4 3 3 3
```

Usage: gridnet filename rows cols resolution > newfile


## Sources of stream networks

There are digital coverages of the rivers for the US. An example is GCIP-Rivers, which can be downloaded from:

```
Hydrologic Units of the Conterminous US
GEWEX Continental-scale International Project GCIP
GCIP Reference Data Set (GREDS)
by: Alan Rea & Joel R Cederstrand
US.Geological Survey
Open-File -Report 94 388
```

Also check the HYDRO1k dataset where you downloaded the DEM  

## Preprocessing Scripts and Programs

There is a suite of programs (such as the flowgen.c program mentioned above) that have been written to facilitate the preparation of flow direction and fraction files for routing VIC output. These programs and scripts, as well as the README.route_prep file documenting a streamlined approach (specifically designed for using a 30 arc-second digital elevation file and aggregating output for a 1/8 degree VIC application) are available in the file [rout_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz), also available under the "Preparing Routing Model Inputs" section of the [download page](../../SourceCode/Code.md).
