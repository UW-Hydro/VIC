# Basin Delineation

The first step in preparing input files for VIC is to create a list of grid cells contained within your basin. This page explains how to delineate a basin and get a list of grid cells at the desired spatial resolution, using ESRI Arc/Info on _**Windows**_. The basin delineation is also used as the basis for [parameter files for the routing model](Routing/PrepRoutingParams.md).

# Find Basin Boundaries

This section describes how to find the basin boundaries. Several information sources are used, i.e. atlas maps, digital scanned basin boundaries and digital elevation models (DEMs).

## Atlas

The first step in basin delineation should always be to locate the river of interest in a good atlas. A recommended atlas for rivers in the United States is:

    Atlas of River Basins of the United States
    U.S. Department of Agriculture
    Soil Conservation Service, WA D.C. 20250
    Second Edition June 1970

It may be obtained from the University of Washington libraries.

The Atlas is used to get a rough feeling for the total area that your river covers. You can also define the region (maxlatitude, minlatitude, maxlongitude, and minlongitude) that covers the entire drainage basin for your river. This box can be used to reduce filesizes by clipping data from large area datasets (like those for the globe, or the entire U.S.) to the smaller region of interest.

**ARC/INFO**: To make a polygon box that covers your basin area use GENERATE in ArcInfo. You can type in the four corners that you need, project the polygon to the right projection and use this to clip the larger datasets. This polygon will be referred to as your _BOX POLYGON_ throughout the rest of this document.

# Delineate a Basin Using HYDRO1k and Arc/Info

Use Arc/Info GIS together with the known gauge location (e.g. from [USGS](http://waterdata.usgs.gov/nwis/discharge)) or predetermined drainage basin boundaries (e.g. from [HYDRO1k](http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30/hydro)).

## Obtain the HYDRO1k dataset

1.  [Navigate to the USGS HYDRO1k page](http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30/hydro)
2.  Select the continent containing your basin (e.g. North America)
3.  Select drainage basins
4.  Download basins dataset as a tar file in Arc/Info export format

## Prepare Dataset for Importing into Arc/Info

On LINUX:

1.  untar: `tar -xvf na_bas.tar`
2.  unzip: `gzip -d na_bas.e00.gz`

On Windows:

*   Use winzip or another utility to extract the contents of the file

## Start Arc/Info

(On Windows) One way this can be done is:

1.  Open a command prompt (start menu -> run -> cmd.exe)
2.  Change directory (cd) to the directory where you have your workspace
3.  (at the prompt) `arc`
4.  (at the prompt) `arctools`
5.  a menu should appear

## Import the Image

1.  Import the image into arcinfo:`Arc:import auto na_bas.e00 na_bas`
2.  You now have a coverage called na_bas that polygons for multi-leveled drainage basins.

## Project the Coverage

Describe this coverage: `Arc: describe na_bas`

You can see that the projection of this basin is Lambert-Azimuthal. Eventually you will want to change the projection of this basin (to geographic in this example). This can be done at this time or after you have picked out your basin. It is easier to determine which basin you are interested in if you change the projection now, but it takes more time to do it this way.

To change the projection, it is easiest to use "arctools":

1.  Arc: arctools
2.  select "Command Tools"
3.  Edit -> Coordinates -> Project Coverage
4.  Input Cover: na_bas (or right click in space and select)
5.  Output Cover: na_bas_geo (name it whatever you want)
6.  PRJ source: click Define..., Hemisphere, Geographic, enter OK
7.  Close out all the windows

## View this Coverage

You can view the coverage either in the module "ArcEdit" (type `ae`) or using arctools:

1.  Arc command: `arctools`
2.  Select "Edit Tools"
3.  File -> Coverage, Open
4.  In right column, select na_bas_geo
5.  In bottom left column, select ARC, hit OK

## Clean the Coverage

Everytime you perform an operation on a coverage, it needs to be "cleaned" in order to re-establish the polygons.

1.  In Edit Tools Menu -> Arctools -> Commands...
2.  ARCEDIT: clean
3.  ARCEDIT: save
4.  Quit arctools and open again, this time plotting the polygons (the polygon option was previously not available)
5.  File -> Coverage, Open
6.  In right column, select na_bas_geo
7.  In bottom left column, select POLYGON, hit OK

## Highlight Basin of Interest

(NOTE: The basins in HYDRO1k are described with up to six digits (or levels) called the Pfafstetter system. For more information see [the README file](http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30/README).)

Basin Numbering:

*   All basins ending in the numbers 2,4,6,and 8 are the four largest tributaries within a basin at that level.
*   All basins ending in the numbers 1,3,5,7, and 9 are the interbasins between the tributaries.
*   All basins ending in 0 are closed basins (no outflow).

Highlighting:

*   LEVEL 1 is the the largest scale, and LEVEL 6 is the finest.
*   In the Edit Tools menu, go to Edit -> Attribute Selection...
*   A window titled Logical Expresssion pops up:
    *   Our objective is to highlight the Columbia River Basin in the Current expression line, enter the following and then hit Apply expression:
    *   LEVEL1 = 2 (This highlights the Mackenzie River Basin)
    *   LEVEL1 = 4 (This highlights the Nelson River Basin)
    *   LEVEL1 = 6 (This highlights the St. Lawrence River Basin)
    *   LEVEL1 = 8 (This highlights the Mississippi River Basin)
*   So, the columbia river is not one of the four largest basins in this dataset.
*   Now, try the "interbasins":
    *   LEVEL1 = 1 (This highlights a large area of Northwestern North America - including the Columbia)
    *   LEVEL2 = 12 (This highlights the Columbia River Basin)
    *   LEVEL2 = 14 (This highlights the Fraser River Basin)
    *   LEVEL2 = 16 (This highlights the Kuskokwim River Basin)
    *   LEVEL2 = 18 (This highlights the Yukon River Basin)
*   Continuing in this manner, smaller basins can be highlighted.
*   Highlight the Columbia basin again:
*   LEVEL2 = 12

## Select Basin of Interest

Once you have the basin highlighted, merge all polygons in basin:

1.  In the Edit Polygon menu: MERGE
2.  In the Feature Selection menu: hit the Switch Selected Set button (two arrows)
3.  In the Edit Polygon menu: DEL
4.  Now clean the basin as before `arctools\command` then `ARCEDIT: clean`
5.  NOTE: After _clean_ how many polygons you have:
    *   "Built 2 polygon(s) 1 of which have newly created labels" is GOOD (the basin and the background) , go the the next step
    *   "Built 3 polygon(s) 1 of which have newly created labels" is NOT good; you have to manually select the basin, then do as previously; delete the rest and then clean until you have 2 polygons remaining
6.  Manually select the basin: click the arrow on the FEATURE SELECTION menu and select the basin. Then type 9 to get out of the selction mode.
    1.  In the Feature Selection menu: hit the Switch Selected Set button (two arrows)
    2.  In the Edit Polygon menu: DEL
    3.  Clean, 2 polygons? if not,
    4.  Click ALL on FEATURE SELECTION
    5.  In the Edit Polygon menu: MERGE
    6.  Clean (now it should work)

Save as something different; ARCEDIT: `save "newname"` (e.g. for the Columbia Basin, you could do `save "colum_dem")`

Quit arctools. You now have a delineation for your basin.

* * *

## Create a Mask File

1.  Go into grid, then use the following commands:
2.  setwindow xmin ymin xmax ymax (look at the arcinfo help)
3.  then setcell your desired value, 0.125, etc
4.  Use the polygrid command to convert the delineation into a mask
    *   `basin_grd = polygrid ( basin_delineation , # , # , # , cell_size )`
5.  Quit arctools
6.  Use gridascii to export it: `gridascii basin_grd basin_mask.asc`

* * *

## Clip Out a (Fine Resolution) DEM for the Basin

For subsequent use in creating parameter files for VIC and the routing model, you will need a DEM covering the basin that you have just delineated. For the routing model, it is important to have a fine resolution (1-2 km) DEM. You can obtain a fine-resolution DEM from any source (most notably GTOPO30; HYDRO1k is based on GTOPO30).

Obtaining a GTOPO30 DEM:

You can download the 1-km GTOPO30 DEM for free from the [USGS](http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30_info).

1.  [Click here to goto the GTOPO30 DEM from USGS](http://eros.usgs.gov/#/Find_Data/Products_and_Data_Available/gtopo30_info)
2.  Download the DEM for the areas you need. The files from the website above are usualy named. w140n40_tar.gz, etc.
3.  Use `gunzip w140n40_tar.gz` to uncompress the files
4.  Use `tar -xvf w140n40_tar` to get the actual datafiles.

Preparing the DEM:

For example, let's assume you have a DEM for North America called "na_dem"; you have delineated the Columbia basin ("colum_del")

1.  Note: GTOPO30/HYDRO1k DEMs are in raster format. Users of Arc/Info or ArcView can display the DEM data directly after simply renaming the file extension from .DEM (or .HGT) to .BIL. However, if a user needs access to the actual elevation values for analysis in Arc/Info the DEM must be converted to an Arc/Info grid with the command _imagegrid_:

    `ARC: imagegrid na_dem.bil na_dem na_dem.clr`

    or, in GRID:

    `GRID: imagegrid na_dem.bil na_dem`

    _imagegrid_ does not support conversion of signed image data, therefore the negative 16-bit DEM values will not be interpreted correctly. After running _imagegrid_, an easy fix can be accomplished using the following formula in GRID:

    `GRID: na_dem_fix = con(na_dem >= 32768, na_dem - 65536, na_dem)`

2.  The converted grid will then have the negative values properly represented, and the statistics of the grid should match those listed in the .stx file. If desired, the -9999 ocean mask values in the grid could then be set to NODATA with the _setnull_ function.
3.  Look in the \*.prj files to see which projection to use. Used this information when you _projectdefine_ in ArcInfo.
4.  If you have more than one grid downloaded from the web, use the _merge_ command in Arc/Info to merge them together.
5.  Clip the DEM using the basin delineation:

    `gridclip na_dem colum_dem cover colum_del`

6.  View the clipped DEM:
    1.  Bring up the graphics window: display 9999 or &stat 9999
    2.  Set the map extent: `mape colum_dem`
    3.  Set the color scheme: `shadeset rainbow`
    4.  Plot the DEM: `grids colum_dem # linear`
7.  Plot the coverage: `arcs colum_del`
8.  Quit grid: `quit`
9.  export your basin DEM to an ascii file: `gridascii colum_dem colum_dem.asc`

# Other Delineation Methods

## HUC-GCIP-CD

1.  There is a CD available that has ArcInfo coverages of all the large basinsin the US. This CD is called:

    Hydrologic Units of the Conterminous US
    GEWEX Continental-scale International Project GCIP
    GCIP Reference Data Set (GREDS)
    by: Alan Rea & Joel R Cederstrand
    US.Geological Survey
    Open-File -Report 94 388

    Copies of this CD set are available to the local working group, otherwise contact the [GEWEX program](http://www.gewex.org). The file xportarc/usbasin.gz on the CD-ROM contains the US basin boundaries. Transfer this file from the CD to the UNIX system and import it into ARC/INFO as f ollows:

    *   UNIX: _gzip -d usbasins.gz_ to uncompress the file
    *   UNIX: _mv usbasins usbasins.e00_ to rename the file
    *   ARC/INFO: _import auto usbasins.e00 basin1_ to import the file
2.  Remember to _CLEAN_ and _DEFINE PROJECTION_ before using the cover in ARC/INFO. The projection for this file is Lambert-Azhimutal, deta ils are given in the documentation on the CD.
3.  After defining the projection, the file should be reprojected to geographic coordinates (in degrees) so that it matches the coordinates used by t he VIC model.
4.  Next clip the file using the _BOX POLYGON_ created previously. This reduces the file size, and helps eliminate basin delineations which are not part of your watershed.
5.  Finally, open the clipped file in ARCEDIT. You can now manually select and delete basin boundaries which do not fall within the watershed of your river. To help in this selection process use the map from the atlas, the DEM (use it as the background for your editing session), or any other map ima ge which can aid you in determining where the watershed boundary lies. It may also be possible to use the basin HUC codes to generate Another useful pr operty is if the basin HUC codes that follow the basin polygons you can filter out the polygons which have a specific coding.

* * *

## Delineating from Digital Elevation Models (DEMs)

Traditionally basins have been delineated from DEMs simply by following the instructions in ARC/INFO. However, more recently it has been established that if maps of the actual river network are available they can be used to "burn" the river network into the DEM. This results in significantly improv ed basin delineations.

To burn the DEM using an ARC line file of the river network do the following:

*   ARC/INFO:_grid01 = LINEGRID_ converts the line file into a grid file. Use the smae projection and resolution as the DEM.
*   ARC/INFO:_grid02 = CON(ISNULL(river_file), 0, -100)_ creates a new grid file where the river channel has values of -100, everything else isset to 0.
*   ARC/INFO: _grid03 = DEM + grid02_ creates a new grid where the river network has been burned in, so DEM values in the channel are now decreased by 100 (units).

Make sure that the value you use to burn the rivers is big enough to force the proper drainage channels, otherwise watershed may force the basin todrain the wrong way, or include other basins (this is especially true in flat regions).

Once the rivers have been "burned" into the DEM, follow the ARC/INFO instructions to complete the basin delineation.


# Useful Arc/Info Commands and Tips

*   Smooth out a grid: `outfile = RESAMPLE()` in GRID
*   Transform a delineation into polygon : _outfile = POLYGRID ( delineation , # , # , # , cellsize )_ in GRID
*   For the preparation of VIC routing files, one needs the clipped DEM. The `gridclip` command will unfortunately change the size of the DEM window to the often weird one of the delineation (cover file). The trick is to:
    *   `setwindow` the DEM window (if the DEM window was correct for you, i.e. the number of cols/rows % number of cells to be aggregated for the new resolution (.25, .5 etc) is correct)
    *   Use `polygrid (cover, # , #, #, dem res)`
    *   Use: `clipped_dem = con (cover_grid == 2 , dem)`
*   Transform a grid into polygons : `outfile = GRIDPOLY ( grid , # )`
*   Get a delineation in ready-to-plot GMT format:
    *   In Arc/Info, use the UNGENERATE POLY command. If you encounter bad results, use LINE instead
    *   In LINUX: `awk '{ if (NF!=2) print ">" ; else print $0;}' basin.txt >! basin_ready.txt`
    *   In GMT: `psxy basin_ready.txt -M -R$COORD $PROJ -K -O -W3 -V >> $OUTF`
    *   The -M allows it to close the polygons.
*   You have a GRID and want to know the value : in COMMANDS: `cellvalue gridname *` then select the cell. Type 9 to stop.
*   Manually select a basin (in a basin coverage): in EDIT, FILE/open coverage in arcs or polygon, then click on the arrow. 9 to get out of the mode.
*   When doing the delineation, make sure there are only 2 polygons remaining after merging

*   In GRID:

    *   `mape` defines the map extent
    *   `grids` plots the grid
    *   `mape *` will give a different map extent, and use grids to replot with the new map extent
    *   `arcs` or `polys` will plot the basin delination of the previous grid
    *   `gridclip backgroundmap outfile cover basin_delination`
    *   Quit GRID and then gridascii to export it
    *   use arcview to create nice maps
*   Copy an Arc/Info file from someone else:

    You not only need to copy the directory (if grid), but also the info directory and create a log file. Therefore you need to be in _>arc_ to do it, then `>arc: copy "entire path and name" "newname"`

*   In order to get or sent an Arc/Info or Arc/GIS files to/from somebody, in my experience it is better to use the _Interchange File Format_ that you use by entering `EXPORT OPTION(AUTO, COVER, GRID, etc.), FILENAME OUTFILENAME` or by `IMPORT OPTION(AUTO, COVER, GRID, etc) FILENAME`
*   `lg` states the grid in the directory
*   `lc` states the coverage in the directory
*   `kill name` removes the file (including the entire file directory if grid)
*   `asciigrid` imports ascii format (used in UNIX) to grid (ArcInfo/ArcView)
*   `gridascii` exports grid files into ascii files
*   Create a buffer for a grid: `outgrid = EXPAND (grid, # of cell desired for the buffer, LIST , value1, value2, value3, etc.)`
*   Look at the help for this command for mode detaild. In particular, one can use FILE instead of LIST
*   `CON` in the grid menu: it transforms a value by another one using a conditional statement: `outfile = con (isnull(infile),0,infile)` for example to change -9999 into 0 in a maskfile
*   Note that Arcview reads grid files and coverages. It is easier to transform the projection in Arc/Info than in Arcview, (I found anyway, especially for large files)
*   [http://gis.washington.edu/phurvitz/arc_lms/](http://gis.washington.edu/phurvitz/arc_lms/) has help for GIS Spatial Analyst and aml language
