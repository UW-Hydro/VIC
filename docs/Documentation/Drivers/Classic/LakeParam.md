# VIC Lake/Wetland Parameter File

If desired, VIC can simulate lakes and seasonal inundation of wetlands. These lakes and wetlands are fed by runoff from within the current grid cell, and optionally (as of VIC release 4.1.2) by channel inputs from other grid cells upstream. Simulation of lakes and wetlands can be enabled by setting the **LAKES** option in the [global parameter file](GlobalParam.md) to the name of a lake parameter file.

The lake/wetland model consists of a landcover tile in which some (time-varying) portion of the land surface is covered by water. The portion of the tile covered by standing water is called the _**lake**_, and the portion of the tile where soil and vegetation are exposed is called the _**wetland**_. Thus, the lake can represent either a "permanent" lake or seasonal flooding of vegetated land. The size of the lake as a function of time depends on the hydrologic budget and on the underlying lake bathymetry (or topography of the land being flooded), which is supplied by the user. The vegetation covering the exposed wetland portion of the tile must be defined in the [vegetation parameter file](VegParam.md); i.e. if you wish to simulate a lake/wetland in a given grid cell, you must either choose an existing landcover tile to associate with the lake/wetland, or define a new tile to contain the lake/wetland. The index number of the landcover tile that contains the lake/wetland must be indicated in the lake parameter file (see the table below).

NOTE: when running VIC over multiple grid cells (which is true of almost all applications), inevitably some grid cells will not contain lakes and the user will not want to run the lake model in those cells. **To denote a cell that does not contain a lake, you should set the landcover tile index to -1 in the lake parameter file**.

By default, the lake/wetland tile will receive some fraction of the runoff from the other landcover tiles within the grid cell as inputs; this fraction is the parameter _**rpercent**_. However, larger lakes often require water inputs from a larger area than just the current grid cell in order to balance evaporation and outflow. This is especially true in arid regions. To allow for this, we created a new forcing input variable, **CHANNEL_IN**, that can be included in the meteorologic forcing file or as a separate forcing file (VIC can read two separate forcing files per grid cell). Please see the documentation on [forcing files](ForcingData.md) for more information.

The lake/wetland parameters mentioned above, plus all other necessary lake parameters, are described in the table below.

NOTE: all depths are measured from the water/air surface to the lake bottom, at the lake's deepest point (approximately the geographic center of the lake). This is similar to using the word "height" to refer to the elevation of a mountain's peak.

## Lake/Wetland Parameter File Format

| Variable Name     | Units     | Description                                                                                                                                                                                                                                                   |
|---------------    |---------- |-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| gridcel           | N/A       | Grid cell number                                                                                                                                                                                                                                              |
| lake_idx          | N/A       | Index of the veg tile containing the lake/wetland, with respect to the list of veg tiles given in the [veg param file](VegParam.md) for the current grid cell. Numbering starts at 0 for the grid cell's first veg tile, up to (Nveg-1) for the grid cell's final veg tile. <br><br>For example, if the lake/wetland is contained in the 2nd landcover tile in the list of landcover tiles for this grid cell, then lake_idx = 1 (=2-1). <br><br>**NOTE:*** **To denote a grid cell that does not contain any lakes or seasonal flooding, set lake_idx to -1.** <br><br>*NOTE*: this may require changing your veg param file (if it already exists) or preparing it in a special way. |
| numnod            | N/A       | Number of lake layers (or "nodes")                                                                                                                                                                                                                            |
| min_depth         | m         | Lake depth below which channel outflow is 0.                                                                                                                                                                                                                  |
| wfrac             | N/A       | Width of lake outlet, as a fraction of the lake perimeter                                                                                                                                                                                                     |
| depth_in          | m         | Initial lake depth                                                                                                                                                                                                                                            |
| rpercent          | fraction  | Fraction of grid cell runoff that enters lake (instead of going directly to channel network)                                                                                                                                                                  |

If LAKE_PROFILE is FALSE in the global parameter file, then the lake depth-area profile will be computed assuming that the lake basin has a circular cross-section at all depths, with a radius that depends on depth to the power BETA, which is set in LAKE.h. In this case, the 2nd line of parameters for the grid cell are as follows:

| Variable Name  | Units     | Description                                                                      |
|-------------- |---------- |---------------------------------------------------------------------------------- |
| basin_depth   | m         | Maximum allowable depth of lake                                                   |
| basin_area    | fraction  | Fraction of the grid cell covered by lake when its depth is at its maximum value  |

If LAKE_PROFILE is TRUE in the global parameter file, then the 2nd line of parameters for the grid cell consists of _numnod_ (depth, area) pairs describing the profile of the lake basin, starting from the maximum depth and area and descending to the minimum depth and area:

| Variable Name     | Units     | Description                                                           |
|---------------    |---------- |-------------------------------------------------------------------    |
| basin_depth       | m         | Depth of lake                                                         |
| basin_area        | fraction  | Fraction of the grid cell covered by lake when its depth is depth     |

## Example Lake/Wetland Parameter File:

    ###  In this example, there are 10 lake nodes, the lake/wetland occurs in the
    ###  1st veg tile, and LAKE_PROFILE is TRUE.

    7 0 10 1 0.005 10 0.5
    15 0.12 13.5 0.11 12 0.10 10.5 0.09 9 0.08 7.5 0.06 6 0.04 4.5 0.02 3 0.01 1.5 0.0075
