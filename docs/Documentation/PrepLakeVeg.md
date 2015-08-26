# Advice on Preparing the VIC Veg Parameter File When Simulating Lakes

As of VIC 4.1.1, if you want to simulate lakes, then for each grid cell, you MUST designate a veg tile to be the lake/wetland tile. This has some important implications for the veg parameter file. Here we describe two possible strategies for preparing the veg parameter file to accommodate lakes.

## Background

VIC's lake model allows lakes to change area as their volume changes. This is not required; you can define a lake's depth-area relationship so that lake area is constant for all depths. However, if you do desire your lake area to change with lake storage, you will need to consider: what type of vegetation will the lake cover (via flooding) if it expands in area? And, what type of vegetation will be exposed if the lake shrinks?

At the moment, VIC is not complex enough to handle dynamic vegetation. In other words, the lake is viewed as a temporary landcover ON TOP of some CONSTANT vegetation cover. When the lake expands, the vegetation that gets covered by the lake stops functioning. Similarly, as the lake shrinks, any newly-exposed vegetation immediately begins functioning, as if it were mature vegetation.

For this reason, users may decide that it is bad to have a forest tile contain a dynamic lake. But what if the lakes in your region of interest are surrounded by forest? Or a mix of vegetation cover? Which tile do we put the lake in?

Additionally, if you are trying to use an existing VIC veg param file that was created for an earlier version of VIC (without lakes), how can you adapt this file to work with the addition of lakes?

We have considered two strategies for dealing with lake-veg coordination:

## Strategy 1: Add lakes to an existing tile

One advantage of this strategy is that, if you have a veg parameter file already, from an earlier version of VIC (pre-lake), you may not need to change the veg parameter file. Presumably, this veg parameter file has been created in such a way as to ensure that, for each grid cell, the areas of all of the veg tiles add up to approximately 1.0\. This may have been accomplished by assigning portions of the area covered by open water to the other vegetation tiles (thereby increasing the areas of these other tiles).

If at least one existing veg tile is large enough to contain the lake area, then you may designate this tile to be the lake/wetland tile. This means that you do NOT need to modify your existing veg param file. However, keep in mind that this will result in lakes covering part of that veg tile. If this causes the resulting landcover to be VERY different from the actual proportions observed in the landscape, you may want to try strategy 2.

If you do not already have an existing veg param file, then you can accomplish strategy 1 by simply modifying the program you use to create the veg param file to:

1.  Compute the total areas of each of the land cover classes in the grid cell, including a new "open water" class
2.  Find the largest of the non-open-water ("dry") classes.
3.  Add the open water area to the area of this largest "dry" class. The tile containing this class will now be your lake/wetland tile.
4.  When writing the veg param file, only create tiles for the "dry" classes. The largest "dry" tile now contains the area of the lakes.
5.  When writing your lake param file, use the index of this largest tile as the lake_idx.

## Strategy 2: Add a new lake/wetland tile to the grid cell

If you do not already have an existing veg param file, you can decide in advance (before you create the veg param file) what veg class will be "next to" the lakes, i.e. what type of vegetation will be flooded/exposed as the lakes grow/shrink.

NOTE: it is OK to define a veg tile that has the same veg class as another veg tile within the same grid cell.

You can then create your veg param file via the following steps (in a veg-param-file-creating program):

1.  Compute the total areas of each of the land cover classes in the grid cell, including the "open water" class
2.  Create tiles for each of these classes, but for the "open water" tile, assign to it the veg class that you earlier decided to associate with lakes.
3.  When writing your lake param file, use the index of this tile as the lake_idx.

If you already have an existing veg param file, this strategy unfortunately requires changing it. One way to do this is to:

1.  Rescale the areas of the existing tiles in the grid cell: multiply them by (1.0-Awater), where Awater = fractional area of open water; this will create space for the lakes
2.  Create a new tile, with area = Awater; veg class = the class you decided earlier (this can be the same class as in some other tile in the grid cell)
3.  When writing your lake param file, use the index of this tile as the lake_idx.
