# Preparation of the Soil Parameter File

The soil parameter file is used by VIC to describe the unique soil properties (in addition to several other variables) for each grid cell in the model domain. VIC 4.2 supports an Ascii (column) Format - this has one file with a column for each soil parameter and a seperate row for each grid cell.

The complete list of parameters, as well as sources of data to estimate them and what ranges of reasonable values would be, are available at the [soil parameter file structure page](SoilParam.md).

For the LDAS domain (-125 to -67 degrees longitude; 25 to 53 degrees latitude) at 1/8 degree resolution, a useful source of data is at the [LDAS Parameters web site](http://ldas.gsfc.nasa.gov/nldas/NLDASsoils.php). Raw soil texture data for the continental US can be obtained from [Penn State's Earth System Science Center](http://www.soilinfo.psu.edu/). Certain parameters can be estimated from a soil-texture index ([see summary table](soiltext.shtml)).

* * *

## Procedure for Generating Soil Parameter Files for "Global" Runs

This procedure was taken directly from the [tutorial page](http://www.hydro.washington.edu/~nathalie/VIC_FILES/VIC_SoilVeg_processing.html) of Nathalie Voisin.

The soil data are taken from the 5 minute FAO soil map of the world. This data was downloaded from the web as part of a prototype MS windows program to query this data, aggregate it, etc.

The windows based program works in the following way:

1.  Generate a pedon attribute (`*.pop`) and map unit file for the variable and horizons of interest: soilprogram -> Generate -> Map Unit Data

    You can specify the soil layers in which you are interested, as well as the variable of interest (e.g. porosity, wilting point, etc).

    The program then samples the pedon data base, which contains information about these properties for a large number of soil types. This results in a file that contains the appropriate soil information for each map unit (`*.map`), as well as a file with the statistics for the various soil pedons (`*.pop`).

2.  Generate a data surface file (`*.sfc`): soilprogram -> Generate -> Data Surface



    After selecting the appropriate image file (in my case world.img) a window will pop up (after pressing OK) which will ask you to specify the required area further.

    You can also aggregate further, but I wrote my own aggregation program (Scalesrf.c), because the soil program was very slow. However, if your area is small maybe this will work OK.

    Assuming that you will use Scalesrf, leave the resolution at 5'.

    Two files are created, a `*.srf` file with the data and a `*.doc` file which describes the `*.srf` file.

3.  Scale the `*.srf` file to a coarser resolution (assuming that you did not yet do this in the soilprogram) using Scalesrf. This program assumes that the input resolution is 5' and that the output resolution is a multiple of 5'.

    Usage:

    Scalesrf<infile><outfile>

    <output resolution="">infile and outfile are the input and output files WITHOUT the .srf and .doc extensions.</output>

    </outfile></infile>

    Check the program to make sure that no adjustments are needed for your particular application.

4.  The resulting files can be displayed by the soilprogram: soilprogram -> View Map -> Data Surface

For the global application I used a top layer of 30 cm and a second layer of 70 cm (for a total depth of 100 cm). All the files in the SURFACES directory are based on these assumptions, with the filenames containing the string '30' indicating the top layer and '100' indicating the bottom layer. The string `2` indicates the files aggregated to 2x2 degrees, the other files are the 5' files before aggregation.

In the first iteration I used the van Genuchten parameters provided by the soilprogram as soil hydraulic parameters. This requires editing of the VIC code to actually use this data. Because of some problems (which in hindsight may not have been associated with these parameters), I later switched to texture-based soil hydraulic parameters, using a lookup table based on Cosby et al. [1984] (Cosby, B. J., G. M. Hornberger, R. B. Clapp, and T. R. Ginn, A statistical exploration of the relationships of soil moisture characteristics to the physical properties of soils, Wat. Resour. Res., 20(6), 682-690, 1984). The files associated with this are identified by the string "corby" (because I cannot spell).

The lookup table is based on soil texture classes. To determine the soil texture class for each of the grid cells, use the usda.triangle program in the soilprogram directory (you may need to recompile using the source in source.usda.triangle.zip.gz).

Once you have the soil texture type for each of the soil type, you can use the script texture2prop.awk in the FAO Soil CD to convert to actual soil hydraulic properties. This script uses the FAO Soil CD, which contains the actual lookup table.

So basically, all you need from the soilprogram are:

*   %sand
*   %silt %clay
*   bulk density,

The rest comes from the lookup table. Some of the soil thermal properties are available from the soilprogram as well.

* * *

## Preparing the Soil Parameter File

The steps described here required a Windows machine for the SoilProgram, and a UNIX/LINUX/ETX machine in order to run the C-shell scripts. You also need awk, c, and f90\. This page is certainly and unfortunately not complete. Some digging is expected (windows to dos and vice versa, code changes for your own use, format changes, etc.).

You need:

*   Soilprogram in Global Soil Data Products CD-ROM (IGBP-DIS)
*   A few other scripts and programs: `usda.triangle,texture2prop.awk,(van_GN.txt), Scalesrf.c, cosby.1984.csv, aggr_soil_2vic.c, (get_van_GN.c), make_arc30_ll_data.scr (make_5min_ll_data.scr), run_awk.scr (run_awk_nv.scr), clipSoilVariable.c, clipSoilVariableClass.c, create_soilf.c, make_clipsoil.scr`.
*   Read about the Soil Parameter File structure [here](SoilParam.md).

## 1.1 The Soil Program

Use the soil program to retrieve the sand, clay, bulk density, field capacity, Ksat, WaterN, WiltPoint, and ThetaS as described in the VIC documentation. Before using the program, if you do not want to go through this process several times, think about the following:

*   How many soil layers do you want? ( 2 or 3 in general)
*   Think about the corresponding depths. Usually layer 1 is 10cm, layer 2 and 3 will be calibrated but give an estimate, say 30cm and 150cm for example.
*   What are the boundaries of the mask file that you will be using. For example in the Zambeze basin, xmin=18, xmax=36, ymin=-21, ymax=-9.
*   Know the resolution, say 0.0833 degrees

The SoilProgram runs on a Windows machine. The executable is in `D:\progam:soildata`, choose "soildatasystem" for the starting option.

*   Create map unit files ( .pop and .map)
    *   Plan to do it for each soil layer. (I originally wanted to do layers deeper than 150cm but too much data was missing deeper than 150 cm.)

        soilprogram -> Data -> Generate Map Unit Data

    *   Create both boxes: -> Create Pedon Attribute file AND -> Create map unit file ( if relatively small basin (not for the globe) Select (.pop). AND (.map).
    *   In the box to the bottom right, enter the depth of the soil ( 0 -10cm, then 10-30cm, then 30cm-150cm).In the scroll down, select the variables (sand, clay, bulk density, field capacity, Ksat), some are in the secondary variables.
    *   Note that Ks(cm/day) is the saturated water content, not the saturated conductivity (CondSat). In my area of interest, there were too many missing data to create the CondSat map unit files.
*   Create the surface data files (.srf) AND one `*.doc` file
    *   (optional step; do it for small basins else go to next bullet).
    *   Extract an image file for the basin of interest. ( go to tools -> extract image, then choose Maps ->YourContinent.img, choose output image name (basin.img), enter. Then choose the box to be clipped ( (-9,18) and (-21,36), resolution is 5 minutes, so that I did not need to change it).
    *   NOTE that F. Su had a bigger domain than I had and she chose not to change the resolution here but to aggregate the parameters later, after the Soil Program step. If you need the entire continent or world, you do not need to do this step.

    *   Then create data surface using the extracted image; go to Data -> create data surface; image file is the .img you just created or YourContinent.img, map unit file are all the `*.map` files that you created in the previous step. This will create the `*.srf` files. Do this step as many times as you have `*.map` files. Use the Soil Data option.
    *   Then repeat the same step (create a datasurface) for ONLY ONE variable only (say Clay10.map) using the IDRISI option this time, creating another clay10.img file and most importantly clay10.doc (ascii) file that describes the `*.srf` files and will be used afterwards.

The resulting files can be displayed by the soilprogram: soilprogram -> View Map -> Data Surface.

You are hopefully done with the soil program. The Soil program also offers products from the UMD.

## 1.2 Disaggregate

The programs in this section are found in the file [spatial_disagg.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/spatial_disagg.tgz), also available under the "Post-processing" section of the [download page](../SourceCode/Code.md). Note that each .srf file needed to be treated for the undesired carriage return. Either copy each file to .txt, then vi , then :1,%s/.$//, then save :wq OR here is a script of mine that will do it automatically for all your .srf files: convert_srf_to_text.scr.

IMPORTANT: The following scripts and codes are ONLY VALID for a ending spatial resolution that is a multiple of 5 minute ( 1/8th of a degree for example, 0.5 degree, etc). If you are aiming at say a 0.3 degree resolution, use the SoilProgram aggregation step. Do the following if you DID NOT change the resolution in the soilprogram before. See F. Su large domain case.

*   `make_arc30_ll_data.scr` is Fengge's script do dissagregate `*.srf` file from 5 minutes to 30 arc seconds.
*   Edit the script in order to make it for each variable and depth layer : foreach (BulkDens10 BulkDens30 BulkDens150 Sand10 Sand30 etc) with the names you gave to your `*.srf` files. You also need to give to right `*.doc` file created as explained above. Finally, edit for the number of rows and columns you expect to have in you "box".
*   output: lonlat30_* (Sand/Clay/BulkDens)
*   `run_soil_aggr.scr` will aggregate 30 arc seconds to 1/8 degree.
*   I did not need this one so you still may need to update it to your basin boundaries, hard coded names and cell size
     /**manually change NaN to some values**/ and need to download and compile `aggr_soil_2vic.c`
*   output: `vic_plata_BuldDens SandClayBuld.txt`
*   Note that Laura used `Scalesrf.c` to perform this step ( see her readme file).Your call at this point.

Do the following if your `*.srf` files have the desired resolution (you ask the soil program to put it to the right resolution).

*   `make_my_ll_data.scr` will create for each variable and depth a file with the lat lon and the corresponding variable value.
*   In the script, you do not need to change manually (anymore) to which resolution you want it. Read the script for more details, but as an example, usage is: make_my_ll_data.scr 1/4 outputfiles_soilprgm/Clay1.doc outputfiles_soilprgm/ box_latlon</dl

## 1.3 Make SandClay.txt

Based on the amount of clay and sand, the soil texture will be chosen later ( texture based soil file ...). More on this later.

*   Here is what to do:
*   `_tail +2 Sand10.txt >! junk_`
*   `_tail +2 Clay10.txt >! junk2_`
*   `_paste junk junk2 > sandclay10.txt_`
*   Repeat for each soil layer

## 1.4 Assign texture to your soil based on % sand and % clay

You need fortan90 for this one so either you have it or go to plane2.

The `usda.triangle` program calls a function that classifies soil in the USDA textural triangle using sand and clay %. (Created by: aris gerakis)

*   Download the program.
*   Unzip and recompile: `f90 -o usda.triangle INPOLY.F90 TRIANGLE.F90 WHAT_TEX.F90`
*   usda.triangle
*   input: sandclay.txt (no "NaN")

*   output: by default this is _soilclas.out_ . Change name before doing it for each soil layer.

## 1.5 Assign hydraulic parameter to your soil texture

Create the hydraulic parameters according the "cosby.1984.csv" look up table. You need run_awk_nv.scr,texture2prop.awk, cosby.1984.csv.

*   Edit `run_awk_nv.scr` to make sure that the path is correct ( and do it foreach layer).
*   Output is `para_hydro.txt` ( change name for each layer afterwards)
*   Column: Ks, Porosity, Field Cap, WiltPoint, b

## 1.6 Retrieve lonlat_PT.txt

You need to get for each cell the annual average precipitation and if you plan on using the energy balance mode, the annual mean temperature as well, alse can be kep at -9999\. At this point, I just used the 0.5 degree global soil file I had in hands and use the _regrid program_ (Symap) to downscale it to my desired resolution.

*   `awk '{ if ($1==1) printf(" %.2f " , $26)}' soilfile >! T_tobegridded.txt`
*   `awk '{ if ($1==1) printf(" %.2f " , $49)}' soilfile >! P_tobegridded.txt`
*   `awk '{ if ($1==1) printf(" %.2f %.2f\n",$3,$4 )} soilfile >! station_list` ( then add number of row in first line_
*   Then regrid using the station_list, T_tobegridded.txt then P_tobegridded.txt, the mask file of your new basin.
*   Paste `<lat><lon><precip><temp>on` `lonlat_PT.txt.</temp></precip></lon></lat>`

Those guidelines are specific to our lab ( regrid prgm, global soilfile etc). Since you should have those data from your own VIC input files, you should derive them yourself, your way.

## 1.7 Finish up the soil file

Retrieve last parameters and put it in the right format:

1.  At this point, paste all the soil layer file for each parameters together:
    *   paste KSlayer1 KSlayer2 Kslayer3 >! Ks.txt, etc
    *   As an example see my script make_parametersfiles.scr.
    *   One last note: if you do not want to have to put 0 and 1 in the first column later on, you can do it here.
    *   I used the following scripts: clipSoilVariable.c, clipSoilVariableClass.c, and make_clipsoil.scr. Input files are the pasted files ( 3 columns in each parameters file.
    *   As before, some edits have to be done for your personal use ( hardcoded lat_lon_files and number of cells to be read, format of files, etc)...
2.  get the Exponent (expt) parameter n, from the Brooks-Corey relationship using the [van_GN.txt](http://www.hydro.washington.edu/~nathalie/VIC_FILES/van_GN.txt). ( needs editing of the hard coded paths and filenames)
    *   _do appropriate changes in path and hard coded names, number of soil layers if needed (I personally updated it to read the clipped output files from usda.trianle (soilclass.out) and the porosity file that is clipped as well and has the three soil layers in one file (step a above).It could be updated to use files from the SoilProgram (Theta in particular ...). NOTE that the porosity file must be in fraction units, not percent, whereas sand and clay are expected in percent, else edit ..._
    *   `gcc -lm -o get_van_GN get_van_GN.c`
    *   `Usage: get_van_GN sandclay1.txt sandclay2.txt sandclay3.txt latlon.txt van_GN.txt`
    *   `Note: assume soil texture are the same in one column`
3.  The present code (create_soilf.c) has been used last by Laura (I believe) to create the soil file in its final format ( with 2 soil layers only). Some files that are used by the code are internal to the lab (ARNO parameters, Bart's calibrated parameters). This version is then just a template to create your own version based on the VIC documentation for the soil file. My own version (N. Voisin) will be posted soon. This version does not use ARNO parameters and uses only the parameters described in this web page. In any cases:
    *   Are hardcoded in the code :
    *   the path to the (clipped) files
    *   the number of cells in the basin
    *   the units of the soilprogram variable
    *   the depth of the three layers
    *   `gcc -o create_soilf create_soilf.c -lm`
    *   `create_soilf lonlatelev.txt yoursoilfile`

NOTE: the ARNOS PARAMS option in the global file should be set up to FALSE. Read Nijssen 2001 to see how the ARNO params are used in the soil file and the VIC update log to see if this is still a default or not. Else find the calc_arno.scr script.

NOTE: At this point you may wish to apply the following two steps:

*   fix_resid.pl - This does the quality control on your soil file.
*   set_moist.pl - This sets your initial soil moisture based on soicharacteristics.

These are both available in the file [soil_param_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/soil_param_prep.tgz), also available under the "Preparing VIC Model Inputs" section of the [download page](../SourceCode/Code.md).
