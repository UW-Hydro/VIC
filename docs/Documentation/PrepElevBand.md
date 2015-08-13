h1 id="ElevationBand">Elevation Band File

The program `elevband.c`, found in the file [snowband_prep.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/snowband_prep.tgz) in the Preparing VIC Model Inputs section of the [download page](../SourceCode/Code.md), will prepare the snow elevation band input file based on the ASCII column format soil input file and an ASCII raster elevation file (with Arc/Info header) at a resolution finer than the resolution of the VIC model run. The soil input file is used simply because it provides a list of grid cell IDs along with grid cell latitude and longitude, both of which are needed by `elevband.c`.

For each pixel in the VIC simulation, `elevband.c` sets an analysis window around the area of the DEM which corresponds to the current pixel. The analysis window will contain elevation values z<sub>i</sub>z<sub>n,</sub>, where

![](../img/Image2.gif)

The range of the elevation bands within each analysis window is calculated as follows:

![](../img/Image3.gif)

Where;

![](../img/Image4.gif)

and NUMBANDS is the number of elevation bands to be modeled for each pixel.

Each elevation, z<sub>i</sub>, is then placed into one of NUMBAND bins of width B<sub>width</sub>, as follows:

![](../img/Image5.gif)

The mean elevation and area fraction are then calculated for each bin:

![](../img/Image6.gif)

Where E<sub>j</sub> and A<sub>j</sub> are the mean elevation and area fraction of each elevation band, _j_, and n<sub>j</sub> are the number of observations in bin _j_.

If the mean elevation of two consecutive bands is less than MINDELTA apart, the bands will be merged and a new mean elevation is calculated. The area fraction of the upper band becomes the sum of the two and the area fraction of the lower band becomes zero. If the area fraction is less than MINFRACTION for any band it is merged with the adjacent band.

Precipitation fraction is calculated in each band based on the band elevation relative to the pixel mean elevation.
