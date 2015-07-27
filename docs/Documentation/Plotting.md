# Plotting

VIC model output may be plotted with available Splus (or R) or GMT scripts.

* * *

## Plotting VIC Outputs with Splus (or R)

A few generic Splus scripts (which can also be run in R) have been developed to allow the user to quickly plot VIC model output. The complete library of scripts is available in the file [R_plot_scripts.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Plotting_Scripts/R_plot_scripts.tgz), also available under the "Plotting Scripts" section of the [download page](../downloads/Code.md).

The main script, compare.results.s, compares the default output files (including fluxes, snow, and any other default output fiels depending on global parameter settings) from two model runs, with time series plots of the data, time series plots of the differences, and scatter plots. These scripts produce multi-page postscript-format plots (with .ps extension), which can be viewed with programs such as _ghostview_ and _ghostscript_. Postscript files can also be converted to more common graphic formats such as jpg and gif via tools such as _convert_ and _ps2raster_.

Descriptions of the Splus functions are included in the table below.

| Function Name                 | Description                                                                                                                                                                                                                                   |
|-----------------------------  |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| compare.results.s             | This function reads results from two model simulations and produces a postscript file with time series of simulated values, time series of the differences between simulated values, and scatter plots of the simulated values.               |
| func.plot.compare.nottime.s   | This function takes two time series of data and produces two plots: the first is a time series of both sets of data, the second is a time series of their differences (first - second).                                                       |
| func.plot.scatter.compare.s   | This function takes two time series of data and produces a scatter plot of the data values.                                                                                                                                                   |
| func.read.flux.files.s        | This function reads a VIC model output flux file, and returns the time series of each variable as components in a list. e.g. the first 10 evaporation values for data set tmp.data can be indexed as tmp.data$evap[1:10].                     |
| func.read.snow.files.s        | This function reads a VIC model output snow file, and returns the time series of each variable as components in a list. e.g. the first 10 snow depth values for data set tmp.data can be indexed as tmp.data$depth[1:10].                     |
| func.read.fdepth.files.s      | This function reads a VIC model output fdepth file, and returns the time series of each variable as components in a list. e.g. values of the top-most freeze front for the first 10 time steps can be indexed as "tmp.data$fdepth[1][1:10]."  |

* * *

## GMT Plotting Scripts

The [GMT package](http://gmt.soest.hawaii.edu/) of plotting scripts is another handy way to plot VIC outputs. An example plot script, [plot_flow_STEHE.scr](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Plotting_Scripts/plot_flow_STEHE.scr), that plots observed and simulated hydrographs for the [example Stehekin dataset](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Datasets/vic.sample.stehekin.tgz), is given on the [download page](../downloads/Code.md) under the Plotting Scripts heading. It is written in c-shell, and it calls perl and gmt commands. The output is in postscript format (see note above for viewing and converting postscript files).
