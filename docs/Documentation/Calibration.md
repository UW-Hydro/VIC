# Calibration of VIC and Routing Model

While many of the parameters for these models are based on satellite observations or geological surveys, some of them are either so hetergeneous in space that in situ measurements cannot capture the large-scale "effective" values, or are more conceptual (such as soil layer boundaries) and do not correspond to actual physically-observable quantities. For these parameters, we either must make assumptions about their values or calibrate them (find optimal values for them) that minimize the differences between model output and observations.

## Contents:

*   [General Notes on Calibration](#General)
*   [VIC Model Calibration](#VICCal)
*   [Routing Model Calibration](#RoutCal)

# General Notes on Calibration

We can think of the VIC model and routing model as two parts of one larger modeling framework. Both models contain parameters that require either calibration or _a priori_ estimates of their values.

The total set of calibratable parameters from the two models can be quite large, so here at UW we typically make educated guesses about some of them to reduce the number of parameters to calibrate. Routing model parameters usually fall into the category of educated guesses. For some reasonable parameter values, see the section on [Routing Model Calibration](#RoutCal), below.

The most common observation to use for calibration is streamflow. To use streamflow to calibrate VIC model parameters, we must convert VIC's simulated runoff and baseflow into simulated discharge using the routing model (presumably using reasonable estimates for the routing model parameters).

It is also possible, of course, to calibrate the model to match other observations such as _in situ_ measurements of soil moisture, ET or snow depth, or satellite observations of soil moisture, inundation, or snow coverage/depth. For these, obviously the routing model need not be run.

Dividing the record: for most applications it is necessary to ensure that you have an independent record for calibration and validation (or evaluation) of the streamflow hydrograph. i.e. DO NOT calibrate on your entire period of record, save approximately one half of the time series for an independent evaluation of the calibration.

# VIC Model Calibration

Although the VIC model contains many parameters, it is more appropriate to adjust some of these parameters during calibration than others. Often the distinction is based on the degree to which the parameter values can actually be measured or observed. The parameters most often adjusted during calibration of the VIC model include, [b_infilt](./Info/definitions.md), [Ds](./Info/definitions.md), [Ws](./Info/definitions.md), [Dsmax](./Info/definitions.md) and [soil depth](./Info/definitions.md).

## Soil Parameters

For calibration information involving soil parameters, [click here](CalibrateSoil.md).

## Snow Parameters

For calibration information involving snow parameters, [click here](CalibrateSnow.md).

## Other Parameters

For calibration information involving parameters other than soil and snow, [click here](CalibrateOther.md).

Care should always be taken to ensure that the parameters are within physically-realistic ranges. In many cases these ranges are provided by clicking on the variable name in the [parameter file description](./Inputs.md), but they must always be tailored to the specific application area.

* * *

## Calibration Methods

There are a variety of methods that may be used to optimize the parameters in the VIC model. The hydrology lab at UW most commonly uses the MOCOM-UA method.

*   [MOCOM-UA](MOCOM.md)
*   [Other Methods](CalibrateMethodOther.md)

* * *

## Tips

*   Vegetation Classes and Snowbands

    Vegetation classes and number of snowbands strongly controls computational expense. Gridcells with vegetation classes that cover less than 1-2 % of the gridcell can be removed. This will save time when running the model. If you are looking at change in streamflow as a function of change in vegetation, this should not be done.

*   Compiler Optimazation

    The GNU gcc compiler has three levels of optimization: -O1, -O2, -O3\. Use of these optimizations must be done carefully. As long as the results are the same as for the non-optimized code you can use any optimization level you want.

*   Filter on Precipitation

    Another way to save time during calibration is to calibrate only the gridcells that contribute to about 75% of the basin's streamflow. Some sub-basins in the Columbia River had areas with very little precipitation and some areas with high precipitation. We filtered out the dry areas and ran only the wet ones.

*   Aggregation to One Degree

    If you are working at resolutions as low as 1/8 or 1/16 of a gridcell it might help to run the basin at 1 degree. When you have found a parameter set that gives you reasonable results, the parameters can be extracted out to the closest 1/8 degree gridcells. This can get you close to good results very fast. Currently, there is no program available for this.

# Routing Model Calibration

Our advice regarding routing model calibration is less specific than for the VIC model. Most of our applications here have focused on monthly discharge from large basins, which does not require high accuracy in the routing model parameters. For this reason, we generally do not calibrate the routing model. For parameters, such as flow direction and contributing fraction, that can be obtained automatically from a DEM, we have tools to generate these files. For other parameters such as flow velocity, diffusivity, and the grid cell unit hydrograph, we tend to choose physically reasonable values without further calibration. For details on the calibration process, we refer interested users to the papers by [Lohmann, et al. (1996, 1998)](References.md#Routing).

## Parameters

The routing model contains many parameters, but it is more appropriate to adjust some of these parameters during calibration than others. Often the distinction is based on the degree to which the parameter values can actually be measured or observed. The routing model parameters most often adjusted during calibration are velocity, diffusivity, and the unit hydrograph.

[Lohmann et al. (1996)](References.md#Routing) suggest velocity values in the range of 1 to 3 m/s and diffusivity of 200 to 4000 m<sup>2</sup>/s for the Weser basin in Germany.

[Nijssen et al. (1997)](References.md#Routing) quote velocity values of 0.5 to 2.0 m/s for the Columbia basin and 1.0 m/s for the Delaware.

If only monthly mean flows are required then diffusivity of 800 m<sup>2</sup>/s and velocity of 1.5 m/s are deemed acceptable. If daily values are required then the calibration methodology outlined in [Lohmann et al. (1996, 1998a and 1998b)](References.md#Routing) should be followed.

The example routing parameter files in the Stehekin dataset ([vic.sample.stehekin.tgz](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Datasets), available under Sample Data Sets on the [download page](../SourceCode/Download.md#SampleData)) should serve as a good starting point for your routing parameters.
