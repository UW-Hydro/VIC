# Sample Data Sets

*   Stehekin Basin ([VIC_sample_data](https://github.com/UW-Hydro/VIC_sample_data)): example of meteorological forcings, and VIC parameters for learning to use the VIC classic and image drivers.

# Test Data Sets

The [VIC test suite](../Development/Testing.md) includes a series of datasets routinely used for testing and validating VIC model output. Reference the test suite [documentation](../Development/Testing.md) for more information on these datasets.

# VIC Input and Output Data Sets

Here are some notable input and output data sets for VIC. If you use these data sets, please cite the corresponding references.

!!! note
    Some of the following datasets are formatted for earlier versions of VIC. Users should not expect that they are can be used with the current VIC version without modification. If you do modify one of these datasets for use with the current version of VIC, we would appreciate it if you alerted the VIC developers so we can update the publicly available dataset here.

## Global Datasets

*   [**Global VIC Input Parameters at 0.5-Degree Resolution**](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/vic_global_0.5deg.html)
    *   NOTE: There are some grid cells in this dataset for which total soil depth exceeds the thermal damping depth of 4 meters. This causes VIC 4.1.2 (and later) to crash on those cells. The solution is to either increase the thermal damping depth (dp in the soil parameter file) to be >= total soil depth (= sum of layer depths 1, 2, 3). We are working on fixing this. For most users, this will not be a problem in their basins.
    *   Spatial Resolution: 0.5-degree
    *   Available data:
        *   Soil Parameter File
        *   Vegetation Parameter File
        *   Vegetation Library File
        *   Snowbands File
    *   References: [Nijssen et al. (2001b)](../Documentation/References.md)
*   [**Global Meteorological Forcing Files at 0.5-Degree Resolution**](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/met_global_0.5deg.html)
    *   Spatial Resolution: 0.5-degree
    *   Temporal Resolution: daily
    *   Time period: 1948-2007
    *   Available data:
        *   Precipitation
        *   Max Temperature
        *   Min Temperature
        *   Wind Speed
    *   References: [Adam and Lettenmaier (2003); Adam et al (2006)](../Documentation/References.md)
*   [**Global Simulations at 2-Degree Resolution**](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/vic_global.html)
    *   Spatial Resolution: 2-degree
    *   Temporal Resolution: daily
    *   Time period: 1980-1993
    *   Model Version: VIC 4.0.3
    *   Simulation Mode: daily water balance
    *   Available data:
        *   Daily Meteorological Forcings (Precip, Tmax, Tmin, Wind, Vapor Pressure, Downward Shortwave, Net Longwave)
        *   Daily Hydrologic Budget Terms (EvapoTranspiration, Runoff, Snow Water Equivalent, Soil Moisture Storage, Total Storage)
        *   Maps of Gridded Parameters (Land/sea mask, depth of top 2 soil layers, total storage capacity)
    *   References: [Nijssen et al, 2001a,b,c](../Documentation/References.md)
*   [**Global River Routing Networks at Various Resolutions**](http://secure.ntsg.umt.edu/publications/2011/WKMS11/)
    *   Spatial Resolution: 2-, 1-, 1/2-, 1/4-, 1/8-, and 1/16-degree
    *   Coarser resolutions were upscaled from finer resolutions using DRT algorithm
    *   Available data:
        *   Grids of flow direction and flow accumulation
    *   References: [Wu et al., 2011](http://www.ntsg.umt.edu/node/724)

## Continental US Datasets

*   [**Simulations of Continental US at 1/8-degree resolution**](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/VIC_retrospective/index.html)
    *   Spatial Resolution: 1/8-degree
    *   Time period: 1950-2000
    *   Model Version: VIC 4.0.3
    *   Simulation Mode: 3-hour energy balance
    *   Available data:
        *   Daily Meteorological Forcings (Precip, Tmax, Tmin, Wind)
        *   3-hourly/Daily/Monthly Water and Energy Budget Terms (EvapoTranspiration, Runoff, Snow Water Equivalent, Soil Moisture, Net Shortwave and Longwave Radiation, Latent and Sensible Heat Fluxes, Ground Flux, Surface Temperature)
        *   Maps of Gridded Parameters (Elevation, Soil Properties, Land Cover)
    *   References: [Maurer et al, 2002](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/VIC_retrospective/index.html)
*   [**Gridded Meteorological Forcings Over Continental US at 1/8-degree resolution**](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/gridded/index.html)
    *   Spatial Resolution: 1/8-degree
    *   Fields: Daily Precip, Tmax, Tmin, Wind
    *   Available data:
        *   Colorado River Basin, 1915-2006, using method of [Maurer et al. (2002)](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/gridded/index.html)
        *   Continental US, 1949-2005, using method of [Maurer et al. (2002)](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/gridded/index.html)
        *   Continental US, 1915-2003, using method of [Hamlet and Lettenmaier (2005)](http://www.hydro.washington.edu/SurfaceWaterGroup/Data/gridded/index.html)
        *   Continental US, Near-real-time, 1915-present (1-day lag), as used by the [UW Surface Water Monitor](http://www.hydro.washington.edu/forecast/monitor/), using method of Wood (2008)
