# VIC Model Overview

# Main Features

The VIC model ([Liang et al., 1994](../Documentation/References.shtml#VIC)) is a large-scale, semi-distributed hydrologic model.  As such, it shares several basic features with the other land surface models (LSMs) that are commonly coupled to global circulation models (GCMs):

1. The land surface is modelled as a grid of large (>1km), flat, uniform cells
    - Sub-grid heterogeneity (e.g. elevation, land cover) is handled via statistical distributions
1. Inputs are time series of daily or sub-daily meteorological drivers (e.g. precipitation, air temperature, wind speed)
1. Land-atmosphere fluxes, and the water and energy balances at the land surface, are simulated at a daily or sub-daily time step
1. Water can only enter a grid cell via the atmosphere
    - Non-channel flow between grid cells is ignored
        -  The portions of surface and subsurface runoff that reach the local channel network within a grid cell are assumed to be >> the portions that cross grid cell boundaries into neighboring cells
    - Once water reaches the channel network, it is assumed to stay in the channel (it cannot flow back into the soil)

This last point has several consequences for VIC model implementation:

1. Grid cells are simulated independently of each other
    - Entire simulation is run for each grid cell separately, 1 grid cell at a time, rather than, for each time step, looping over all grid cells
    - Meteorological input data for each grid cell (for the entire simulation period) are read from a file specific to that grid cell
    - Time series of output variables for each grid cell (for the entire simulation period) are stored in files specific to that grid cell
1. Routing of stream flow is performed separately from the land surface simulation, using a separate model (typically the routing model of Lohmann et al., 1996 and 1998)

## Land Cover and Soil

![VIC Grid Cell Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_grid_cell_schematic_thumb.jpg)

### Land Cover

- can subdivide each grid cell's land cover into arbitrary number of "tiles", each corresponding to the fraction of the cell covered by that particular land cover (e.g. coniferous evergreen forest, grassland, etc.)
- geographic locations or configurations of land cover types are not considered; VIC lumps all patches of same cover type into 1 tile
- includes a lake/wetland cover type
- fluxes and storages from the tiles are averaged together (weighted by area fraction) to give grid-cell average for writing to output files
- for a given tile, jarvis-style veg stomatal response used in computing transpiration
- considers canopy energy balance separately from ground surface
- accounts for partial veg cover fraction, allowing for bare soil evaporation from between the individual plants
- supports optional input of daily timeseries of LAI, albedo, and partial veg cover fraction from forcing files instead of using the monthly climatology specified in the veg library or veg parameter files

#### Note Regarding Evapotranspiration and Time Step

        In order to compensate for the inaccuracies in simulating canopy interception and evaporation at a 24-hour time step, VIC makes an exception for the 24-hour case: in this case, canopy evaporation is allowed to encompass not only the water in the canopy at the beginning of the time step, but also any precipitation, up to the atmospheric demand for water.  At smaller time steps, canopy evaporation is limited to just the amount of water stored in the canopy at the beginning of the time step.  This can result in a) inaccurate apportioning of total ET between canopy evaporation and transpiration, and b) different behavior between VIC simulations at 24 hour time steps and simulations at smaller time steps (with the biggest differences occurring between 12-hour and 24-hour time steps).  For more information, see Haddeland et al (2006a).

## Soil

- arbitrary number of soil layers, but typically 3
- infiltration into the top-most layers controlled by variable infiltration capacity (VIC) parameterization
- top-most layers can lose moisture to evapotranspiration
- gravity-driven flow from upper layers to lower layers
- ARNO baseflow formulation for drainage from bottom layer
- each land cover tile has its own soil temperature profile
- simulation of frozen soil (Cherkauer and Lettenmaier, 1999)
- Simulates soil temperature spatial heterogeneity within a land cover tile (Cherkauer et al., 2003)
- simulates permafrost-specific processes such as melting of excess ground ice (Adam and Lettenmaier, 2008)
- computes the water table depth as a function of soil moisture and soil texture, as described in Bohn et al. (2013b).
- accounts for the thermal properties of organic soil, as described in Farouki (1981).

## Snow Model

VIC considers snow in several forms: ground snow pack, snow in the vegetation canopy, and snow on top of lake ice.  Main features:

VIC considers snow in several forms: ground snow pack, snow in the vegetation canopy, and snow on top of lake ice. Main features:

- ground snow pack is quasi 2-layer; the topmost portion of the pack is considered separately for solving energy balance at pack surface (Andreadis et al., 2009)
- considers partial snow coverage
- considers blowing snow sublimation (Bowling et al, 2004)

For more information about the snow pack formulation, [click here](SnowModelText.shtml).

![VIC Snow Model Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_snow_model_schematic_thumb.jpg)

# Meteorology

## Meteorological Input Data

VIC can use any combination of daily or sub-daily meteorogolical forcings, from point observations, gridded observations, or reanalysis fields.  At minimum, VIC requires daily {precipitation, max/min air temperature, and wind speed}.  VIC will derive all other needed forcings via the approach described in [Bohn et al., 2013a](../Documentation/References.shtml#VIC), which includes:


- If incoming shortwave radiation or humidity are not supplied as forcings, VIC can estimate their daily average values via the algorithms of Kimball et al. (1997), Thornton and Running (1999), and Thornton et al. (2000). These algorithms are part of a stand-alone system called MTCLIM, produced by Steve Running's Numerical Terradynamics Simulation Group at U. Montana. VIC 4.1.2 contains MTCLIM version 4.3; previous versions of VIC contained MTCLIM 4.2.
- If incoming longwave radiation is not supplied, VIC can estimate this via the Tennessee Valley Authority algorithm (TVA, 1972) or the Prata (1996) algorithm (before VIC 4.1.2, TVA was the default; for VIC 4.1.2 onwards, Prata is the default)
- VIC can disaggregate daily forcings to sub-daily as needed, using a cubic spline to interpolate between min and max temperatures, and deriving the other variables from that (Bohn et al., 2013a)

## Non-Meteorological Input Data

Can read daily timeseries of land cover information such as albedo, LAI, and partial vegetation cover fraction as forcing variables.

## Elevation Bands

VIC can consider spatial heterogeneity in precipitation, arising from either storm fronts/local convection or topographic heterogeneity.  Here we consider the influence of topography, via elevation bands.  This is primarily used to produce more accurate estimates of mountain snow pack.  This functionality is controlled by the SNOW_BAND option in the [global parameter file](../Documentation/GlobalParam.shtml).  Main features:


- Can subdivide the grid cell into arbitrary number of elevation bands, to account for variation of topography within cell
- Within each band, meteorologic forcings are lapsed from grid cell average elevation to band's elevation
- Geographic locations or configurations of elevation bands are not considered; VIC lumps all areas of same elevation range into 1 band
- Fluxes and storages from the bands are averaged together (weighted by area fraction) to give grid-cell average for writing to output files
- However, the band-specific values of some variables can be written separately in the output files

For more information about the snow/elevation band formulation, [click here](SnowBandsText.shtml).

![VIC Snow Bands Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_snow_bands_schematic_thumb.jpg)
](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_snow_bands_schematic.gif)

## Frozen Soil Formulation

### Soil Thermal Solution

VIC can use either the approximate soil temperature profile of [Liang et al. (1999)](../Documentation/References.shtml#VIC) or a finite difference solution that takes soil ice content into account, vis a vis [Cherkauer and Lettenmaier (1999)](../Documentation/References.shtml#VIC).

- Liang et al. (1999): set QUICK_FLUX to TRUE in global parameter file; this is the default for FULL_ENERGY = TRUE and FROZEN_SOIL = FALSE.
- Cherkauer et al. (1999): set QUICK_FLUX to FALSE in global parameter file; this is the default for FROZEN_SOIL = TRUE.
    - By default, the finite difference formulation is an explicit method.
    - By default, the nodes of the finite difference formulation are spaced linearly.

For more information about the frozen soil formulation, [click here](FrozenSoilText.shtml).

![VIC Frozen Soil Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_frozen_soil_schematic_thumb.gif)

### Permafrost

These apply to the case QUICK_FLUX = FALSE and FROZEN_SOIL = TRUE, i.e. the formulation of [Cherkauer and Lettenmaier (1999)](../Documentation/References.shtml#VIC).

- global parameter file option: IMPLICIT: uses an implicit scheme to solve the soil thermal profile. This is the default scheme, as of release 4.1.2.
- global parameter file option: EXP_TRANS: uses exponential node spacing (dense node spacing near soil surface; sparse node spacing at depth)

![VIC Permafrost Enhancements Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_permafrost_enhancements_thumb.jpg)

### Temperature Heterogeneity: "Spatial Frost"

- Details described in Cherkauer et al., 2003
- Linear (uniform) distribution of soil temperature around a mean
- Allows some moisture movement in soiil when the average temperature is below freezing

![VIC Spatial Frost Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_spatial_frost_thumb.jpg)

## Dynamic Lake/Wetland Model

The lake/wetland model ([Bowling and Lettenmaier, 2010](../Documentation/References.shtml#VIC)) handles the impoundment of surface water within a grid cell.  Each grid cell is allowed to have a lake/wetland system contained within one of its landcover tiles.  Here, a **_lake_** refers to any impounded surface water, including permanent lakes and seasonal flooding of vegetated land.  The lake's area can vary with time as a function of storage and topography (bathymetry).  In this context, **_wetland_** refers to the exposed portion of the landcover tile.

### Lake Model


- Multi-layer lake model, based on the model of Hostetler and Bartlein (1990), Hostetler (1991), and Hostetler et al (2000)
- Energy balance model
- Mixing, radiation attenuation, variable ice cover
- Dynamic lake area as a function of storage
    - This relationship must be specified as an input parameter
    - For the case of a grid cell containing only one lake, the storage-area relationship == lake bathymetry (plus nearby topography)
    - For multiple lakes within one grid cell, the model considers a single, composite lake; storage-area relationship in this case != lake bathymetry (algorithms for this are under development)
- As of release 4.1.2, lakes can be linked directly to channel network. Lakes can receive inflows from both a) runoff from the surrounding upland within the same grid cell and b) channel flows from upstream grid cells (see Gao et al., 2011 for an example application).
- Lakes drain directly into the channel network. Lake outflows consist of:
    - channel flow (modeled as flow over a broad-crested weir)
        subsurface flow (which will flow into the wetland if the wetland is dry)
- To run the lake model, user must set LAKES to the name of a suitable lake parameter file in the global parameter file. Cells that do not contain lakes can be denoted within the lake parameter file.
- To turn the lake model off completely, the user must either set LAKES to FALSE or omit any mention of LAKES in the global parameter file.

![VIC Dynamic Lake Model Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_dynamic_lake_model_thumb.jpg)

### Wetland Model


- Lakes and wetlands exist in their own "tile" within the grid cell
- Dynamic wetland area = tile area - dynamic lake area
- Allows seasonal inundation of wetlands as lake grows and shrinks
- Wetland moisture/energy flux computations are similar to those of upland tiles
- Wetland soils will tend to be wetter than upland soils due to frequent inundation and recharge by lake

![VIC Dynamic Lake/Wetland Model Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_dynamic_lake_wetland_model_thumb.jpg)

# Carbon Cycle Processes

- VIC optionally simulates photosynthesis, autotrophic respiration, and heterotrophic respiration, as described in Bohn et al. (2013b)
- Plant phenology is NOT dynamic. Vegetative biomass is not simulated. LAI and other vegetation characteristics are prescribed in the same way as when carbon cycling is turned off.
- There are three soil carbon reservoirs: litter (residence time of 2.86 years), intermediate (residence time of 33.3 years), and slow (residence time of 1000 years).
- The previous year's total net primary productivity (NPP; photosynthesis minus autotrophic respiration) is added to the litter pool during the current year, with a constant flux equal to 1/365 * (previous year's total NPP).
- Respiration from the three soil carbon pools is proportional to (the amount of carbon stored in the pool) * exp(-residence time) * (Lloyd-Taylor temperature dependence) * (function of soil moisture).
- A constant fraction of respiration from each pool enters the atmosphere as CO2. For the litter pool, the remainder of respired carbon is sent to the intermediate and slow pools. For the intermediate pool, the remainder is sent to the slow pool. For the slow pool, all respired carbon is sent to the atmosphere.

# Flow Routing

- Routing of stream flow is performed separately from the land surface simulation, using a separate model, typically the routing model of Lohmann, et al. (1996; 1998)
- Each grid cell is represented by a node in the channel network
- The total runoff and baseflow from each grid cell is first convolved with a unit hydrograph representing the distribution of travel times of water from its points of origin to the channel network
- Then, each grid cell's input into the channel network is routed through the channel using linearized St. Venant's equations

![VIC Routing Model Schematic Link](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/images/VIC_river_routing_model_thumb.jpg)
