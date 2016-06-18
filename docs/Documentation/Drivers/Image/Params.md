# VIC Model Parameters

To convert ASCII parameter files that were used in previous versions of VIC, you may follow instructions laid out in [ASCII to NetCDF Parameter Conversion](Ascii_to_NetCDF_params.md).

We have also put together an example [Ipython notebook](../../../../samples/notebooks/example_reformat_vic4_parameters_to_vic5image.ipynb), that provides the steps necessary to do this conversion.

After going through these steps and setting up the tonic virtual environment, to look at the parameter file, you can use 'nco' to show the file contents:

```shell
ncdump -h example.params.vic5.nc
```

To access parameter classes, for example, snow bands:

```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.snow_band
```

# Parameter File Information


The VIC model parameter file for the Image Driver, [Parameters File](Params.md), contains spatially distributed parameters describing the land surface.

# Soil Parameters

Soil parameter file serves three main purposes:

*   Define the cell ID number of each grid cell. This ID number is essentially a database key that links a grid cell to its parameters in the various parameter files.
*   Define the grid cell soil parameters
*   Define initial soil moisture conditions, to be used in the absence of an initial state file.

The soil parameters are supplied to VIC in a NetCDF file, with a separate variable for each soil parameter.

For a list of soil parameters, see the [Classic Driver Soil Parameter File](../Classic/SoilParam.md).

To access individual soil parameters, for example, bulk density,  
```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.bulk_density
```

# Vegetation Parameters

For a list of vegetation file parameters, see the [Classic Driver Vegetation Parameter File](../Classic/VegParam.md).

To access individual vegetation parameters, for example, gridcell number,

```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.gridcell
```

# Vegetation Library Parameters

Vegetation parameters needed for each vegetation type used in the VIC model are listed in the [Classic Driver Vegetation Library File](../Classic/VegLib.md). Parameters are given for different vegetation types. Information includes the number of vegetation types per grid cell, and their fractional coverage.

To access individual vegetation library parameters, for example, the trunk ratio for the third vegetation class,

```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.trunk_ratio.sel(veg_class=2)
```


# Elevation Bands (Optional)

By default, VIC assumes each grid cell is flat. This assumption can lead to inaccuracies in snow pack estimates in mountainous grid cells. For this reason, the option exists to have VIC break each grid cell up into a number of _elevation bands_ (also called _snow bands_) and simulate them separately. Each band's mean elevation is used to lapse the grid cell average temperature, pressure, and precipitation to a more accurate local estimate.

The [Parameters File](Params.md) contains information needed to define the properties of each elevation band used by the snow model. Snow elevation bands are used to improve the model's performance in areas with pronounced topography, especially mountainous regions, where the effects of elevation on snow pack accumulation and ablation might be lost in a large grid cell.

The number of snow elevation bands (_option.SNOW_BAND_) to be used with the model is defined in the [global parameter file](GlobalParam.md). The elevation band information from the [Parameters File](Params.md) is only read if the number of snow elevation bands is greater than 1.

It is not necessary that all grid cells have the same number of elevation bands. _SNOW_BAND_ is simply the maximum number of elevation bands specified anywhere in the domain given by the domain file. For relatively flat grid cells, some of the elevation bands will have _AreaFract_ values of 0\. For these zero-area bands, a value of 0 may be supplied for _elevation_ and _Pfactor_.

To access snow band attributes:

```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.snow_band
```

For a list of variable names, see the [Classic Driver Snowbands File](../Classic/SnowBand.md).

To access individual snow band variables, such as elevation for the 1st snowband:

```Python
import xarray as xr
filename = '~workdir/example.params.vic5.nc'
ds = xr.open_dataset(filename)
ds.elevation.sel(snow_band=0)
```
