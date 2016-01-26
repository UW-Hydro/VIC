# VIC Model Parameters

!!! Warning
    Docs not written yet.

# Converting from VIC 4 to VIC 5

We are developing a few scripts to aid users in moving from VIC 4 to VIC 5.  These scripts are currently distributed via [`tonic`](https://github.com/UW-Hydro/tonic).

### Installing Tonic Using Anaconda

The easiest way to install Tonic and its dependencies is to use the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda/).

To install Anaconda, follow these two simple steps (check to make sure the installer version is the most current)

1.  download and run the Anaconda installer:  [http://continuum.io/downloads](http://continuum.io/downloads)

2.  setup a virtual environment for Tonic
```shell
conda create -n tonic anaconda
source activate tonic
```

*Note:  you'll need to do the `source activate tonic` to activate the Tonic virtual environment in any new shells.*

Now, download the Tonic source code:

```shell
git clone git@github.com:UW-Hydro/tonic.git
```

From the Tonic source code repository, Tonic can be installed using Python's `distutils`:

```shell
python setup.py install
```

This installs a top level script, `tonic`, into your bin/ directory and the `tonic` package into your Python path.

If you don't want to use the Anaconda installation I've shown above, you can build the package in your local python installation using:
```python
python setup.py develop
```

### Option 2a:  Using a local Python Install (With Write Permissions)

If you have write permissions to the location of your Python distribution, you can just run

```shell
python setup.py install
```

from the top level Tonic directory.  This will install Tonic into your `$PYTHONPATH`.

### Option 2b:  Using a local Python Install (Without Write Permissions)

If you do not have write permissions, you can install Tonic in your local `$PYTHONPATH` by following these steps:

Create a `lib/python` directory in your `$HOME` directory:

```shell
mkdir -p $HOME/lib/python/
```

Add this library path to your `$PYTHONPATH` in your `.bashrc`:

```shell
export PYTHONPATH=$HOME/lib/python:$PYTHONPATH
```

Run `setup.py`:

```shell
python setup.py install --home=$HOME
```

### Testing your install

From the following command line:

```shell
vic_utils -h
python -c 'import tonic'
```

If you don't get any errors, you should be ready to use `tonic`.

## Using Tonic to convert VIC 4 parameters to netCDF format for VIC 5

### Option 1: Using the `vic_utils` command-line utility

Tonic includes a command-line utility for converting VIC style parameters to gridded netCDF. Usage of the utility can be found by running:

```shell
vic_utils grid_params --help
```

### Option 2: Using the `tonic` api

This option allows the user to specify the format of the individual VIC parameter files.  This can be quite useful if the format does not match the assumptions used by the command-line utility.

Example usage:

```Python
from tonic.models.vic.grid_params import soil, snow, veg, veg_class, Cols, Desc

n_veg_classes = 4
root_zones = 3
months_per_year = 12

# Read the soil parameters
soil_dict = soil('~/workdir/Stehekin_VIC.4.1.2_soil.txt', c=Cols(nlayers=3))

# Read the snow parameters
soil_dict = soil('~/workdir/Stehekin_VIC.4.1.2_snow.txt', c=Cols(snow_bands=5))

# Read the veg parameter file
veg_dict = veg('~/workdir/Stehekin_VIC.4.1.2_vegparam.txt',
               soil_dict,
               lai_index=False,
               veg_classes=n_veg_classes)

# Read the veg library file
veg_lib = veg_class('~/workdir/Stehekin_VIC.4.1.2_veglib.txt',
                    skiprows=1)

# Determine the grid shape
target_grid, target_attrs = calc_grid(soil_dict['lats'], soil_dict['lons'])

# Grid all the parameters
grid_dict = grid_params(soil_dict, target_grid, version=version
                        veg_dict=veg_dict, veglib_dict=veg_lib, snow_dict=snow_dict)

# Write a netCDF file with all the parameters
write_netcdf('~/workdir/example.params.vic5.nc', target_attrs,
             target_grid=target_grid,
             soil_grid=grid_dict['soil_dict'],
             snow_grid=grid_dict['snow_dict'],
             veglib_dict=veg_lib,
             veg_grid=grid_dict['veg_dict'],
             version='5.0.dev')
```
