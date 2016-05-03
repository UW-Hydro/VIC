#!/usr/bin/env python 
''' VIC Image Driver testing '''

import xarray as xr
import numpy as np 
import numpy.testing as npt 


# In[28]:

file = '/Users/diana/workdir/VIC_tests_20160331/examples/Example-Image-Stehekin-base-case/results/Stehekin.history.nc'
domainfile = '/Users/diana/Dropbox/UW/VIC_sample_data/image/Stehekin/parameters/domain.stehekin.20151028.nc'


# In[29]:

ds = xr.open_dataset(file)
ds_domain = xr.open_dataset(domainfile)


# In[56]:

def assert_nan_equal(ds_domain, ds_output): 
    """ 
    test to see if nans in two data arrays are the same in order to check domain 
    first dataarray is domain file, second dataarray is a variable in the output
    dataset for image driver 
    """
    
    # extract mask from domain file 
    mask = ds_domain['mask']
    
    # check that lats and lons are the same 
    try: 
        npt.assert_array_equal(ds.lat, ds_domain.lat) and npt.assert_array_equal(ds.lon, ds_domain.lon)
    except AssertionError: 
        print("Domains for domain file and output files are different")
    
    for da in ds_output.data_vars: 
    
        # slice first time step of output DataArray 
        da_sel = ds_output[da].isel(time=0)

        if (len(ds[da].dims) > 3): 
            # check each layer 
            for layer in range(0, len(ds[da].nlayer)):
                try:
                    npt.assert_array_equal(np.isnan(da_sel.isel(nlayer=layer).values), np.isnan(mask.values))
                except AssertionError:
                    print("NaN values in output files and domain files do not match up in %s layer %1f" % (da, layer))
                    
        else: 
            # check that nans are occurring in the same place in the two arrays 
            try:
                npt.assert_array_equal(np.isnan(da_sel.values), np.isnan(mask.values))
            except AssertionError:
                print("NaN values in output file and domain file do not match up in %s" % da)


# In[59]:

# calling function 
assert_nan_equal(ds_domain, ds)

